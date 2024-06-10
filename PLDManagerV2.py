import subprocess
import re
import time
import os


def run_julia_script(script_path, inputfile, args,  timeout=180, output_file="output.txt", method = "sym"):

    num_processes = []

    try:
        # Construct the command to run Julia script
        command = ["julia", script_path] + inputfile

        # Run the command, redirecting output to a file
        with open(output_file, "w") as output_file_handle:
            main_process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.PIPE, text=True)

        start_time = time.time()

        symTasksDone = False

        symTaskFinishedLoading = False

        while True:
            time.sleep(10)  # Avoid busy watching 

            if symTasksDone != True:

                elapsed_time = time.time() - start_time

                if elapsed_time > float(timeout) and symTaskFinishedLoading:
                    print("Timeout reached. Julia script is taking too long. Restarting...")

                    with open("output.txt", "r") as file:
                        lines = file.readlines()

                        if not lines:
                            print(f"File {file_path} is empty. The Julia process must not have started correctly...")
                            return None, None

                        last_line = lines[-1]

                    # Update parameters as needed
                    match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', last_line)

                    if match == False: #Symbolic couldn't complete the first face
                        codim_start, face_start, next_codim_start, next_face_start = incrementFace(str(args[4]) + ".dat", next_face_start, next_codim_start)
                    else:
                        codim_start, face_start, next_codim_start, next_face_start = read_codim_face_from_file(str(args[4]) + ".dat")

                    if codim_start == None or next_face_start == None:
                        break

                    main_process.terminate()  # Terminate the existing process

                    print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))
                    print("Starting a new process to try numerical method.")
                    print("-----------")

                    args[5] = codim_start
                    args[6] = face_start
                    args[7] = "num"

                    with open("PLDinputs.txt", 'w') as file: 
                        for arg in args:
                            file.write(f"{arg}\n")

                    #Create a new process to try the numeric method in the background
                    with open("numOutput" + str(len(num_processes)+1) + ".txt", "w") as output_file_handle:
                        num_processes.append(subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.PIPE, text=True))

                    #Continue with symbolic calculations 
                    
                    print("Continuing on with symbolic calculation in parallel.")

                    args[5] = next_codim_start
                    args[6] = next_face_start
                    args[7] = "sym"

                    with open("PLDinputs.txt", 'w') as file: 
                        for arg in args:
                            file.write(f"{arg}\n")

                    with open(output_file, "w") as output_file_handle:
                        main_process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.PIPE, text=True)               

                    symTaskFinishedLoading = False

                file = open(output_file, "r")

                lines = file.readlines()

                if len(lines) > 0:

                    # Check the last modification time of the output file
                    last_modification_time = os.path.getmtime(output_file)

                    last_line = lines[-1]

                    if last_line == "Finished.":
                        symTasksDone = True

                    

                file.close()

                if last_modification_time > start_time:
                    symTaskFinishedLoading = True
                    # New output detected, update the start time
                    start_time = last_modification_time           
                else:
                    # No new output after the last check, continue waiting
                    continue

            numTasksDone = True

            for task in num_processes:
                if task.poll() == None:
                    numTasksDone = False
                    
            if symTasksDone and numTasksDone:
                break

    except subprocess.CalledProcessError as e:
        print(f"PLD.jl encountered an error! Error message: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)  # Clean up the output file

    return 0

def incrementFace(file_path, prev_face, prev_codim):

    with open(file_path, "r") as file:
        lines = file.readlines()

        last_line = lines[-1]

        # Update parameters as needed
        match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', last_line)

        next_face = prev_face + 1

        if next_face > match.group(3):
            next_codim = prev_codim - 1
            next_face = 1

    return next_face, next_codim

def read_codim_face_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            # Read all lines from the file
            lines = file.readlines()

            # Check if the file is not empty
            if not lines:
                print(f"File {file_path} is empty.")
                return None, None

            # Get the last line
            last_line = lines[-1]

            # Use regular expression to extract codim and face values
            match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', last_line)

            if match:
                codim = int(match.group(1))
                face = int(match.group(2)) + 1
                next_codim = codim
                next_face = face + 1

                if codim == 0:
                    return None, None, None

                if int(match.group(2)) == int(match.group(3)):
                    codim -= 1
                    face = 1
                    next_codim = codim
                    if next_codim != 0:
                        next_face = 2
                    else:
                        next_face = None
                elif int(match.group(2)) == int(match.group(3)) -1:
                    next_codim -= 1
                    next_face = 1

            
                return codim, face, next_codim, next_face
            else:

                print(f"Unable to extract codim and face from the last line of {file_path}.")
                return None, None

    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return None, None


def main():
    # Specify the path to your Julia script
    julia_script_path = "PLDJob.jl"

    # Initial parameters
    edges =  [[1,2],[2,3],[3,4],[4,1]]
    nodes =  [1,2,3,4]
    #masses = [sym.Symbol("m1"),sym.Symbol("m2"),sym.Symbol("m3"),sym.Symbol("m4"),sym.Symbol("m5"),sym.Symbol("m6")]  
    #internalM = [sym.Symbol("M1"),sym.Symbol("M2"),sym.Symbol("M3"),sym.Symbol("M4"),sym.Symbol("M5"),sym.Symbol("M6")]  
    internal_masses_strings =  ["m1", "m2", "m3", "m4"]
    internal_masses = "[" + ",".join(internal_masses_strings) + "]"
    external_masses_strings =  ["p1", "p2", "p3", "p4"]
    external_masses = "[" + ",".join(external_masses_strings) + "]"
    save_output = 'square'
    codim_start = -1
    face_start = 1
    method = "sym"

    args = [edges, nodes, internal_masses, external_masses, save_output, codim_start, face_start, method]

    #write arguments to file
    with open("PLDinputs.txt", 'w') as file: 
        for arg in args:
            file.write(f"{arg}\n")
    
    run_julia_script(julia_script_path, ["PLDinputs.txt"], args)

    print("Sucessfully executed PLD.jl for the provided diagram(s)")

if __name__ == "__main__":
    main()