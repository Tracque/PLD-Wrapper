import subprocess
import re
import time
import os

#TODO: Julia Script to find and output all face/codim numbers DONE
#Thus, fix the stepping over to the next face/codim DONE
#Clean up output files properly (delete them)
#Investigate if the numeric threads are even running properly?
#Queue the extra processes to protect system resources

def run_julia_script(script_path, inputfile, args, codims, faces, timeout=60, output_file="output.txt", method = "sym"):

    num_processes = []

    save_output = args[4] + ".dat"
    face_start = args[6]
    codim_start = args[5]

    if codim_start < 0:
        codim_start = codims[0]

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

            #Done want to continue checking and potentially starting new processes when finished
            if symTasksDone != True:

                elapsed_time = time.time() - start_time

                if elapsed_time > float(timeout) and symTaskFinishedLoading:
                    print("Timeout reached. Julia script is taking too long. Restarting...")

                    with open(save_output, "r") as file:
                        lines = file.readlines()

                        if not lines:
                            print(f"File {save_output} is empty. The Julia process must not have started correctly...")
                            return None, None

                        last_line = lines[-1]

                    # Use regex to match the last line (most recent output)
                    match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', last_line)

                    if int(match.group(2)) < face_start or int(match.group(1)) > codim_start: #Symbolic couldn't complete the first face

                        #If couldn't complete first face, then need to retry it numerically and increment the face
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

                        codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    else: #Symbolic method computed at least one face

                        #Read what the failed face was
                        codim_start, face_start = read_codim_face_from_file(str(args[4]) + ".dat")

                        #Now retry the failed face numerically and increment the face
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

                        codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    if codim_start == None:
                        break

                    main_process.terminate()  # Terminate the existing process

                    #Continue with symbolic calculations 

                    if os.path.exists(output_file):
                        os.remove(output_file)  # Clean up the output file
                    
                    print("Continuing on with symbolic calculation in parallel.")

                    args[5] = codim_start
                    args[6] = face_start
                    args[7] = "sym"

                    with open("PLDinputs.txt", 'w') as file: 
                        for arg in args:
                            file.write(f"{arg}\n")

                    with open(output_file, "w") as output_file_handle:
                        main_process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.PIPE, text=True)               

                    start_time = time.time()
                    symTaskFinishedLoading = False

                # Check the last modification time of the output file
                last_modification_time = os.path.getmtime(output_file)

                if main_process.poll() != None:
                    symTasksDone = True

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

def increment_face(codim, face, codim_array, face_array):

    try:
        index = codim_array.index(codim)

        if face == face_array[index]: #if this is the last face of this codim
            codim -= 1
            face = 1
        else:
            face += 1
    except ValueError:
        print("Something has gone wrong, and the codim is not in the allowed set!")

    return codim, face

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

                if codim == 0:
                    return None, None, None

                if int(match.group(2)) == int(match.group(3)):
                    codim -= 1
                    face = 1

                return codim, face
            else:

                print(f"Unable to extract codim and face from the last line of {file_path}.")
                return None, None

    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return None, None
    
def get_faces_codims(script_path, input_file_path):

    command = ["julia", script_path] + input_file_path

    with open("output.txt", "w") as output_file_handle:
            main_process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.PIPE, text=True)

    while True:
        time.sleep(10) #Avoid busy watching

        if main_process.poll() != None:

            with open("output.txt", "r") as file:

                lines = file.readlines()

                os.remove("output.txt")  # Clean up the output file
                codim_string, face_string =  lines[-2], lines[-1] #Codim and face respectively

                return convert_string_to_array(codim_string), convert_string_to_array(face_string)
                
def convert_string_to_array(string):

    elements = string.split(",")

    return [int(e) for e in elements]


def main():
    # Specify the path to your Julia script
    julia_script_path = "PLDJobV2.jl"

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

    #Find codims/faces
    get_faces_path = "PLDGetFaces.jl"
    print("Extracting faces and codimensions")
    codim_array, face_array = get_faces_codims(get_faces_path, ["PLDinputs.txt"])
    print("Starting calculation of singularities")
    
    run_julia_script(julia_script_path, ["PLDinputs.txt"], args, codim_array, face_array)

    print("Sucessfully executed PLD.jl for the provided diagram(s)")

if __name__ == "__main__":
    main()
