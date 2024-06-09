import subprocess
import re
import sympy as sym
import time
import os


def run_julia_script(script_path, inputfile, args,  timeout=18000, output_file="output.txt", method = "sym"):
    try:
        # Construct the command to run Julia script
        command = ["julia", script_path] + inputfile

        # Run the command, redirecting output to a file
        with open(output_file, "w") as output_file_handle:
            process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.PIPE, text=True)

        start_time = time.time()
        restarted = False

        while True:
            time.sleep(10)  # Adjust the interval as needed 

            elapsed_time = time.time() - start_time

            if elapsed_time > float(timeout):
                print("Timeout reached. Julia script is taking too long. Restarting...")
                process.terminate()  # Terminate the existing process

                # Update other parameters as needed
                codim_start, face_start = read_codim_face_from_file(str(args[5]) + ".dat", restarted)

                if codim_start == None:
                    break

                print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))

                # Julia script failed, modify parameters and try again
                method = "num" if method == "sym" else "sym"  # Switch method

                print("Retrying this face using numerical method.")
                print("-----------")

                if restarted == False:
                    if len(args) == 6:
                        args.append(codim_start)
                        args.append(face_start)
                        args.append(method)
                    else:
                        args[6] = codim_start
                        args[7] = face_start
                        args[8] = method

                with open("PLDinputs.txt", 'w') as file: 
                    for arg in args:
                        file.write(f"{arg}\n")

                restarted = True

                run_julia_script(script_path, ["PLDinputs.txt"], args, method = method)

                if method == "num":

                    print("Restarting with the symbolic method.")

                    #Update starting point
                    codim_start, face_start = read_codim_face_from_file(str(args[5]) + ".dat", restarted)
                    args[6] = codim_start
                    args[7] = face_start
                    method = "sym"

                    run_julia_script(script_path, ["PLDinputs.txt"], args, method = method)

                break

            # Check the last modification time of the output file
            last_modification_time = os.path.getmtime(output_file)

            file = open(output_file, "r")

            lines = file.readlines()

            if len(lines) > 0:
                if lines[-1] == "Finished.":
                    break

            file.close()

            if last_modification_time > start_time:
                # New output detected, update the start time
                start_time = last_modification_time           
            else:
                # No new output after the last check, continue waiting
                continue
    except subprocess.CalledProcessError as e:
        print(f"PLD.jl encountered an error! Error message: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)  # Clean up the output file

    return 0

def read_codim_face_from_file(file_path, restarted):
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
                    return None, None

                if int(match.group(2)) == int(match.group(3)):
                    codim -= 1
                    face = 1
                return codim, face
            else:
                print(f"Unable to extract codim and face from the last line of {file_path}.")
                return None, None

    except FileNotFoundError:
        print(f"File {file_path} not found.")
        if restarted == True:
            return None, None
        else:
            return 3, 1

def main():
    # Specify the path to your Julia script
    julia_script_path = "PLDJob.jl"

    # Initial parameters
    edges =  [[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,1],[2,5]]
    nodes =  [1,3,4,6,7]
    masses = [sym.Symbol("m1"),sym.Symbol("m2"),sym.Symbol("m3"),sym.Symbol("m4"),sym.Symbol("m5"),sym.Symbol("m6")]  
    internalM = [sym.Symbol("M1"),sym.Symbol("M2"),sym.Symbol("M3"),sym.Symbol("M4"),sym.Symbol("M5"),sym.Symbol("M6")]  
    numberOfMasses = 5
    internal_masses =  [0,0,0,0,0,0,0,0]
    external_masses =  [0,masses[1],0,0,masses[4]]
    save_output = 'FivePoint_2Loop_2mass_squarepent'

    args = [edges, nodes, numberOfMasses, internal_masses, external_masses, save_output, 1, 10, "num"]

    #write arguments to file
    with open("PLDinputs.txt", 'w') as file: 
        for arg in args:
            file.write(f"{arg}\n")
    
    run_julia_script(julia_script_path, ["PLDinputs.txt"], args)

    print("Sucessfully executed PLD.jl for the provided diagram(s)")

if __name__ == "__main__":
    main()
