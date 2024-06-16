import subprocess
import psutil
import re
import time
import os
import copy

#TODO: Julia Script to find and output all face/codim numbers DONE
#Thus, fix the stepping over to the next face/codim DONE
#Clean up output files properly (delete them) DONE
#Investigate if the numeric threads are even running properly? STILL TODO
#Queue the extra processes to protect system resources DONE?

def run_julia_script(script_path, inputfile, args, codims, faces, timeout=60, load_timeout=120, output_file="output.txt"):

    num_processes = []
    num_queue = []
    last_num_start_time = time.time()

    save_output = args[4]
    face_start = args[6]
    codim_start = args[5]

    #write arguments to file
    with open("PLDinputs.txt", 'w') as file: 
        for arg in args:
            file.write(f"{arg}\n")

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

            #Don't want to continue checking and potentially starting new processes when already finished
            if symTasksDone != True:

                elapsed_time = time.time() - start_time

                if (elapsed_time > float(timeout) and symTaskFinishedLoading) or elapsed_time > float(load_timeout):
                    print("Timeout reached. Julia script is taking too long. Restarting...")

                    with open(save_output + ".dat", "r") as file:
                        lines = file.readlines()

                        if not lines:
                            print(f"File {save_output + ".dat"} is empty. The Julia process must not have started correctly...")
                            return None, None

                        last_line = lines[-1]

                    # Use regex to match the last line (most recent output)
                    match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', last_line)

                    if (int(match.group(2)) < face_start and int(match.group(1)) == codim_start) or int(match.group(1)) > codim_start: #Symbolic couldn't complete the first face

                        #If couldn't complete first face, then need to retry it numerically and increment the face
                        print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))
                        print("Starting a new process to try numerical method.")
                        print("-----------")

                        args[5] = codim_start
                        args[6] = face_start
                        args[7] = "num"

                        num_queue.append(copy.copy(args))

                        codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    else: #Symbolic method computed at least one face

                        #Read what the failed face was
                        codim_start, face_start = read_codim_face_from_file(save_output + ".dat")

                        #Now retry the failed face numerically and increment the face
                        print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))
                        print("Starting a new process to try numerical method.")
                        print("-----------")

                        args[5] = codim_start
                        args[6] = face_start
                        args[7] = "num"

                        num_queue.append(copy.copy(args))
                        
                        codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    if codim_start == None:
                        break

                    main_process.terminate()  # Terminate the existing process

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

                    #Now check if we have the resources to start a numeric process 
                    #(with a timer to ensure that processes have started fully i.e. are close to peak resource usage)
                    if time.time() - last_num_start_time > 60 and len(num_queue) > 0 and psutil.cpu_percent(interval=1) < 80:
                        #and psutil.virtual_memory().percent < 80

                        #Adjust the inputs to avoid race conditions
                        num_inputs = "PLDinputs" + str(len(num_processes) + 1) + ".txt"
                        num_queue[0][4] = save_output+ "_num_" + str(len(num_processes) + 1)
                        num_queue[0][7] = "num"

                        with open(num_inputs, 'w') as file: 
                            for arg in num_queue[0]:
                                file.write(f"{arg}\n")

                        #Create a new process to try the numeric method in the background
                        with open("numOutput" + str(len(num_processes)+1) + ".txt", "w") as output_file_handle:
                            num_processes.append(subprocess.Popen(["julia", script_path] + [num_inputs], stdout=output_file_handle, stderr=subprocess.PIPE, text=True))

                        num_queue.pop(0)

                        last_num_start_time = time.time()
                            

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

                os.remove("PLDinputs.txt")

                for i in len(num_processes):
                    os.remove("PLDinputs" + str(i+1) + ".txt")
                    os.remove("numOutput" + str(i+1) + ".txt") #Clean up output files

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
        print("Something has gone wrong, the codim is not in the allowed set!")

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

                if not lines:

                    continue

                else:

                    codim_string, face_string =  lines[-2], lines[-1] #Codim and face respectively
                    os.remove("output.txt")  # Clean up the output file

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
