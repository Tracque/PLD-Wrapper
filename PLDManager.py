import subprocess
import psutil
import re
import time
import os
import copy
import glob

#TODO: Benchmark the wrapper. This had better be faster!

def run_julia_script(script_path, inputfile, args, codims, faces, timeout=90, output_file="output.txt"):

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

                if (elapsed_time > float(timeout) and symTaskFinishedLoading):
                    print("Timeout reached. Julia script is taking too long. Restarting...")

                    if os.path.exists(save_output + ".txt"):

                        with open(save_output + ".txt", "r") as file:
                            lines = file.readlines()

                            if not lines:
                                print("Output file is empty. The Julia process must not have started correctly...")
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

                            if codim_start == codims[-1]:
                                symTasksDone = True
                            else:
                                codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                        else: #Symbolic method computed at least one face

                            #Read what the failed face was
                            codim_start, face_start = read_codim_face_from_file(save_output + ".txt")

                            #Now retry the failed face numerically and increment the face
                            print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))
                            print("Starting a new process to try numerical method.")
                            print("-----------")

                            args[5] = codim_start
                            args[6] = face_start
                            args[7] = "num"

                            num_queue.append(copy.copy(args))
                            
                            if codim_start == codims[-1]:
                                symTasksDone = True
                            else:
                                codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    else: #No output file exists, we got stuck immediately. (same code here as above)

                        #If couldn't complete first face, then need to retry it numerically and increment the face
                        print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))
                        print("Starting a new process to try numerical method.")
                        print("-----------")

                        args[5] = codim_start
                        args[6] = face_start
                        args[7] = "num"

                        num_queue.append(copy.copy(args))

                        if codim_start == codims[-1]:
                            symTasksDone = True
                        else:
                            codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    if codim_start == None:
                        break

                    main_process.terminate()  # Terminate the existing process

                    if os.path.exists(output_file):
                        os.remove(output_file)  # Clean up the output file

                    #Don't keep going if we only wanted one face
                    if args[8] == True:
                        symTasksDone = True

                    if symTasksDone == False:
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
                    if time.time() - last_num_start_time > num_delay and len(num_queue) > 0 and psutil.cpu_percent(interval=1) < 80 and psutil.virtual_memory().percent < 80:

                        #Adjust the inputs to avoid race conditions
                        num_inputs = "PLDinputs" + str(len(num_processes) + 1) + ".txt"
                        num_queue[0][4] = save_output+ "_num_" + str(len(num_processes) + 1)
                        num_queue[0][7] = "num"

                        with open(num_inputs, 'w') as file: 
                            for arg in num_queue[0]:
                                file.write(f"{arg}\n")

                        #Create a new process to try the numeric method in the background
                        num_processes.append(subprocess.Popen(["julia", script_path] + [num_inputs], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))

                        num_queue.pop(0)

                        last_num_start_time = time.time()

                        for task in num_processes:
                            if task.poll() != None:
                                stdout, stderr = task.communicate() #This should clean up output pipelines
                                        

                # Check the last modification time of the output file
                if symTasksDone == False:
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

                for i in range(len(num_processes)):
                    os.remove("PLDinputs" + str(i+1) + ".txt")
                    #Clean up output files

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
                
def compile_diagram_data(diagram_name):

    all_output = []

    #First open the main (symbolic) output file
    if os.path.exists(diagram_name + ".txt"):
        with open(diagram_name + ".txt", "r") as file:

            lines = file.readlines()

            for line in lines:
                match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', line)
                codim = int(match.group(1))
                face = int(match.group(2))
                all_output.append(["(sym) " + line, codim, face])

        #We have extracted all the output now, so can clean up the output files
        os.remove(diagram_name + ".txt")

    num_files = glob.glob(diagram_name + "_num_*.txt")

    for file in num_files:

        with open(file, "r") as f:

            line = f.readlines()[0]

            match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', line)

            if match != None:
                codim = int(match.group(1))
                face = int(match.group(2))
                all_output.append(["(num) " + line, codim, face])

    sorted_output = sorted(all_output, key=lambda o: (o[1], o[2]), reverse=True)

    for file in num_files:
        os.remove(file)

    final_output =  [output[0] for output in sorted_output]

    return final_output
                
def convert_string_to_array(string):

    elements = string.split(",")

    return [int(e) for e in elements]


def main():
    # Specify the path to your Julia script
    julia_script_path = "PLDJob.jl"

    # Initial parameters
    edges =  [[1, 2], [2, 3], [3, 6], [4, 6], [5, 6], [5, 7], [1, 7], [4, 7]] #formatted like [[a,b],[c,d],...] with a,b,c,d being integer labels for the vertices of the diagram. 
    #MAKE SURE THAT EACH EDGE IS IN ASCENDING ORDER. (That is, [i,j] s.t. i <= j)
    nodes =  [1, 2, 3, 4, 5] #formatted like [a,b,c,d,...] with a,b,c,d being integer labels for the vertices of the diagram
    internal_masses =  "[0, 0, 0, m2, m2, m2, 0, m2]" #formatted like [m1,m2,...]. See the GUI or PLDJob.jl to see/modify the allowed variable symbols.
    external_masses =  "[0, 0, 0, 0, p2]" #note that all masses label the SQUARED masses
    save_output = 'Hj-npl-pentb' #give either a file path or a file name (if you want the file to appear in this directory) WITHOUT the file extension
    codim_start = -1 #integer. Make this <0 if you want to do everything
    face_start = 1 #integer. Make this 1 if you want to do everything in and past the starting codim
    method = "sym" #"sym" or "num". DON'T TOUCH THIS. (The whole point of the wrapper is that it will take care of which method is best on its own)
    single_face = False #Set this to True if you only want to find the discriminant associated with just one face.

    args = [edges, nodes, internal_masses, external_masses, save_output, codim_start, face_start, method, single_face]

    #write arguments to file
    with open("PLDinputs.txt", 'w') as file: 
        for arg in args:
            file.write(f"{arg}\n")

    #Find codims/faces
    get_faces_path = "PLDGetFaces.jl"

    if single_face == False:
        print("Extracting faces and codimensions")
        codim_array, face_array = get_faces_codims(get_faces_path, ["PLDinputs.txt"])
    else:
        codim_array = []
        face_array = []

    print("Starting calculation of singularities")
    
    run_julia_script(julia_script_path, ["PLDinputs.txt"], args, codim_array, face_array)

    print("Sucessfully executed PLD.jl for the provided diagram(s)")

    print("-----------")

    print("Compiling output...")

    lines = compile_diagram_data(save_output)

    with open(save_output + ".txt", "w") as file:
        for line in lines:
            file.write(f"{line}")

    print("Output cleaned up, sorted and printed to file: " + save_output + ".txt")

    print("-----------")

    print("Finally, creating an info file with extra output...")

    #Print out final info to file

    args[4] = save_output + ".txt"
    args[5] = save_output + "_info.txt"

    with open("ExtraInputs.txt", "w") as file:
        for arg in args:
            file.write(f"{arg}\n")

    with open("output.txt", "w") as output_file_handle:
            extra_info_process = subprocess.Popen(["julia", "PLDExtraInfo.jl", "ExtraInputs.txt"], stdout=output_file_handle, stderr=subprocess.PIPE, text=True)

    while True:
        if extra_info_process.poll() == None:
            continue
        else:
            break

    os.remove("output.txt")
    os.remove("ExtraInputs.txt")

    print("Extra info printed to file: " + save_output + "_info.txt")

if __name__ == "__main__":

    num_delay = 60

    main()
