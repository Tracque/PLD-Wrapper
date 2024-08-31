import subprocess
import psutil
import re
import time
import os
import copy
import glob


def run_julia_script(script_path, inputfile, args, codims, faces, timeout=90, num_delay = 60, output_file="output/output", output_dir = "output/", max_mem = 0, proc_num = "main"):

    if max_mem == 0:
        max_mem = psutil.virtual_memory().total
    
    program_proc = psutil.Process(os.getpid())

    baseline_mem_usage = psutil.virtual_memory().used

    num_retry_cap = 10
    num_retries = []
    num_cpu_times = []
    num_idle_cycles = []
    inactive_num_processes = 0

    num_processes = []
    num_queue = []
    last_num_start_time = time.time()

    if len(args) == 10:
        save_output = args[4]
        face_start = args[6]
        codim_start = args[5]
    else:
        save_output = args[3]
        face_start = args[5]
        codim_start = args[4]

    #write arguments to file
    with open(output_dir + "PLDinputs_proc_" + proc_num + ".txt", 'w') as file: 
        for arg in args:
            file.write(f"{arg}\n")

    if codim_start < 0 and codims != []:
        codim_start = codims[0]

    try:
        # Construct the command to run Julia script
        command = ["julia", script_path] + inputfile

        # Run the command, redirecting output to a file
        with open(output_file, "w") as output_file_handle:
            main_process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.DEVNULL, text=True)

        baseline_mem_usage += 2147483648 #2GB
        current_estimated_mem_usage = baseline_mem_usage

        start_time = time.time()

        allFacesStarted = False

        symTaskFinishedLoading = False

        while True:
            time.sleep(10)  # Avoid busy watching 

            #Don't want to continue checking and potentially starting new processes when already finished
            if allFacesStarted != True:

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

                            if len(args) == 10:
                                args[5] = codim_start
                                args[6] = face_start
                                args[7] = "num"
                            else:
                                args[4] = codim_start
                                args[5] = face_start
                                args[7] = "num"

                            num_queue.append(copy.copy(args))

                            if codim_start == codims[-1]:
                                allFacesStarted = True
                            else:
                                codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                        else: #Symbolic method computed at least one face

                            #Read what the failed face was
                            codim_start, face_start = read_codim_face_from_file(save_output + ".txt")

                            #Now retry the failed face numerically and increment the face
                            print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))
                            print("Starting a new process to try numerical method.")
                            print("-----------")

                            if len(args) == 10:
                                args[5] = codim_start
                                args[6] = face_start
                                args[7] = "num"
                            else:
                                args[4] = codim_start
                                args[5] = face_start
                                args[7] = "num"

                            num_queue.append(copy.copy(args))
                            
                            if codim_start == codims[-1]:
                                allFacesStarted = True
                            else:
                                codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    else: #No output file exists, we got stuck immediately. (same code here as above)

                        #If couldn't complete first face, then need to retry it numerically and increment the face
                        print("Symbolic method got stuck at codim " + str(codim_start) + " and on face " + str(face_start))
                        print("Starting a new process to try numerical method.")
                        print("-----------")

                        if len(args) == 10:
                            args[5] = codim_start
                            args[6] = face_start
                            args[7] = "num"
                        else:
                            args[4] = codim_start
                            args[5] = face_start
                            args[7] = "num"

                        num_queue.append(copy.copy(args))

                        if codim_start == codims[-1]:
                            allFacesStarted = True
                        else:
                            codim_start, face_start = increment_face(codim_start, face_start, codims, faces)

                    if codim_start == None:
                        break

                    main_process.terminate()  # Terminate the existing process

                    if os.path.exists(output_file):
                        os.remove(output_file)  # Clean up the output file

                    #Don't keep going if we only wanted one face
                    if len(args) == 10 and args[8] == True:
                        allFacesStarted = True
                    elif args[6] == True:
                        allFacesStarted = True

                    if allFacesStarted == False:
                        print("Continuing on with symbolic calculation in parallel.")

                        if len(args) == 10:
                            args[5] = codim_start
                            args[6] = face_start
                            args[7] = "sym"
                        else:
                            args[4] = codim_start
                            args[5] = face_start
                            args[7] = "sym"

                        with open(output_dir + "PLDinputs_proc_" + proc_num + ".txt", 'w') as file: 
                            for arg in args:
                                file.write(f"{arg}\n")

                        with open(output_file, "w") as output_file_handle:
                            main_process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.DEVNULL, text=True)               

                        start_time = time.time()
                        symTaskFinishedLoading = False


                    #Now check if we have the resources to start a numeric process 
                    #(with a timer to ensure that processes have started fully i.e. are close to peak resource usage)
                    #psutil.cpu_percent(interval=1) < 80 
                    current_measured_mem_usage = program_proc.memory_full_info().uss
                    if time.time() - last_num_start_time > num_delay and len(num_queue) > 0 and float(current_estimated_mem_usage/max_mem) < 70 and float(current_measured_mem_usage / max_mem) < 70:

                        #Adjust the inputs to avoid race conditions
                        num_inputs = output_dir + "PLDinputs_proc_" + proc_num + "_" + str(len(num_processes) + 1) + ".txt"
                        if len(args) == 10:
                            num_queue[0][4] = save_output+ "_num_" + str(len(num_processes) + 1)
                            num_queue[0][7] = "num"
                        else:
                            num_queue[0][3] = save_output+ "_num_" + str(len(num_processes) + 1)
                            num_queue[0][7] = "num"

                        with open(num_inputs, 'w') as file: 
                            for arg in num_queue[0]:
                                file.write(f"{arg}\n")

                        #Create a new process to try the numeric method in the background
                        num_processes.append(subprocess.Popen(["julia", script_path] + [num_inputs], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True))
                        num_retries.append(0)
                        num_cpu_times.append(0)
                        num_idle_cycles.append(0)

                        num_queue.pop(0)

                        last_num_start_time = time.time()


                inactive_num_processes = 0



                for i in range(len(num_processes)):
                    if num_processes[i].poll() != None: #Process finished
                        if os.path.exists(save_output + "_num_" + str(i+1) + ".txt"): #Finished with output
                            inactive_num_processes += 1
                        else: #Retry a few times if failed
                            if num_retries[i] < num_retry_cap:
                                print("One of the numeric processes encountered an error! Restarting it...")
                                num_processes[i] = subprocess.Popen(["julia", script_path] + [output_dir + "PLDinputs_proc_" + proc_num + "_" + str(i+1) + ".txt"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)
                                num_retries[i] += 1
                            else:
                                print("Warning: the retry cap of " + str(num_retry_cap) + " has been exceeded.")

                                with open(output_dir + "PLDinputs_proc_" + proc_num + "_" + str(i+1) + ".txt", "r") as file:

                                    lines = file.readlines()

                                    print("The contribution from codim " + lines[5].strip() + ", at face " + lines[6].strip() + " will therefore be missing.")

                                    num_retries[i] += 1

                        
                    else:
                        pid = num_processes[i].pid
                        proc = psutil.Process(pid)
                        cpu_time = proc.cpu_times()
                        if cpu_time.user + cpu_time.system > num_cpu_times[i]:
                            num_cpu_times[i] = cpu_time.user + cpu_time.system
                            num_idle_cycles[i] = 0
                        else:
                            num_idle_cycles[i] += 1

                        if num_idle_cycles[i] > 30: #5 mins idle
                            print("One of the numeric processes has been idling for too long. Restarting it...")
                            num_processes[i].kill()
                            if num_processes[i].stdout:
                                num_processes[i].stdout.close()
                            if num_processes[i].stderr:
                                num_processes[i].stderr.close()
                            num_processes[i] = subprocess.Popen(["julia", script_path] + [output_dir + "PLDinputs_proc_" + proc_num + "_" + str(i+1) + ".txt"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)
                            num_idle_cycles[i] = 0
                            num_retries[i] += 1

                                    
                current_estimated_mem_usage = baseline_mem_usage + (len(num_processes) - inactive_num_processes) * 4294967296 * 2 #8GB per num process

                # Check the last modification time of the output file
                if allFacesStarted == False:
                    last_modification_time = os.path.getmtime(output_file)
                    if main_process.poll() != None: #If we are only waiting on numeric processes, make sure we don't timeout
                        last_modification_time = time.time()

                if main_process.poll() != None and num_queue == []:
                    allFacesStarted = True

                if last_modification_time > start_time:
                    symTaskFinishedLoading = True
                    # New output detected, update the start time
                    start_time = last_modification_time           
                else:
                    # No new output after the last check, continue waiting
                    continue

            numTasksDone = True

            for i in range(len(num_processes)):
                    if num_processes[i].poll() != None: #Process finished
                        if os.path.exists(save_output + "_num_" + str(i+1) + ".txt"): #Finished with output
                            inactive_num_processes += 1
                        else: #Retry a few times if failed
                            if num_retries[i] < num_retry_cap:
                                numTasksDone = False
                                print("One of the numeric processes encountered an error! Restarting it...")
                                num_processes[i] = subprocess.Popen(["julia", script_path] + [output_dir + "PLDinputs_proc_" + proc_num + "_" + str(i+1) + ".txt"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)
                                num_retries[i] += 1
                            else:
                                print("Warning: the retry cap of " + str(num_retry_cap) + " has been exceeded.")

                                with open(output_dir + "PLDinputs_proc_" + proc_num + "_" + str(i+1) + ".txt", "r") as file:

                                    lines = file.readlines()

                                    print("The contribution from codim " + lines[5].strip() + ", at face " + lines[6].strip() + " will therefore be missing.")

                                    num_retries[i] += 1

                        
                    else:
                        numTasksDone = False

                        #Check that numeric processes aren't idling
                        pid = num_processes[i].pid
                        proc = psutil.Process(pid)
                        cpu_time = proc.cpu_times()
                        if cpu_time.user + cpu_time.system > num_cpu_times[i]:
                            num_cpu_times[i] = cpu_time.user + cpu_time.system
                            num_idle_cycles[i] = 0
                        else:
                            num_idle_cycles[i] += 1

                        if num_idle_cycles[i] > 10: #5 mins idle
                            print("One of the numeric processes has been idling for too long. Restarting it...")
                            num_processes[i].kill()
                            if num_processes[i].stdout:
                                num_processes[i].stdout.close()
                            if num_processes[i].stderr:
                                num_processes[i].stderr.close()
                            num_processes[i] = subprocess.Popen(["julia", script_path] + [output_dir + "PLDinputs_proc_" + proc_num + "_" + str(i+1) + ".txt"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)
                            num_idle_cycles[i] = 0
                            num_retries[i] += 1

                    
            if allFacesStarted and numTasksDone:

                os.remove(output_dir + "PLDinputs_proc_" + proc_num + ".txt")

                num_input_files = glob.glob(output_dir + "PLDinputs_proc_" + proc_num + "*.txt")
    
                for file in num_input_files:
                    os.remove(file)
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
    
def get_faces_codims(script_path, input_file_path, output_dir="output/", proc_num = "main"):

    command = ["julia", script_path] + input_file_path

    with open(output_dir + "output_proc_" + proc_num + ".txt", "w") as output_file_handle:
            main_process = subprocess.Popen(command, stdout=output_file_handle, stderr=subprocess.DEVNULL, text=True)

    while True:
        time.sleep(10) #Avoid busy watching

        if main_process.poll() != None:

            with open(output_dir + "output_proc_" + proc_num + ".txt", "r") as file:

                lines = file.readlines()

                if not lines:

                    continue

                else:

                    codim_string, face_string =  lines[-2], lines[-1] #Codim and face respectively
                    #os.remove("output_proc_" + proc_num + "")  # Clean up the output file

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
                
                #Don't falsely relabel an already labeled line
                if line[0] == "(":
                    all_output.append([line, codim, face])
                else:
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
                
                #Don't falsely relabel an already labeled line
                if line[0] == "(":
                    all_output.append([line, codim, face])
                else:
                    all_output.append(["(num) " + line, codim, face])

    sorted_output = sorted(all_output, key=lambda o: (o[1], -o[2]), reverse=True)

    for file in num_files:
        os.remove(file)

    final_output =  [output[0] for output in sorted_output]

    return final_output
                
def convert_string_to_array(string):

    elements = string.split(",")

    return [int(e) for e in elements]


