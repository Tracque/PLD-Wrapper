import re
import os
import glob
import subprocess

def cleanup_diagram_data(diagram_names):

    input_files = glob.glob("PLDinputs*.txt")
    
    for file in input_files:
        os.remove(file)

    all_output = []

    #First open the main (symbolic) output file
    for diagram_name in diagram_names:
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

        num_files = glob.glob(diagram_name + "_num_*.txt")

        for file in num_files:

            with open(file, "r") as f:

                line = f.readlines()[0]

                match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', line)

                if match != None:
                    codim = int(match.group(1))
                    face = int(match.group(2))

                    #No need to worry about false relabelling since num files are never sorted
                    all_output.append(["(num) " + line, codim, face])

        #We have extracted all the output now, so can clean up the output files
        os.remove(diagram_name + ".txt")

        for file in num_files:
            os.remove(file)

    sorted_output = sorted(all_output, key=lambda o: (o[1], -o[2]), reverse=True)

    final_output =  [output[0] for output in sorted_output]

    return final_output
                
def convert_string_to_array(string):

    elements = string.split(",")

    return [int(e) for e in elements]

if __name__ == "__main__":

    # Initial parameters
    edges =  [[1,7], [1,6], [2,7], [2,3], [3,6], [4,5], [4,6], [5,7]] #formatted like [[a,b],[c,d],...] with a,b,c,d being integer labels for the vertices of the diagram. 
    #MAKE SURE THAT EACH EDGE IS IN ASCENDING ORDER. (That is, [i,j] s.t. i <= j)
    nodes =  [1, 2, 3, 4, 5] #formatted like [1,2,3,...,n] for an n-point diagram 
    internal_masses =  "[0, 0, 0, 0, 0, 0, 0, 0]" #formatted like [m1,m2,...]. See the GUI or PLDJob.jl to see/modify the allowed variable symbols.
    external_masses =  "[0, 0, 0, 0, 0]"

    #If you needed to use mutltiple calculations (and thus different file names to avoid overwriting) then make the first element of this list the name you want in the end
    output_file_names = ['dpent-s23-0-lim'] 
    save_output = output_file_names[0]

    subs = "[s23 => 0]" #Set this to "[]" if you do not need to make any specific substitutions

    if subs == "[]":
        args = [edges, nodes, internal_masses, external_masses, save_output + ".txt", save_output + "_info.txt"]
    else:
        #A few dummy arguments here so that I can reuse code
        args = [edges, nodes, internal_masses, external_masses, save_output + ".txt", save_output + "_info.txt", "a", "a", "a", subs]

    print("Manually compiling output...")

    lines = cleanup_diagram_data(output_file_names)

    with open(save_output + ".txt", "w") as file:
        for line in lines:
            file.write(f"{line}")

    print("Output cleaned up, sorted and printed to file: " + save_output + ".txt")

    print("-----------")

    print("Finally, creating an info file with extra output...")

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

    print("REMINDER: This output was manually compiled and may not represent the complete Principal Landau Determinant.")