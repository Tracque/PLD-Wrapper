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
                codim = match.group(1)
                face = match.group(2)

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
                    codim = match.group(1)
                    face = match.group(2)

                    #No need to worry about false relabelling since num files are never sorted
                    all_output.append(["(num) " + line, codim, face])

    sorted_output = sorted(all_output, key=lambda o: (o[1], o[2]))

    #We have extracted all the output now, so can clean up the output files
    os.remove(diagram_name + ".txt")

    for file in num_files:
        os.remove(file)

    final_output =  [output[0] for output in sorted_output]

    return final_output
                
def convert_string_to_array(string):

    elements = string.split(",")

    return [int(e) for e in elements]

if __name__ == "__main__":

    # Initial parameters
    edges =  [[1,2],[2,3],[3,4],[4,1]]
    nodes =  [1,2,3,4] 
    internal_masses =  "[m1, m2, m3, m4]"
    external_masses =  "[p1, p2, p3, p4]"

    #If you needed to use mutltiple calculations (and thus different file names to avoid overwriting) then make the first element of this list the name you want in the end
    output_file_names = ["box", "box-2"] 
    save_output = output_file_names[0]

    args = [edges, nodes, internal_masses, external_masses, save_output + ".txt", save_output + "_info.txt"]

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