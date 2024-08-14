import subprocess
import os
import PLDUtils
import sys

def main(mem_limit=0, proc_num="main"):

    if type(mem_limit) == str:
        #Check for units
        if "k" in mem_limit[-3:] or "K" in mem_limit[-3:]:
            mem_limit = int(mem_limit[:-1]) * 1024
        elif "m" in mem_limit[-3:] or "M" in mem_limit[-3:]:
            mem_limit = int(mem_limit[:-1]) * 1048576
        elif "g" in mem_limit[-3:] or "G" in mem_limit[-3:]:
            mem_limit = int(mem_limit[:-1]) * 1073741824

    # Specify the path to your Julia script
    julia_script_path = "PLDJob.jl"

    # Initial parameters
    edges =  [[1,2],[2,3],[3,4],[4,5],[1,5]] #formatted like [[a,b],[c,d],...] with a,b,c,d being integer labels for the vertices of the diagram. 
    #MAKE SURE THAT EACH EDGE IS IN ASCENDING ORDER. (That is, [i,j] s.t. i <= j)
    nodes =  [1, 2, 3, 4, 5] #formatted like [1,2,3,...,n] for an n-point diagram 
    internal_masses =  "[0, 0, 0, 0, 0]" #formatted like [m1,m2,...].
    external_masses =  "[0, 0, 0, 0, 0]" #note that all masses label the SQUARED masses

    output_dir = "output/"
    save_output = "pentagon-massless" #give either a file path or a file name (if you want the file to appear in this directory) WITHOUT the file extension

    codim_start = -1 #integer. Make this <0 if you want to do everything
    face_start = 1 #integer. Make this 1 if you want to do everything in and past the starting codim
    method = "sym" #"sym" or "num". DON'T TOUCH THIS. (The whole point of the wrapper is that it will take care of which method is best on its own)
    single_face = False #Set this to True if you only want to find the discriminant associated wit  h just one face.

    subs = "[]" #Set this to "[]" if you do not need to make any specific substitutions
    #Format your substitutions as "[s => 0, t => a]" etc.

    args = [edges, nodes, internal_masses, external_masses, save_output, codim_start, face_start, method, single_face, subs]

    #write arguments to file
    with open(output_dir + "PLDinputs_proc_" + proc_num + ".txt", 'w') as file: 
        for arg in args:
            file.write(f"{arg}\n")

    #Find codims/faces
    get_faces_path = "PLDGetFaces.jl"

    if single_face == False:
        print("Extracting faces and codimensions")
        codim_array, face_array = PLDUtils.get_faces_codims(get_faces_path, [output_dir + "PLDinputs_proc_" + proc_num + ".txt"], output_dir, proc_num=proc_num)
    else:
        codim_array = []
        face_array = []

    print("Starting calculation of singularities")
    
    PLDUtils.run_julia_script(julia_script_path, [output_dir + "PLDinputs_proc_" + proc_num + ".txt"], args, codim_array, face_array, output_file=output_dir + "output_proc_" + proc_num + ".txt", output_dir = output_dir, max_mem = mem_limit, proc_num = proc_num)

    print("Sucessfully executed PLD.jl for the provided diagram(s)")

    print("-----------")

    print("Compiling output...")

    lines = PLDUtils.compile_diagram_data(save_output)

    with open(save_output + ".txt", "w") as file:
        for line in lines:
            file.write(f"{line}")

    print("Output cleaned up, sorted and printed to file: " + save_output + ".txt")

    print("-----------")

    print("Finally, creating an info file with extra output...")

    #Print out final info to file

    args[7] = save_output + ".txt"
    args[8] = save_output + "_info.txt"

    with open(output_dir + "ExtraInputs_proc_" + proc_num + ".txt", "w") as file:
        for arg in args:
            file.write(f"{arg}\n")

    extra_info_process = subprocess.Popen(["julia", "PLDExtraInfo.jl", output_dir + "ExtraInputs_proc_" + proc_num + ".txt"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)
    while True:
        if extra_info_process.poll() == None:
            continue
        else:
            break


    if extra_info_process.poll() == 0:
        os.remove(output_dir + "ExtraInputs_proc_" + proc_num + ".txt")
        print("Extra info printed to file: " + save_output + "_info.txt")
    else:
        print("An error occured when trying to create the extra output.")
        print("If you want to retry, you can run the command 'julia PLDExtraInfo.jl " + output_dir + "ExtraInputs.txt'.")

if __name__ == "__main__":

    if sys.argv[0] == "./PLDManager.py":
        if len(sys.argv) == 1:
            main()
        elif len(sys.argv) == 2:
            main(mem_limit=sys.argv[1])
        elif len(sys.argv) == 3:
            main(mem_limit=sys.argv[1], proc_num=sys.argv[2])
        else:
            print("WARNING: You may have passed too many arguments to the program!")
            main(mem_limit=sys.argv[1], proc_num=sys.argv[2])
    else: #Command was python3 PLDManager.py
        if len(sys.argv) == 2:
            main()
        elif len(sys.argv) == 3:
            main(mem_limit=sys.argv[2])
        elif len(sys.argv) == 4:
            main(mem_limit=sys.argv[2], proc_num=sys.argv[3])
        else:
            print("WARNING: You may have passed too many arguments to the program!")
            main(mem_limit=sys.argv[2], proc_num=sys.argv[3])
