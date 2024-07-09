import subprocess
import os
import PLDUtils

def main():
    # Specify the path to your Julia script
    julia_script_path = "PLDJob.jl"

    # Initial parameters
    edges =  [[1,2], [1,6], [2,3], [3,7], [4,6], [4,7], [5,6], [5,7]] #formatted like [[a,b],[c,d],...] with a,b,c,d being integer labels for the vertices of the diagram. 
    #MAKE SURE THAT EACH EDGE IS IN ASCENDING ORDER. (That is, [i,j] s.t. i <= j)
    nodes =  [1, 2, 3, 4, 5] #formatted like [1,2,3,...,n] for an n-point diagram 
    internal_masses =  "[0, 0, 0, 0, 0, 0, 0, 0]" #formatted like [m1,m2,...]. See the GUI or PLDJob.jl to see/modify the allowed variable symbols.
    external_masses =  "[0, 0, 0, 0, 0]" #note that all masses label the SQUARED masses

    output_dir = "output/"
    save_output = "output/heptabox-s23-eq-0-lim" #give either a file path or a file name (if you want the file to appear in this directory) WITHOUT the file extension

    codim_start = -1 #integer. Make this <0 if you want to do everything
    face_start = 1 #integer. Make this 1 if you want to do everything in and past the starting codim
    method = "sym" #"sym" or "num". DON'T TOUCH THIS. (The whole point of the wrapper is that it will take care of which method is best on its own)
    single_face = False #Set this to True if you only want to find the discriminant associated with just one face.

    subs = "[s23 => 0]" #Set this to "[]" if you do not need to make any specific substitutions

    if subs != "[]":
        args = [edges, nodes, internal_masses, external_masses, save_output, codim_start, face_start, method, single_face, subs]
    else:
        args = [edges, nodes, internal_masses, external_masses, save_output, codim_start, face_start, method, single_face, subs]

    #write arguments to file
    with open(output_dir + "PLDinputs.txt", 'w') as file: 
        for arg in args:
            file.write(f"{arg}\n")

    #Find codims/faces
    get_faces_path = "PLDGetFaces.jl"

    if single_face == False:
        print("Extracting faces and codimensions")
        codim_array, face_array = PLDUtils.get_faces_codims(get_faces_path, [output_dir + "PLDinputs.txt"], output_dir)
    else:
        codim_array = []
        face_array = []

    print("Starting calculation of singularities")
    
    PLDUtils.run_julia_script(julia_script_path, [output_dir + "PLDinputs.txt"], args, codim_array, face_array, output_dir = output_dir)

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

    with open(output_dir + "ExtraInputs.txt", "w") as file:
        for arg in args:
            file.write(f"{arg}\n")

    with open(output_dir + "output.txt", "w") as output_file_handle:
            extra_info_process = subprocess.Popen(["julia", "PLDExtraInfo.jl", output_dir + "ExtraInputs.txt"], stdout=output_file_handle, stderr=subprocess.PIPE, text=True)

    while True:
        if extra_info_process.poll() == None:
            continue
        else:
            break

    os.remove(output_dir + "output.txt")
    os.remove(output_dir + "ExtraInputs.txt")

    print("Extra info printed to file: " + save_output + "_info.txt")

if __name__ == "__main__":
    main()
