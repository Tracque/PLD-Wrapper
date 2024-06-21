import subprocess
import psutil
import re
import time
import os
import copy
import glob
import webview

#TODO: Implement GUI (keep a non-GUI version)
#TODO: Benchmark the wrapper. This had better be faster!

class interaction_API():
    def __init__(self):
        self.stop = False

    #dummy arg here because the api interaction seems to insist on passing an extra garbage (as far as our purposes) argument 
    def execute_PLD(dummy, edges, nodes, internal_masses, external_masses, save_output, codim_start=-1, face_start=1, single_face=False):

        # Specify the path to the main Julia script
        julia_script_path = "PLDJob.jl"

        # Initial parameters
        codim_start = int(codim_start)
        face_start = int(face_start)
        method = "sym"

        args = [edges, nodes, internal_masses, external_masses, save_output, codim_start, face_start, method, single_face]

        #write arguments to file
        with open("PLDinputs.txt", 'w') as file: 
            for arg in args:
                file.write(f"{arg}\n")

        #Find codims/faces
        get_faces_path = "PLDGetFaces.jl"

        if single_face == False:
            window.evaluate_js('appendToOutput("Extracting faces and codimensions")')
            codim_array, face_array = get_faces_codims(get_faces_path, ["PLDinputs.txt"])
        else:
            codim_array = []
            face_array = []

        window.evaluate_js(r'appendToOutput("Starting calculation of singularities\n")')
        
        run_julia_script(julia_script_path, ["PLDinputs.txt"], args, codim_array, face_array)

        window.evaluate_js('appendToOutput("Sucessfully executed PLD.jl for the provided diagram(s)")')

        window.evaluate_js('appendToOutput("-----------")')

        window.evaluate_js('appendToOutput("Compiling output...")')

        lines = compile_diagram_data(save_output)

        with open(save_output + ".txt", "w") as file:
            for line in lines:
                file.write(f"{line}")

        window.evaluate_js(f'appendToOutput("Output cleaned up, sorted and printed to file: {save_output}.txt")')

        window.evaluate_js('appendToOutput("-----------")')

        window.evaluate_js('appendToOutput("Finally, creating an info file with extra output...")')

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

        window.evaluate_js(f'appendToOutput("Extra info printed to file: {save_output}_info.txt")')

        return 0

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
                    window.evaluate_js('appendToOutput("Timeout reached. Julia script is taking too long. Restarting...")')

                    if os.path.exists(save_output + ".txt"):

                        with open(save_output + ".txt", "r") as file:
                            lines = file.readlines()

                            if not lines:
                                window.evaluate_js('appendToOutput("Output file is empty. The Julia process must not have started correctly...")')
                                return None, None

                            last_line = lines[-1]

                        # Use regex to match the last line (most recent output)
                        match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', last_line)

                        if (int(match.group(2)) < face_start and int(match.group(1)) == codim_start) or int(match.group(1)) > codim_start: #Symbolic couldn't complete the first face

                            #If couldn't complete first face, then need to retry it numerically and increment the face
                            window.evaluate_js(f'appendToOutput("Symbolic method got stuck at codim {codim_start} and on face {face_start}")')
                            window.evaluate_js('appendToOutput("Starting a new process to try numerical method.")')
                            window.evaluate_js('appendToOutput("-----------")')

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
                            window.evaluate_js(f'appendToOutput("Symbolic method got stuck at codim {codim_start} and on face {face_start}")')
                            window.evaluate_js('appendToOutput("Starting a new process to try numerical method.")')
                            window.evaluate_js('appendToOutput("-----------")')

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
                        window.evaluate_js(f'appendToOutput("Symbolic method got stuck at codim {codim_start} and on face {face_start}")')
                        window.evaluate_js('appendToOutput("Starting a new process to try numerical method.")')
                        window.evaluate_js('appendToOutput("-----------")')

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
                        window.evaluate_js('appendToOutput("Continuing on with symbolic calculation in parallel.")')

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
                        num_queue[0][4] = save_output + "_num_" + str(len(num_processes) + 1)
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
                    os.remove("numOutput" + str(i+1) + ".txt") #Clean up output files

                break

    except subprocess.CalledProcessError as e:
        window.evaluate_js(f'appendToOutput("PLD.jl encountered an error! Error message: {e}")')
    except Exception as e:
        window.evaluate_js(f'appendToOutput("An unexpected error occurred: {e}")')

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
        window.evaluate_js('appendToOutput("Something has gone wrong, the codim is not in the allowed set!")')

    return codim, face

def read_codim_face_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            # Read all lines from the file
            lines = file.readlines()

            # Check if the file is not empty
            if not lines:
                window.evaluate_js(f'appendToOutput("File {file_path} is empty.")')
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

                window.evaluate_js(f'appendToOutput("Unable to extract codim and face from the last line of {file_path}.")')
                return None, None

    except FileNotFoundError:
        window.evaluate_js(f'appendToOutput("File {file_path} not found.")')
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
                codim = match.group(1)
                face = match.group(2)
                all_output.append(["(sym) " + line, codim, face])

        #We have extracted all the output now, so can clean up the output files
        os.remove(diagram_name + ".txt")

    num_files = glob.glob(diagram_name + "_num_*.txt")

    for file in num_files:

        with open(file, "r") as f:

            line = f.readlines()[0]

            match = re.search(r'codim: (\d+), face: (\d+)/(\d+)', line)

            if match != None:
                codim = match.group(1)
                face = match.group(2)
                all_output.append(["(num) " + line, codim, face])

    sorted_output = sorted(all_output, key=lambda o: (o[1], o[2]), reverse=True)

    for file in num_files:
        os.remove(file)

    final_output =  [output[0] for output in sorted_output]

    return final_output
                
def convert_string_to_array(string):

    elements = string.split(",")

    return [int(e) for e in elements]



if __name__ == "__main__":
    
    html_content = r"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>PLD User Interface</title>
        <style>
        /* Basic reset */
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: Arial, sans-serif;
            background-color: #f4f4f4;
            display: flex;
            flex-direction: column;
            align-items: center;
            min-height: 100vh;
            margin: 20px;
            line-height: 1.4;
        }

        .container {
        display: flex;
        justify-content: space-between;
        width: 90%;
        max-width: 1600px;
        }

        .left-column {
            display: flex;
            flex-direction: column;
            width: 45%;
        }

        #box1, #box2 {
            width: 100%;
        }

        #box3 {
            width: 45%;
        }

        @media (max-width: 768px) {
            .container {
                flex-direction: column;
                align-items: center;
            }

            .left-column, #box3 {
                width: 100%;
            }
        }

        .box {
            background-color: #fff;
            border: 2px solid #ddd;
            border-radius: 8px;
            padding: 20px;
            flex: 1 1 45%;
            margin: 10px;
            text-align: center;
            transition: all 0.3s ease;
        }


        .box:hover {
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            transform: translateY(-5px);
        }

        #box3 {
            width: 100%;
        }

        #outputArea {
            width: 100%;
            height: 300px; /* Adjust the height as needed */
            background-color: #f5f5f5; /* Light grey background */
            border: 1px solid #ccc; /* Light grey border */
            border-radius: 4px;
            padding: 10px;
            font-family: monospace; /* Monospaced font to mimic console or code block */
            font-size: 14px;
            color: #333; /* Dark text color */
            overflow-y: scroll; /* Vertical scroll */
            resize: none; /* Disable resizing */
        }


        </style>
    </head>
    <body>
        <div class="container">

            <div class="left-column">

                <div class="box" id="box1" style="text-align:left">
                    <h1 style="font-size:30px; text-align:center;">PLD User Interface</h1>
                    <p style="line-height: 1.2">
                        Keep your inputs to the suggested format to ensure compatibility with PLD.jl <br>
                        I have tried to allow for many possible variable labels for the squared masses, including p, P, m, M, q, Q, l and L. (as in m1, m2, m3,... or similar) <br>
                        You should ensure you align the edges/internal masses and the nodes/external masses. For example, this means if your edges are input as [[1,2],...], then the internal mass array [m1,...] will set the mass of the [1,2] edge to m1.
                    </p>
                    <br>
                    <div style="text-align: center;">
                        Edges <input type="text" placeholder="e.g. [[1,2],[2,3],[1,3]]" id="edges"> <br>
                        Internal Squared Masses <input type="text" placeholder="e.g. [m1,0,m1,m2]" id="internals"> <br>
                        Nodes <input type="text" placeholder="e.g. [1,2,3]" id="nodes"> <br>
                        External Squared Masses <input type="text" placeholder="e.g. [p1,p1,p2,0]" id="externals"> <br>
                        <button id="visualiseButton">Generate Diagram Visualisation</button>
                        <br><br>
                    </div>
                    <p style="line-height: 1.2">
                        Note: The provided visualisation is algorithmically generated and is not supposed to produce literature ready diagrams.<br>
                        Hopefully, it will at least give you an idea of if you are inputting what you think you are.
                        The solid lines are massive, the dashed ones massless.
                    </p>
                
                    <br>
            
                    <p style="line-height: 1.2">
                        If you wish to calculate the full Principal Landau Determinant, leave the starting codim and face input fields below blank. (note: if you leave EITHER blank, the program will default to calculating the full PLD)
                    </p>
            
                    <br>
            
                    <div style="float:left">
                    <input type="checkbox" id="startCheckbox" name="runType" value="singleFace">
                    <label for="startCheckbox"> Calculate singularities on a single face?</label><br>
                    Starting codim <input type="text" id="startC"> <br>
                    Starting face <input type="text" id="startF"> <br>
                    </div>

                    <br><br><br><br>
                    Below, input the desired output file name or path and name for the calculation. (DO NOT INCLUDE A FILE EXTENSION. THIS WILL BE DONE FOR YOU)<br>
                    File Path <input type="text" id="filePath"> <br>
                    <button id="startButton">Start Calculation</button>
                </div>

                <div class="box" id="box2">
                    <h2>Output:</h2> <br>
                    <textarea id="outputArea" readonly>Output will appear here...</textarea>
                </div>

            </div>

            <div class="box" id="box3">
                <canvas id="Visualisation" height="800" width="800"></canvas> <br>
            </div>
        </div>


        <script>
            // JavaScript function to send the input value to Julia
            function sendInputToJulia() {
                var edgeInputs = document.getElementById("edges").value;
                var nodeInputs = document.getElementById("nodes").value;
                var internalInputs = document.getElementById("internals").value;
                var externalInputs = document.getElementById("externals").value;
                var outputFilePath = document.getElementById("filePath").value;
                var startCodim = document.getElementById("startC").value;
                var startFace = document.getElementById("startF").value;

                if (startCodim == "" || startFace == "") {
                    pywebview.api.execute_PLD(edgeInputs, nodeInputs, internalInputs, externalInputs, outputFilePath);
                } else if(document.getElementById("startCheckbox").checked == true) {
                    pywebview.api.execute_PLD(edgeInputs, nodeInputs, internalInputs, externalInputs, outputFilePath, startCodim, startFace, true);
                } else {
                    pywebview.api.execute_PLD(edgeInputs, nodeInputs, internalInputs, externalInputs, outputFilePath, startCodim, startFace, false);
                }

                // Show that we have started successfully
                document.getElementById("outputArea").value = "PLD is starting. Output should appear here soon...\n\n";
            }

            function appendToOutput(newLine) {
                var outputArea = document.getElementById("outputArea");
                outputArea.value += newLine + "\n";
                outputArea.scrollTop = outputArea.scrollHeight; // Scroll to the bottom
            }

            function generatePolygonCoords(edges, nodes) {
                // Find the number of points we need
                points = 0;

                for (i=0; i < edges.length; i++) {
                    if (Math.max(...edges[i]) > points) {
                        points = Math.max(...edges[i]);
                    }
                }

                for (i=0; i < points; i++) {
                    coords.push([xCentre + size * Math.cos(i * 2 * Math.PI / points), yCentre + size * Math.sin(i * 2 * Math.PI / points)]);
                    if (nodes.includes(String(i+1))) { 
                        externalLines.push([xCentre + externalSize * Math.cos(i * 2 * Math.PI / points), yCentre + externalSize * Math.sin(i * 2 * Math.PI / points)]);
                    } else {
                        externalLines.push([xCentre + size * Math.cos(i * 2 * Math.PI / points), yCentre + size * Math.sin(i * 2 * Math.PI / points)]);
                    }
                }
            }

            // To allow for visualisation of multiple edges connected to the same two vertices
            function checkDegeneracy(array,index) {
                var degeneracy = 0;
                for (i=0; i < index; i++) {
                    if ((array[i][0] == array[index][0] && array[i][1] == array[index][1]) || (array[i][0] == array[index][1] && array[i][1] == array[index][0])) {
                        degeneracy++;
                    }
                }
                return degeneracy;
            }

            // To make the degenerate edges always point outwards (this is just some extra case handling for [n,1] edges)
            function outwardBezier(edges, i) {
                if ((edges[i][0] == 1 && edges[i][1] == points) || (edges[i][1] == 1 && edges[i][0] == points)) {
                    return [xCentre - (size + bezierSize * checkDegeneracy(edges, i)) * Math.cos((edges[i][0] + edges[i][1] -2) * Math.PI / points), yCentre - (size + bezierSize * checkDegeneracy(edges, i)) * Math.sin((edges[i][0] + edges[i][1] -2) * Math.PI / points)];
                } else {
                    return [xCentre + (size + bezierSize * checkDegeneracy(edges, i)) * Math.cos((edges[i][0] + edges[i][1] -2) * Math.PI / points), yCentre + (size + bezierSize * checkDegeneracy(edges, i)) * Math.sin((edges[i][0] + edges[i][1] -2) * Math.PI / points)];
                }
            }

            function drawDiagram(edges, nodes, internal_masses, external_masses) {
                // Clear the canvas
                ctx.setLineDash([]);
                ctx.clearRect(0, 0, canvas.width, canvas.height);
                ctx.strokeStyle = "black";
                ctx.strokeRect(0, 0, canvas.width, canvas.height);

                // Fetch inputs from the page
                var edgeInputs = document.getElementById("edges").value;
                var nodeInputs = document.getElementById("nodes").value;
                var internalInputs = document.getElementById("internals").value;
                var externalInputs = document.getElementById("externals").value;

                // Parse inputs into arrays
                var edges = JSON.parse(edgeInputs);

                if (!Array.isArray(edges) || !edges.every(item => Array.isArray(item) && item.length === 2)) {
                    throw new Error('Invalid input format for edges');
                }

                var nodes = nodeInputs.replace(/^\[|\]$/g, '').split(',');
                nodes = nodes.map(item => item.trim());

                var internal_masses = internalInputs.replace(/^\[|\]$/g, '').split(',');
                internal_masses = internal_masses.map(item => parseInt(item.trim(), 10));

                var external_masses = externalInputs.replace(/^\[|\]$/g, '').split(',');
                external_masses = external_masses.map(item => item.trim());

                coords = [];
                externalLines = [];

                generatePolygonCoords(edges, nodes);

                // Start drawing the diagram.
                for (i=0; i < edges.length; i++) {
                    // If this if the first edge between these two vertices, just draw a line
                    if (checkDegeneracy(edges, i) == 0) {
                        if (internal_masses[i] == "0") {
                            ctx.beginPath();
                            ctx.setLineDash([10, 15]);
                            ctx.moveTo(coords[edges[i][0] -1][0], coords[edges[i][0] -1][1]);
                            ctx.lineTo(coords[edges[i][1] -1][0], coords[edges[i][1] -1][1]);
                            ctx.stroke();
                        } else {
                            ctx.beginPath();
                            ctx.setLineDash([]);
                            ctx.moveTo(coords[edges[i][0] -1][0], coords[edges[i][0] -1][1]);
                            ctx.lineTo(coords[edges[i][1] -1][0], coords[edges[i][1] -1][1]);
                            ctx.stroke();
                        }
                    // If this is a degenerate edge, draw extras as increasingly wide bezier curves
                    } else {
                        bezierCoords = outwardBezier(edges, i);
                        bezierX = bezierCoords[0];
                        bezierY = bezierCoords[1];
                        if (internal_masses[i] == "0") {
                            ctx.beginPath();
                            ctx.setLineDash([10, 15]);
                            ctx.moveTo(coords[edges[i][0] -1][0], coords[edges[i][0] -1][1]);
                            ctx.quadraticCurveTo(bezierX, bezierY, coords[edges[i][1] -1][0], coords[edges[i][1] -1][1]);
                            ctx.stroke();
                        } else {
                            ctx.beginPath();
                            ctx.setLineDash([]);
                            ctx.moveTo(coords[edges[i][0] -1][0], coords[edges[i][0] -1][1]);
                            ctx.quadraticCurveTo(bezierX, bezierY, coords[edges[i][1] -1][0], coords[edges[i][1] -1][1]);
                            ctx.stroke();
                        }
                    }
                }

                // Draw external legs
                for (i=0; i < nodes.length; i++) {
                    if (external_masses[i] == "0") {
                        ctx.beginPath();
                        ctx.setLineDash([10, 15]);
                        ctx.moveTo(coords[nodes[i] - 1][0], coords[nodes[i] - 1][1]);
                        ctx.lineTo(externalLines[nodes[i] - 1][0], externalLines[nodes[i] - 1][1]);
                        ctx.stroke();
                    } else {
                        ctx.beginPath();
                        ctx.setLineDash([]);
                        ctx.moveTo(coords[nodes[i] - 1][0], coords[nodes[i] - 1][1]);
                        ctx.lineTo(externalLines[nodes[i] - 1][0], externalLines[nodes[i] - 1][1]);
                        ctx.stroke();
                    }
                }
            }

            var canvas = document.getElementById("Visualisation");
            var ctx = canvas.getContext("2d");
            var coords = []; 
            var externalLines = [];
            var points = 0;

            var xCentre = canvas.width / 2;
            var yCentre = canvas.height / 2;
            var size = 150;
            var externalSize = 300;
            var bezierSize = 50;
            var bezierX, bezierY;
            var bezierCoords;

            // Add event listeners to the buttons
            document.getElementById('visualiseButton').addEventListener('click', drawDiagram);
            document.getElementById('startButton').addEventListener('click', sendInputToJulia);
        </script>
    </body>
    </html>

    """

    api = interaction_API()
    window = webview.create_window("PLD GUI", html = html_content, js_api = api)
    webview.start()
