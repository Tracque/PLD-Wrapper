using Distributed

@everywhere begin

    push!(LOAD_PATH,string(pwd() * "/src")) #Adjust this to be the correct path to the PLD source files

    using ThreadingUtilities
    using Base.Threads: @spawn, fetch, nthreads, threadid
    using Base.Threads: Atomic, atomic_add!
    using Pkg; Pkg.activate()

    using Oscar
    using HomotopyContinuation
    using PLD

    using Blink #For GUI

    #Struct to help keep track of available threads
    mutable struct ThreadPool
        max_threads::Int
        available_threads::Atomic{Int}
        
        #Constructor for the struct
        function ThreadPool(max_threads::Int)
            return new(max_threads, Atomic{Int}(max_threads))
        end
    end

    #This struct stores some important variables we need to carry through the program
    mutable struct Ctx
        discriminants::Array
        codim_cones::Array
        codim_start::Int
        face_start::Int
        min_codim::Int
        single_face::Bool
        save_output::String

        #Contructor for the struct
        function Ctx(discriminants, codim_cones, codim_start, face_start, min_codim, single_face, save_output)
            return new(discriminants, codim_cones, codim_start, face_start, min_codim, single_face, save_output)
        end

    end

    #Method to check if the threads are available
    function acquire!(pool::ThreadPool, n::Int)
        available = atomic_add!(pool.available_threads, -n)
        if available >= 0
            return true
        else
            atomic_add!(pool.available_threads, n)
            return false
        end
    end

    #Method to release threads back into the pool
    function release!(pool::ThreadPool, n::Int)
        atomic_add!(pool.available_threads, n)
    end

    #This struct will allow us to manage the symbolic and numeric tasks
    mutable struct TaskPairManager
        pool::ThreadPool
        sym_tasks::Array
        num_tasks::Array
        task_start_times::Array
        task_not_aborted::Array
        task_codim::Array
        task_face::Array
        timeout::Int

        #Constructor of the struct
        function TaskPairManager(sym_tasks::Array, num_tasks::Array, task_start_times::Array, task_not_aborted::Array, task_codim::Array, task_face::Array; timeout::Int)
            return new(sym_tasks, num_tasks, task_start_times, task_not_aborted, task_codim, task_face, timeout)
        end

    end

    function start_pair(manager::TaskPairManager, IF, vars, pars, high_prec, codim, verbose, homogeneous, face)
        #Create and keep track of 2 new tasks (symbolic and numeric for this face)
        push!(manager.sym_tasks, @spawn getSpecializedDiscriminant(IF,vars,pars; method = :sym, high_prec = high_prec, codim = codim, verbose = verbose, homogeneous = homogeneous))
        push!(manager.num_tasks, @spawn getSpecializedDiscriminant(IF,vars,pars; method = :num, high_prec = high_prec, codim = codim, verbose = verbose, homogeneous = homogeneous))

        push!(manager.task_not_aborted, true)

        #Keep note of when we started each task pair
        push!(manager.task_start_times, time())

        #Keep note of the current face and codim
        push!(manager.task_codim, codim)
        push!(manager.task_face, face)
    end

    function check_tasks(manager::TaskPairManager, unique_discriminants::Array, Context::Ctx)
        for i = 1 : 1 : length(manager.sym_tasks)
            if istaskdone(manager.sym_tasks[i])

                if manager.task_not_aborted[i] #Only try to fetch results if the process finished normally, not if we killed it ourselves

                    #Get, store and display output
                    disc = fetch(manager.sym_tasks[i])
                    unique_discriminants = store_and_output_discs(disc, manager.task_codim[i], manager.task_face[i], unique_discriminants, Context)

                    #Kill the other thread
                    Base.throwto(manager.num_tasks[i], InterruptException())

                    #Remove the now dead tasks from the tracking arrays
                    deleteat!(manager.sym_tasks, i)
                    deleteat!(manager.num_tasks, i)
                    deleteat!(manager.task_start_times, i)
                    deleteat!(manager.task_not_aborted, i)

                    release!(manager.pool, 2)

                    faces_done += 1
                
                end
            
            #The symbolic thread shouldn't take long if it succeeds, so we time it out after a while
            elseif manager.task_not_aborted[i] && (time() - manager.task_start_times[i] > manager.timeout)

                #Kill the symbolic thread if it has been too long
                Base.throwto(manger.sym_tasks[i], InterruptException())

                #Release the symbolic thread
                manager.pool.release(1)

                #Update the status of this thread
                manager.task_not_aborted[i] = false

            elseif istaskdone(manager.num_tasks[i])
                
                #Get, store and display output
                disc = fetch(manager.num_tasks[i])
                unique_discriminants = store_and_output_discs(disc, manager.task_codim, manager.task_face, unique_discriminants, Context)

                #Kill the other thread
                Base.throwto(manager.sym_tasks[i], InterruptException())

                #Remove the now dead tasks from the array
                deleteat!(manager.sym_tasks, i)
                deleteat!(manager.num_tasks, i)
                deleteat!(manager.task_start_times, i)

                #If the symbolic task has already been killed, then it has already been released
                if manager.task_not_aborted[i]
                    release!(manager.pool, 1)
                else
                    release!(manager.pool, 2)
                end

                deleteat!(manager.task_not_aborted, i)

                faces_done += 1

            end
        end

        return unique_discriminants
    end

    #Slightly modified version of the output code in getSpecializedPAD to allow for verbose output to the GUI
    #(and to account for some variable name changes!)
    function store_and_output_discs(disc, codim, face, unique_discriminants, Context)

        if codim > Context.min_codim
            weight = sum(Oscar.rays(Context.codim_cones[codim][face]))
        else
            weight = zeros(Int,Oscar.ambient_dim(P))
        end

        len_before = length(unique_discriminants)
        push!(Context.discriminants, disc)

        disc_string = [replace(s, "//" => "/") for s in string.(vcat(disc...))]

        append!(unique_discriminants, disc_string)
        unique_discriminants = sort(unique(unique_discriminants))
        if codim > Context.min_codim
            noFaces = length(Context.codim_cones[Context.codim_start - codim + 1])
        else
            noFaces = 1
        end

        #printstyled("codim: $(codim), face: $(face)/$(noFaces), weights: [$(join(string.(vcat(weight...)),", "))], discriminant: $(join(disc_string,", "))\n"; color = :red)
        #if !Context.single_face && isnothing(single_weight) && verbose println("") end
        if Context.save_output != ""
            open(Context.save_output, "a") do file
                print(file, "codim: $(codim), face: $(face)/$(noFaces), weights: [$(join(string.(vcat(weight...)),", "))], discriminant: $(join(disc_string,", "))\n")
            end
        end

        len_after = length(unique_discriminants)

        #if !Context.single_face && isnothing(single_weight) && len_after > len_before
        if !Context.single_face && len_after > len_before
            #printstyled("New discriminants after codim $(codim), face $(face)/$(noFaces). The list is: $(join(unique_discriminants, ", "))\n"; color = :yellow)
            #if verbose println("") end
        end

        flush(stdout)

        weight_string = string([(join(string.(vcat(weight...)),", "))])
        disc_string_out = (join(disc_string,", "))

        #Display also to Blink Window
        #output_string = string("Most recent PLD calculation was codim: ", string(codim), ", face: ", string(face), "/", string(noFaces), ", weights: ", weight_string,  ", discriminant: ", disc_string_out)
        #@js display_window displayJuliaOutput($output_string)

        return unique_discriminants

    end
end

#This function creates the web-based GUI for using PLD
function openWindow()

    global display_window = Window()

    # Define the HTML content
    html_content = raw"""
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
        justify-content: space-around;
        width: 90%;
        max-width: 1600px;
        margin-bottom: 20px; /* Add some space between the container and the third box */
    }
    
    .box {
        background-color: #fff;
        border: 2px solid #ddd;
        border-radius: 8px;
        padding: 20px;
        flex: 1;
        margin: 10px;
        text-align: center;
        transition: all 0.3s ease;
    }
    
    .box:hover {
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        transform: translateY(-5px);
    }
    
    #box3 {
        width: 90%;
    }
    
    @media (max-width: 768px) {
        .container {
            flex-direction: column;
        }
    
        #box3 {
            width: 100%;
        }
    }
    
        </style>
    </head>
    <body>

    <div class="container">

        <div class="box" id="box1">
            <h1 style="font-size:30px">PLD User Interface</h1>

        <p style="line-height: 1.2">
        Keep your inputs to the suggested format to ensure compatibility with PLD.jl <br>
        I have tried to allow for many possible variable labels for the squared masses, including p, P, m, M, q, Q, l and L. (as in m1, m2, m3,... or similar) <br>
        You should ensure you align the edges/internal masses and the nodes/external masses. For example, this means if your edges are input as [[1,2],...], then the internal mass array [m1,...] will set the mass of the [1,2] edge to m1.
        </p>
        <br>

            Edges <input type="text" placeholder="e.g. [[1,2],[2,3],[1,3]]" id="edges"> <br>
            Internal Squared Masses <input type="text" placeholder="e.g. [m1,0,m1,m2]" id="internals"> <br>
            Nodes  <input type="text" placeholder="e.g. [1,2,3]" id="nodes"> <br>
            External Squared Masses <input type="text" placeholder="e.g. [p1,p1,p2,0]" id="externals"> <br>
            <button id="visualiseButton">Generate Diagram Visualisation</button>
            <br><br>

        <p style="line-height: 1.2">
            Note: The provided visualisation is algorithmically generated and is not supposed to produce literature ready diagrams.<br>
            The visualisation places all vertices at the vertices of a regular polygon and will not avoid intersecting lines for non-planar diagrams. <br>
            Dashed lines represent massless propagators. Solid lines represent massive propagators.<br>
        The vertex labelled "1" will be the rightmost vertex, with subsequent vertices being labelled clockwise. <br>
        </p>
    
        <br>

        <div style="float:left">
        <input type="radio" id="continuingCheckbox" name="runType" value="continue">
        <label for="continuingCheckbox"> Continue an existing calculation?</label><br>
        Starting codim <input type="text" id="startC"> <br>
        Starting face <input type="text" id="startF"> <br>
        </div>

        <div style="float:right">
        <input type="radio" id="startCheckbox" name="runType" value="singleFace">
        <label for="startCheckbox"> Calculate singularities on a single face?</label><br>
        Codim <input type="text" id="startC"> <br>
        Face <input type="text" id="startF"> <br>
        </div>

        <br><br><br><br><br><br>

        <input type="radio" id="normalCheckbox" name="runType" value="normal">
        <label for="normalCheckbox"> Run a full PLD calculation?</label><br>

        Below, input the desired output file name or path and name for the calculation.

        File Path <input type="text" id="filePath"> <br>

        <button id="startButton">Start Calculation</button>

        

    </div>
        

    <div class="box" id="box2">
        <canvas id="Visualisation" height="800" width="800"></canvas> <br>


    </div>

    </div>

    <div class="container-below">

    <div class="box" id="box3">
        <h2>Output:</h2> <br>
        <p style="font-size:20px" id="outputPara">Output will appear here.</p>
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

                
                if (document.getElementById("continuingCheckbox").checked == true) {
                    var startCodim = document.getElementById("startC").value;
                    var startFace = document.getElementById("startF").value;
                    Blink.msg("Start Calc", [edgeInputs, nodeInputs, internalInputs, externalInputs, startCodim, startFace, outputFilePath, "Continue"]);
                } else if (document.getElementById("startCheckbox").checked == true) {
                    var startCodim = document.getElementById("startC").value;
                    var startFace = document.getElementById("startF").value;
                    Blink.msg("Start Calc", [edgeInputs, nodeInputs, internalInputs, externalInputs, startCodim, startFace, outputFilePath, "Single Face"]);
                } else {
                    Blink.msg("Start Calc", [edgeInputs, nodeInputs, internalInputs, externalInputs, outputFilePath]);
                }

                // Show that we have started successfully
                document.getElementById("outputPara").innerText = "PLD is starting. Output should appear here soon...";
            }

            function displayJuliaOutput(outputTxt) {
                document.getElementById("outputPara").innerText = outputTxt;
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

    body!(display_window, html_content)

    handle(display_window, "Start Calc") do args

        if length(args) == 8
            receiveInput(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8])
        else
            receiveInput(args[1], args[2], args[3], args[4], args[5])
        end

    end

    

end

function convertStringToArray(str)
    # Parse the string as Julia code
    parsed_array = Meta.parse(str)

    # Evaluate the parsed expression
    result_array = eval(parsed_array)

    return result_array

end

function receiveInput(edgeInputs, nodeInputs, internalInputs, externalInputs, outputFilePath, startCodim = -1, startFace = 1, runType = "Normal")

    #First, we parse the edges and nodes
    edges = convertStringToArray(edgeInputs)
    nodes = convertStringToArray(nodeInputs)

    #Here, we create the symbolic variables for the masses
    allowed_chars = ["m","M","q","Q","l","L","p","P"];

    for i in 1:length(allowed_chars)
        for j in 1:length(edges)
            var_name = Symbol(allowed_chars[i], j)
            @eval HomotopyContinuation.ModelKit.@var $var_name
        end
    end

    #Now we can parse the rest of our inputs
    internal_masses = convertStringToArray(internalInputs)
    external_masses = convertStringToArray(externalInputs)

    if (isa(startCodim, String))
        codim_start = convertStringToArray(startCodim)
        face_start = convertStringToArray(startFace)
    else
        codim_start = startCodim
        face_start = startFace
    end

    #With all our inputs handled, we pass over to the main function
    main(edges, nodes, internal_masses, external_masses, codim_start, face_start, outputFilePath, runType)
end

function main(edges, nodes, internal_masses, external_masses, codim_start, face_start, outputFilePath, runType)

    #The simplest runType is single face
    if (runType == "Single Face")
        Contx = Ctx([], [], codim_start, face_start, 0, true, outputFilePath)
        PLD, kinematic_vars, schwinger_pars, U, F = getPLDMultithreaded(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, Context = Contx)
    else
        Contx = Ctx([], [], codim_start, face_start, 0, false, outputFilePath)
        PLD, kinematic_vars, schwinger_pars, U, F = getPLDMultithreaded(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, Context = Contx)
    end

end

#This function is an exact copy of getPLD, but using getPADMultithreaded
function getPLDMultithreaded(edges, nodes; internal_masses = :zero, external_masses = :zero, high_prec = false, single_weight = nothing, verbose = false, homogeneous = true, save_output = "", load_output = "", Context = Ctx([],[],-1,1,0,false, "test.txt"))

    U_oscar, F_oscar, pars_oscar, vars_oscar = getUF(edges, nodes; internal_masses = internal_masses, external_masses = external_masses)

    if load_output == ""
        discriminants = getPADMultithreaded(U_oscar + F_oscar, pars_oscar, vars_oscar; high_prec = high_prec, single_weight = single_weight, verbose = verbose, homogeneous = homogeneous, save_output = save_output, Context)
    else
        pars_string = [subscript_to_bracket(str) for str in string.(pars_oscar)]

        if internal_masses == :generic
            global m = Vector{fmpq_mpoly}(undef, length(edges))
        end
        if external_masses == :generic
            global M = Vector{fmpq_mpoly}(undef, length(nodes))
        end

        global R, pars = PolynomialRing(QQ, pars_string)
        eval(Meta.parse("global ("*join(pars_string, ", ")*") = pars"))

        discriminants = []
        open(load_output, "r") do f
            for line in eachline(f)
                # Extract discriminants using regex
                disc = match(r"discriminant: (.+)$", line)
                if disc != nothing
                    push!(discriminants, eval(Meta.parse(replace("["*subscript_to_bracket(disc[1])*"]", "/" => "//"))))
                end
            end
        end
    end

    return discriminants, pars_oscar, vars_oscar, U_oscar, F_oscar

end

#This function is an exact copy of getSpecializedPAD(), but with multithreading implemented for the loop over faces and some extra output code to update the GUI.
function getPADMultithreaded(f, pars, vars; high_prec = false, single_weight = nothing, verbose = false, homogeneous = true, save_output = "", Context = Ctx([],[],-1,1,0,false, "test.txt"))

    addprocs(nthreads() - 2)

    P = newton_polytope(f)
    Σ = normal_fan(P)

    E2 = convert.(Array{Int64},collect(Oscar.AbstractAlgebra.exponent_vectors(f)))
    mons2 = collect(Oscar.monomials(f))
    coeff2 = collect(Oscar.coefficients(f))

    Context.discriminants = []
    unique_discriminants = []

    Context.min_codim = Oscar.ambient_dim(P) - Oscar.dim(P)
    fvector = Oscar.f_vector(P)

    if verbose
        println("f = $(f)");
        println("pars = [$(split(string(pars),'[')[2])");
        println("vars = [$(split(string(vars),'[')[2])");
        println("method = $(method)");
        println("high_prec = $(high_prec)");
        println("codim_start = $(Context.codim_start)");
        println("face_start = $(Context.face_start)");
        println("single_face = $(Context.single_face)");
        println("single_weight = $(single_weight)");
        println("verbose = $(verbose)");
        println("homogeneous = $(homogeneous)");
        println("save_output = $(Context.save_output)");
        println("f_fector = [$(split(string(fvector),'[')[2])\n");
    end

    if Context.codim_start < 0
        Context.codim_start = Oscar.ambient_dim(P)
    end

    if single_weight !== nothing
        Context.single_face = true
        weights = getWeights(f)
        indices = [(i, j) for i in 1:length(weights) for j in findall(x -> x == single_weight, weights[i])]
        if length(indices) > 0
            Context.codim_start = indices[1][1] - 1
            Context.face_start = indices[1][2]
        else
            printstyled("The specified weight $(single_weight) was not found\n"; color = :yellow)
            return
        end
    end

    #Set starting point
    codim = Context.codim_start
    face = Context.face_start

    #Store some info about the diagram
    faces_done = 0
    faces_todo = 0

    #Create a struct to manage multithreading of calculations
    PLD_manager = TaskPairManager(nthreads()-2, [], [], [], [], [], [], timeout=300)

    #CALCULATE TOTAL FACES TO BE FOUND
    for i = Context.codim_start : -1 : Context.min_codim + 1
        push!(Context.codim_cones, cones(Σ, i))
        faces_todo += length(Context.codim_cones[Context.codim_start - i + 1])
    end

    #One extra face for the minimum codim
    faces_todo += 1

    #Loop over all faces we need to compute
    while faces_done < faces_todo

        if acquire!(PLD_manager.pool, 2)

            #Find the arguments of getSpecializedDiscriminant()
            if codim > Context.min_codim
                weight = sum(Oscar.rays(Context.codim_cones[codim][face]))
                vals = convert.(Rational{Int64}, [transpose(weight)*e for e in E2] )
                mininds = findall(weight -> weight == minimum(vals), vals)
                IF = transpose(coeff2[mininds])*mons2[mininds]
            else
                weight = zeros(Int,Oscar.ambient_dim(P))
                IF = f
            end

            start_pair(PLD_manager, IF, vars, pars, high_prec, codim, verbose, homogeneous, face)

            #Display to GUI
            #output_string = "Started calculation on codim: " * string(codim) * ", face: " * string(face)
            #@js display_window displayJuliaOutput($output_string)

            #Move on to the next face
            if (face < length(Context.codim_cones[Context.codim_start - codim + 1]))
                face += 1
            else
                face = 1
                codim -= 1
            end

        end

        #This checks if any calculations have finished and handles their output if they have
        unique_discriminants = check_tasks(PLD_manager, unique_discriminants, Context)

        #Sleep to avoid busy watching
        sleep(0.1)

    end 

    return(discriminants)
    
end

openWindow()