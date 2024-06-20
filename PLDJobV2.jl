push!(LOAD_PATH,string(pwd() * "/src"))
using Pkg; Pkg.activate()

using Oscar
using HomotopyContinuation

using PLD

function convertStringToArray(str)
    # Parse the string as Julia code
    parsed_array = Meta.parse(str)

    # Evaluate the parsed expression
    result_array = eval(parsed_array)

    return result_array

end

#Function to help beautify printed array-like objects. Shamelessly stolen from https://discourse.julialang.org/t/printing-an-array-without-element-type/6731
function show_vector_sans_type(io, v::AbstractVector)
	print(io, "[")
	for (i, elt) in enumerate(v)
			i > 1 && print(io, ", ")
			if elt isa AbstractVector
				show_vector_sans_type(io, elt)
			else
				print(io, elt)
			end
	end
	print(io, "]")
end
show_vector_sans_type(v::AbstractVector) = show_vector_sans_type(stdout, v)

# Expecting filepath to an input file

inputfile = ARGS[1];

args = []

open(inputfile) do file
    for line in eachline(file)
        push!(args, line)
    end
end

edges = convertStringToArray(args[1]);
nodes = convertStringToArray(args[2]);

#Here, we create the symbolic variables for the masses
allowed_chars = ["m","M","q","Q","l","L","p","P"]

for i in 1:length(allowed_chars)
    for j in 1:length(edges)
        var_name = Symbol(allowed_chars[i], j)
        @eval HomotopyContinuation.ModelKit.@var $var_name
    end
end

internal_masses = convertStringToArray(args[3]);
external_masses = convertStringToArray(args[4]);
diagramName = args[5];

codimStart = parse(Int, args[6]);
faceStart = parse(Int, args[7]);

println("edges: $edges")
println("nodes: $nodes")
println("internal_masses: $internal_masses")
println("external_masses: $external_masses")
println("diagramName: $diagramName")
println()

flush(stdout)

if args[9] == "True"
    single_face = true
else
    single_face = false
end

if args[8] == "sym"

    U, F, pars, vars = getUF(edges, nodes, internal_masses=internal_masses, external_masses=external_masses)

    open(string(diagramName, "_info.txt"), "w") do file
        write(file, string(diagramName, "\n\n"))
        write(file, "edges: $(edges)\n")
        write(file, "nodes: $(nodes)\n")
        write(file, "internal_masses: "); show_vector_sans_type(file, internal_masses); write(file, "\n")
        write(file, "external_masses: "); show_vector_sans_type(file, external_masses); write(file, "\n\n")
        write(file, "schwinger parameters: "); show_vector_sans_type(file, vars); write(file, "\n")  #Confusingly, pars are the variables and vars are the parameters!
        write(file, "kinematic variables: "); show_vector_sans_type(file, pars); write(file, "\n\n")
        write(file, "U: $(U)\n")
        write(file, "F: $(F)\n\n")
    end

    PrincipleLandauDet, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, save_output=string(diagramName, ".txt"), codim_start = codimStart, face_start = faceStart, method = :sym, single_face = single_face)
else
    PrincipleLandauDet, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, save_output=string(diagramName, ".txt"), codim_start = codimStart, face_start = faceStart, method = :num, single_face = true)
end 

println("Finished.")