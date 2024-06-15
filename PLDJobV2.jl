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

if args[8] == "sym"
    PrincipleLandauDet, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, save_output=string(diagramName, ".dat"), codim_start = codimStart, face_start = faceStart, method = :sym)
else
    PrincipleLandauDet, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, save_output=string(diagramName, ".dat"), codim_start = codimStart, face_start = faceStart, method = :num, single_face = true)
end 


# Open and write to file
raw"""
open(string(diagramName, ".txt"), "w") do file
    write(file, string(diagramName, "\n\n"))
    write(file, "edges: $(edges)\n")
    write(file, "nodes: $(nodes)\n")
    write(file, "internal_masses: $(internal_masses)\n")
    write(file, "external_masses: $(external_masses)\n\n")
    write(file, "schwinger parameters: $(vars)\n")  #Confusingly, pars are the variables and vars are the parameters!
    write(file, "kinematic variables: $(pars)\n\n")
    write(file, "U: $(U)\n")
    write(file, "F: $(F)\n\n")
    write(file, "discriminants: $(PrincipleLandauDet)\n")
end
"""

println("Finished.")