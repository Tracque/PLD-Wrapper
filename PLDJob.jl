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

function isolateElements(str)

    stripped_array = strip(str, ['[',']'])

    elems = split(stripped_array, ',')

    return strip.(elems)

end

function isolateSubSymbols(str)

    stripped_array = strip(str, ['[',']'])

    subrules = split(stripped_array, ',')

    allsymbols = []

    for i in 1:length(subrules)

        syms = split(subrules[i], "=>")

        leftsym = split.(syms[1])
        rightsym = split.(syms[2])

        push!(allsymbols, strip.(leftsym[1]))
        push!(allsymbols, strip.(rightsym[1]))

    end

    return allsymbols

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

#Create symbol variables for masses/kinematics
symbols_to_define = isolateElements(args[3])
symbols_to_define = vcat(symbols_to_define, isolateElements(args[4]))

if length(args) == 10
    symbols_to_define = vcat(symbols_to_define, isolateSubSymbols(args[10]))
end

for var_name in symbols_to_define
    if isnothing(tryparse(Int, string(var_name)))
        variable = Symbol(string(var_name))
        @eval HomotopyContinuation.ModelKit.@var $variable
    end
end

internal_masses = convertStringToArray(args[3])
external_masses = convertStringToArray(args[4])
diagramName = args[5]

codimStart = parse(Int, args[6])
faceStart = parse(Int, args[7])
if length(args) == 10
    subRules = convertStringToArray(args[10])
else
    subrules = []
end

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

    U, F, pars, vars = getUF(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, substitute = subRules)

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

    PrincipleLandauDet, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, save_output=string(diagramName, ".txt"), codim_start = codimStart, face_start = faceStart, method = :sym, single_face = single_face, substitutions = subRules)
else
    PrincipleLandauDet, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, save_output=string(diagramName, ".txt"), codim_start = codimStart, face_start = faceStart, method = :num, single_face = true, substitutions = subRules)
end 

println("Finished.")