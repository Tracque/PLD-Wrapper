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

#Create symbol variables for masses/kinematics
symbols_to_define = isolateElements(args[1])
symbols_to_define = vcat(symbols_to_define, isolateElements(args[2]))

for var_name in symbols_to_define
    if isnothing(tryparse(Int, string(var_name)))
        variable = Symbol(string(var_name))
        @eval HomotopyContinuation.ModelKit.@var $variable
    end
end

pars = convertStringToArray(args[1])
vars = convertStringToArray(args[2])

poly = eval(Meta.parse(subscript_to_bracket(args[3])))

diagramName = args[4]
codimStart = args[5]
faceStart = args[6]
single_face = args[7]

println("Starting...")
println()

flush(stdout)

poly, S, pars, vars = HC_to_oscar_S_mynames(poly, pars, vars; parnames = string.(pars), varnames = string.(pars))

if args[8] == "sym"

    open(string(diagramName, "_info.txt"), "w") do file
        write(file, "Custom Polynomial\n")
        write(file, string(diagramName, "\n\n"))
        write(file, "variables: "); show_vector_sans_type(file, vars); write(file, "\n")  #Confusingly, pars are the variables and vars are the parameters!
        write(file, "parameters: "); show_vector_sans_type(file, pars); write(file, "\n\n")
        write(file, "polynomial: $(poly)\n\n")
    end

    discs = getSpecializedPAD(poly, pars, vars, save_output=string(diagramName, ".txt"), codim_start = codimStart, face_start = faceStart, method = :sym, single_face = single_face)
else
    discs = getSpecializedPAD(poly, pars, vars, save_output=string(diagramName, ".txt"), codim_start = codimStart, face_start = faceStart, method = :num, single_face = true)
end 

println("Finished.")