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
codimStart = Parse(Int, args[5])
faceStart = Parse(Int, args[6])

if args[7] == "True"
    single_face = true
else
    single_face = false
end

discs_file = args[8]
output_file = args[9]
faces_done = []

poly, S, pars, vars = HC_to_oscar_S_mynames(poly, pars, vars; parnames = string.(pars), varnames = string.(pars))



open(discs_file, "r") do file

    for line in eachline(file)

        matched = match(r"codim: (\d+), face: (\d+)/(\d+)", line)

        push!(faces_done, [parse(Int, matched.captures[1]), parse(Int, matched.captures[2])])

    end
end

faces_done = sort(faces_done, by = x -> (-x[1], x[2]))

discs = []
open(load_output, "r") do f
    for line in eachline(f)
        # Extract discriminants using regex
        disc = match(r"discriminant: (.+)$", line)
        if disc != nothing
            push!(discriminants, eval(Meta.parse(replace("["*subscript_to_bracket(disc[1])*"]", "/" => "//"))))
        end
    end
end

unique_discs, weight_list = discriminants_with_weights(poly, discs, faces_done = faces_done)

genericEuler = maximum([getGenericEuler(poly, pars, vars, :random) for i in 1:10])

disc_eulers = []

for disc in unique_discs
    push!(disc_eulers, maximum([getGenericEuler(poly, pars, vars, disc) for i in 1:10]))
end

open(output_file, "a") do file
    write(file, "Generic Euler Discriminant, χ∗: $(genericEuler)\n\n\n")
    
    for i in length(unique_discs):-1:1
        write(file, "#################\n")
        write(file, "Discriminant $(length(unique_discs) - i + 1)\n")
        write(file, "#################\n\n")

        write(file, "Euler Discriminant, χ: $(disc_eulers[i])\n")
        write(file, "Weights: "); show_vector_sans_type(file, weight_list[i]); write(file, "\n");
        write(file, "Discriminant: $(string.(unique_discs[i]))\n\n\n")
    end
end

