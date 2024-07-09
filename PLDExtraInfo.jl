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

if length(args) == 10
    subRules = convertStringToArray(args[10])
else
    subRules = []
end

internal_masses = convertStringToArray(args[3])
external_masses = convertStringToArray(args[4])
discs_file = args[8]
output_file = args[9]
codim_start = nothing
face_start = nothing

open(discs_file, "r") do file
    firstline = readline(file)

    matched = match(r"codim: (\d+), face: (\d+)/(\d+)", firstline)

    global codim_start = parse(Int, matched.captures[1])
    global face_start = parse(Int, matched.captures[2])
end


discs, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, load_output = discs_file, substitutions = subRules)

unique_discs, weight_list = discriminants_with_weights(U+F, discs, codim_start=codim_start, face_start=face_start)

genericEuler = maximum([getGenericEuler(U+F, pars, vars, :random) for i in 1:10])

disc_eulers = []

for disc in unique_discs
    push!(disc_eulers, maximum([getGenericEuler(U+F, pars, vars, disc) for i in 1:10]))
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

