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
discs_file = args[5];
output_file = args[6]

discs, pars, vars, U, F = getPLD(edges, nodes, internal_masses=internal_masses, external_masses=external_masses, load_output = discs_file)

unique_discs, weight_list = discriminants_with_weights(U+F, discs)

genericEuler = maximum([getGenericEuler(U+F, pars, vars, :random) for i in 1:10])

disc_eulers = []

for disc in unique_discs
    push!(disc_eulers, maximum([getGenericEuler(U+F, pars, vars, disc) for i in 1:10]))
end

open(output_file, "a") do file
    write(file, "Generic Euler Discriminant, χ∗: $(genericEuler)\n\n\n")
    
    for i in 1:length(unique_discs)
        write(file, "#################\n")
        write(file, "Discriminant $(i)\n")
        write(file, "#################\n\n")

        write(file, "Euler Discriminant, χ: $(disc_eulers[i])\n")
        write(file, "Weights: "); show_vector_sans_type(file, weight_list[i]); write(file, "\n");
        write(file, "Discriminant: $(string.(unique_discs[i]))\n\n\n")
    end
end

