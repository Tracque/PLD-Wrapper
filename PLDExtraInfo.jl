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
        write(file, "Weights: $(string(weight_list))\n")
        write(file, "Discriminant: $(string.(unique_discs[i]))\n\n\n")
    end
end

