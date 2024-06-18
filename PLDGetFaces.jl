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
codim_start = parse(Int, args[6]);
face_start = parse(Int, args[7]);

println("edges: $edges")
println("nodes: $nodes")
println("internal_masses: $internal_masses")
println("external_masses: $external_masses")
println()
println("Finding codims and faces")

flush(stdout)

codim_array = []
face_array = []

U_oscar, F_oscar, pars_oscar, vars_oscar = getUF(edges, nodes; internal_masses = internal_masses, external_masses = external_masses)

f = U_oscar + F_oscar

P = newton_polytope(f)
Σ = normal_fan(P)

E2 = convert.(Array{Int64},collect(Oscar.AbstractAlgebra.exponent_vectors(f)))
mons2 = collect(Oscar.monomials(f))
coeff2 = collect(Oscar.coefficients(f))

min_codim = Oscar.ambient_dim(P) - Oscar.dim(P)
fvector = Oscar.f_vector(P)

if codim_start < 0
    codim_start = Oscar.ambient_dim(P)
end

for codim = codim_start : -1 : min_codim

    # Cones in this codimension
    codim_cones = cones(Σ,codim)
    if codim > min_codim && isnothing(codim_cones)
        continue
    end
    
    # Number of faces in this codimension
    noFaces = codim > min_codim ? length(codim_cones) : 1

    push!(codim_array, codim)
    push!(face_array, noFaces)

end

println(replace(string(codim_array), r"\w+\[([^\]]+)\]" => s"\1"))
println(replace(string(face_array), r"\w+\[([^\]]+)\]" => s"\1"))