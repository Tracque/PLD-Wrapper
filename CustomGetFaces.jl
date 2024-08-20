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
        variable = Symbol(string(var_name))push!(LOAD_PATH,string(pwd() * "/src"))
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
        
        codim_start = parse(Int, args[5])
        face_start = parse(Int, args[6])
        
        println("Finding codims and faces")
        
        flush(stdout)
        
        codim_array = []
        face_array = []
        
        f, S, pars_oscar, vars_oscar = HC_to_oscar_S_mynames(poly, pars, vars; parnames = string.(pars), varnames = string.(pars))
        
        
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
        @eval HomotopyContinuation.ModelKit.@var $variable
    end
end

pars = convertStringToArray(args[1])
vars = convertStringToArray(args[2])

poly = eval(Meta.parse(subscript_to_bracket(args[3])))

diagramName = args[4]

codim_start = parse(Int, args[5])
face_start = parse(Int, args[6])

println("Finding codims and faces")

flush(stdout)

codim_array = []
face_array = []

f, S, pars_oscar, vars_oscar = HC_to_oscar_S_mynames(poly, pars, vars; parnames = string.(pars), varnames = string.(pars))

P = newton_polytope(f)
Σ = normal_fan(P)

E2 = convert.(Array{Int64},collect(Oscar.AbstractAlgebra.exponent_vectors(f)))
mons2 = collect(Oscar.monomials(f))
coeff2 = collect(Oscar.coefficients(f))push!(LOAD_PATH,string(pwd() * "/src"))
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
        variable = Symbol(string(var_name))push!(LOAD_PATH,string(pwd() * "/src"))
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
        
        codim_start = parse(Int, args[5])
        face_start = parse(Int, args[6])
        
        println("Finding codims and faces")
        
        flush(stdout)
        
        codim_array = []
        face_array = []
        
        f, S, pars_oscar, vars_oscar = HC_to_oscar_S_mynames(poly, pars, vars; parnames = string.(pars), varnames = string.(pars))
        
        
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
        @eval HomotopyContinuation.ModelKit.@var $variable
    end
end

pars = convertStringToArray(args[1])
vars = convertStringToArray(args[2])

poly = eval(Meta.parse(subscript_to_bracket(args[3])))

diagramName = args[4]

codim_start = parse(Int, args[5])
face_start = parse(Int, args[6])

println("Finding codims and faces")

flush(stdout)

codim_array = []
face_array = []

f, S, pars_oscar, vars_oscar = HC_to_oscar_S_mynames(poly, pars, vars; parnames = string.(pars), varnames = string.(pars))

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