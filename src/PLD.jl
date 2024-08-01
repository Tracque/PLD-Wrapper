module PLD

# Authors: Claudia Fevola, Sebastian Mizera, Simon Telen
# Date: November 14, 2023
# Short description: This code accompanies the papers
# [PLD1] Claudia Fevola, Sebastian Mizera, and Simon Telen, "Landau Singularities Revisited", arXiv:
# [PLD2] Claudia Fevola, Sebastian Mizera, and Simon Telen, "Principal Landau Determinants", arXiv:

using HomotopyContinuation
using LinearAlgebra
using GenericSVD
using Oscar
using Printf

include("Landau.jl")

export normalize_coeffs,
oscar_to_HC_Q,
oscar_to_HC_S,
HC_to_oscar_Q,
HC_to_oscar_S,
HC_to_oscar_S_mynames,
assess,
project_codim1,
discriminants_with_weights,
getSpecializedDiscriminant,
getIF,
getWeights,
getSpecializedPAD,
getUF,
getEuler,
getPLD,
getGenericEuler,
EulerDiscriminantQ

# Computes the initial forms of a polynomial (in Oscar format) for a given set of weight vectors
# -------------  Input:
# f         the polynomial we want the initial forms of
# weights   weights has a list of vectors in the interior of some of the cones in the normal fan of Newt(f).
# -------------  Output:
# initial_forms      a list of all initial forms of f corresponding to faces corresponding to the vectors in weights
function getIF(f,weights)
    initial_forms = []
    E = convert.(Array{Int64},collect(Oscar.AbstractAlgebra.exponent_vectors(f)))
    mons = collect(Oscar.monomials(f))
    coeff = collect(Oscar.coefficients(f))
    for w in weights
        vals = convert.(Rational{Int64}, [transpose(w)*e for e in E] )
        mininds = findall(w -> w == minimum(vals), vals)
        push!(initial_forms,transpose(coeff[mininds])*mons[mininds])
    end
    return initial_forms
end

# Computes the weights of the Newton polytope of a polynomial (in Oscar format)
# ------------ Input:
# f            the polynomial we want the weights of
# ------------ Output:
# weights      weights[k] has a list of vectors in the interior of each cone of codimension k-1 in the normal fan of Newt(f).
function getWeights(f)
    P = newton_polytope(f)
    Σ = normal_fan(P)
    min_codim = Oscar.ambient_dim(P) - Oscar.dim(P)
    weightlist = []
    for codim = 0:Oscar.ambient_dim(P)
        weightlist = push!(weightlist,[])
        if codim > min_codim

            cones_codim = cones(Σ,codim)
            if isnothing(cones_codim)
                continue
            end
            cones_codim = collect(cones_codim)

            noFaces = length(cones_codim)
            for i in 1:noFaces
                weight = sum(Oscar.rays(cones_codim[i]))
                push!(weightlist[codim+1],convert(Vector{Int64},numerator.(lcm(denominator.(weight)).*weight)))
            end
        elseif codim == min_codim
            weight = zeros(Int,Oscar.ambient_dim(P))
            push!(weightlist[codim+1],convert(Vector{Int64},weight))
        else
            push!(weightlist[codim+1],[])
        end
    end
    convert.(Vector{Vector{Int64}},weightlist)
end

# Scale a polynomial f so that its leading term has coefficient 1.
# -------------  Input:
# f         the polynomial we want to normalize
# -------------  Output:
# scaled_f  the scaled version of f
function normalize_coeffs(f)
    cffs = collect(Oscar.coefficients(f))
    scaled_f = f//cffs[1]
    return scaled_f
end

# Converts a polynomial in Oscar.jl format, with coefficients in Q, into a polynomial in HomotopyContinuation.jl format.
# -------------  Input:
# f         the polynomial with coefficients in Q we want to convert
# vars      the HomotopyContinuation variables we want to use for the converted polynomial
# -------------  Output:
# a HomotopyContinuation.jl polynomial representing f
function oscar_to_HC_Q(f,vars)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.AbstractAlgebra.exponent_vectors(f))
    sum([cffs[i]*prod(vars.^exps[i]) for i = 1:length(cffs)])
end

# Converts a polynomial in Oscar.jl format, with coefficients in a polynomial ring, into a polynomial in HomotopyContinuation.jl format.
# -------------  Input:
# f         the polynomial with coefficients in a polynomial ring we want to convert
# x         the HomotopyContinuation.jl variables we want to use for the converted polynomial
# a         the HomotopyContinuation.jl variables we want to use to replace the variables of S
# -------------  Output:
# a HomotopyContinuation.jl polynomial representing f
function oscar_to_HC_S(f,x,a)
    cffs_oscar = collect(Oscar.coefficients(f))
    cffs_HC = [oscar_to_HC_Q(g,a) for g in cffs_oscar]
    exps = collect(Oscar.AbstractAlgebra.exponent_vectors(f))
    mons_HC = [prod(x.^exps[i]) for i = 1:length(exps)]
    sum([cffs_HC[i]*mons_HC[i] for i = 1:length(cffs_HC)])
end


# Converts a polynomial in HomotopyContinuation.jl format, with coefficients in Q, into a polynomial in Oscar.jl format.
# -------------  Input:
# f         the polynomial with coefficients in Q we want to convert
# vars      the Oscar variables we want to use for the converted polynomial
# HC_vars   the HomotopyContinuation.jl variables of f
# -------------  Output:
# an Oscar.jl polynomial representing f
function HC_to_oscar_Q(f,vars,HC_vars)
    E,C = exponents_coefficients(f,HC_vars)
    sum([C[i]*prod(vars.^E[:,i]) for i = 1:length(C)])
end

# Converts a polynomial in HomotopyContinuation.jl format, with coefficients in a polynomial ring, into a polynomial in Oscar.jl format.
# -------------  Input:
# f         the polynomial with coefficients in Q we want to convert
# pars      the HomotopyContinuation.jl parameters of f
# vars      the HomotopyContinuation.jl variables of f
# -------------  Output:
# an Oscar.jl polynomial representing f
# S         the ring in which the new polynomial lives
function HC_to_oscar_S(f,pars,vars)
    E, C = exponents_coefficients(f,vars)
    n = length(vars)
    m = length(pars)
    R, a = polynomial_ring(QQ, ["a$i" for i = 1:m]) # Polynomial ring in the parameters
    S, (x) = laurent_polynomial_ring(R, ["x$i" for i = 1:n]) # Polynomial ring with coefficients in R
    cffs_oscar = [HC_to_oscar_Q(g,a,pars) for g in C]
    sum([cffs_oscar[i]*prod(x.^E[:,i]) for i = 1:length(cffs_oscar)]), S
end

# Converts a polynomial in HomotopyContinuation.jl format, with coefficients in a polynomial ring, into a polynomial in Oscar.jl format.
# -------------  Input:
# f         the polynomial with coefficients in Q we want to convert
# vars      the Oscar variables we want to use for the converted polynomial
# HC_vars   the HomotopyContinuation.jl variables of f
# -------------  Output:
# an Oscar.jl polynomial representing f
# S         the ring in which the new polynomial lives
function HC_to_oscar_S_mynames(f,pars,vars; parnames = [], varnames = [])
    E, C = exponents_coefficients(f,vars)
    C = convert.(Expression, C)
    n = length(vars)
    m = length(pars)
    if !isempty(parnames)
        R, a = polynomial_ring(QQ, parnames)
    else
        R, a = polynomial_ring(QQ, ["a$i" for i = 1:m]) # Polynomial ring in the parameters
    end
    if !isempty(varnames)
        S, (x) = laurent_polynomial_ring(R, varnames)
    else
        S, (x) = laurent_polynomial_ring(R, ["x$i" for i = 1:n]) # Polynomial ring with coefficients in R
    end
    cffs_oscar = [HC_to_oscar_Q(g,a,pars) for g in C]
    sum([cffs_oscar[i]*prod(x.^E[:,i]) for i = 1:length(cffs_oscar)]), S, a, x
end

# Computes the discriminant of a polynomial F with variables x and parameters a. That is, it computes
# the defining equations of the hypersurface in a-space obtained by projecting the incidence variety Y = {F(x,a) = ∇F(x,a) = 0, prod(x)*y - 1 = 0}.
# The last equation, which uses an extra variable y, makes sure we only account for critical points in the torus.
# For each component of Y that projects to a hypersurface in parameter space, we compute the defining equation.
# -------------  Input:
# F                     the polynomial we want the discriminant of
# x                     variables of F
# a                     parameters of F. F is assumed to be a polynomial in both a and x
# method                use symbolic or numerical algorithm (:sym or :num), default = :sym
# high_prec             use high precision, default = false
# codim                 codimension of the Newton polytope of f. This is useful to estimate the fiber dimension, default = 0
# verbose               print intermediate results, default = false
# homogeneous           use homogeneous ansatz for polynomial reconstruction, default = true
# -------------  Output:
# discs                 a list of discriminants in Oscar format. These are polynomials in a.
function getSpecializedDiscriminant(F,x,a; method = :sym, high_prec = false, codim = 0, verbose = false, homogeneous = true)
    coeffs = collect(Oscar.coefficients(F))
    if length(coeffs) == 1
        return [normalize_coeffs(coeffs[1])]
    else
        r = length(x)
        T = polynomial_ring(QQ,vcat(vcat(["α$i" for i = 1:r],["y"]),string.(a)))
        tvars = T[2]
        exps = collect(Oscar.AbstractAlgebra.exponent_vectors(F))
        new_mons = [prod(tvars[1:r].^e) for e in exps]
        new_coeffs = [Oscar.evaluate(c,tvars[r+2:end]) for c in coeffs]
        newF = transpose(new_coeffs)*new_mons
        ders = [derivative(newF, j) for j = 1:r]
        eqs = [newF,ders...,prod(tvars[1:r+1])-1][:]

        # Symbolic method
        if method == :sym
            I = radical(ideal(eqs))
            PD = primary_decomposition(I)
            if verbose && !isempty(PD)
                println("Incidence variety: $(join([string(p[2]) for p in PD],", "))")
            end
            J = [eliminate(p[2],tvars[1:r+1]) for p in PD]
            discs = gens.(J[findall(k->length(gens(k))==1 && gens(k)[1]!=0, J)])
            dominant_components = findall(k->length(gens(k))==1 && gens(k)[1]==0, J)

            if length(dominant_components) >= 1
                if verbose
                    printstyled("This face has dominant components\n"; color = :magenta)
                end
            end
            if length(discs) >= 1

                discs_return = []

                for d in discs
                    
                    exps = collect(Oscar.AbstractAlgebra.exponent_vectors(d[1]))
                    new_mons = [prod(a.^(e[r+2:end])) for e in exps]
                    coeffs = collect(Oscar.coefficients(d[1]))
                    new_d = transpose(coeffs)*new_mons

                    push!(discs_return, normalize_coeffs(new_d))

                end

                return discs_return
            else
                return 1
            end

        # Numerical method
        elseif method == :num
            @var y[1:r+1] p[1:length(a)]
            eqs_HC = [oscar_to_HC_Q(eq,[y;p]) for eq in eqs]
            vv = variables(eqs_HC[1])
            y = [vv[findall(v-> v in y, vv)];y[end]]
            pp = vv[findall(v-> v in p, vv)]
            indsinp = [findfirst(k->k==e,p) for e in pp]
            aa = a[indsinp]
            p = pp
            eqs_HC = [eqs_HC[1]; differentiate(eqs_HC[1],y[1:end-1]); prod(y)-1]
            if isnothing(findfirst(eq -> HomotopyContinuation.degree(eq)==0 && eq !=0, eqs_HC)) && length(p)>0
                interpolants_per_dim, samples, gaps = project_codim1(eqs_HC,p,y; HP = high_prec, first_fd = codim-(r+1-length(y)), verbose = verbose, homogeneous = homogeneous)
                discs = []
                for i = 1:length(interpolants_per_dim)
                    if !isempty(interpolants_per_dim[i])

                        gooddiscs = findall(j -> gaps[i][j] > 1e7, 1:length(gaps[i]))
                        newdiscs = [HC_to_oscar_Q(interpolants_per_dim[i][j],aa,p) for j in gooddiscs]
                        nonfactorized = findall(j -> length(Oscar.factor(newdiscs[j])) == 1, 1:length(newdiscs)) 
                        
                        # Drop reducible components
                        if verbose && length(nonfactorized) < length(newdiscs)
                            printstyled("Dropping incorrectly sampled components\n"; color = :magenta)
                        end
                        newdiscs = newdiscs[nonfactorized]

                        discs = vcat(discs,newdiscs)
                    end
                end
                discs = [normalize_coeffs(dd) for dd in discs]
            else
                discs = [1]
            end
            if isempty(discs)
                discs = [1]
            end
            return discs

        else
            println("choose method = :sym or :num")
        end
    end
end

# Compute a numerical irreducible decomposition (see below) of the variety V defined as follows. Let x = vars, p = pars.
# Consider the incidence variety Y = {F(x,p) = 0}, where F is a system of polynomials.
# For each irreducible component X of Y, consider the projection X -> πₚ(X) onto the p coordinates.
# Collect those components X for which the generic fiber dimension of this projection is fiberdim (an optional input, default = 0),
# and the image πₚ(X) is a hypersurface in p-space. Among those, keep only the components X that are generically smooth (no scheme structure). The union all these components X is V.
# "Numerical irreducible decomposition" means that we compute degree many points on each irreducible component of V.
# The function returns a list of solutions on each component. Together with equations having these points as regular solutions.
# It also returns a the number of singular solutions n_sing of these equations, and a list degs with the degrees of the projections πₚ(X).
# -------------  Input:
# F                 a system of polynomial equations with variables vars and parameters pars
# pars              parameters of F
# vars              variables of F
# fiberdim          fiber dimension as explained above, default = 0
# verbose           print intermediate output, default = false
# -------------  Output:
# groups,           groups of solutions on different components of the incidence variety
# incidence_param,  parameterized system used to compute these solutions
# line,             instance of the line in parameter space that was used
# fibercut,         instance of the linear space in variable space that was used
# degs,             estimated degrees of the projections of the components
# n_sing            number of singular solutions. If nonzero, we should look at higher fiber dimension
function assess(F,pars,vars; fiberdim = 0, verbose = false)

    @var a[1:length(pars)-1,1:length(pars)+1] # parameters for a line in parameter space
    @var b[1:fiberdim, 1:length(vars)+1] # parameters for a linear cut in variable space
    incidence_param = System([F;a*[pars;1];b*[vars;1]], variables = [pars;vars], parameters = [a[:];b[:]])
    
    # A first target instance
    line = randn(ComplexF64,length(pars)-1,length(pars)+1)
    fibercut = randn(ComplexF64,fiberdim,length(vars)+1)

    local regsols, n_sing
    try
        if length(pars)-1 == 0 && fiberdim == 0
            R = HomotopyContinuation.solve(incidence_param; compile = false, threading = true, show_progress = verbose)
        else
            R = HomotopyContinuation.solve(incidence_param, target_parameters = [line[:];fibercut[:]]; compile = false, threading = true, show_progress = verbose)
        end
        regsols = solutions(R)
        n_sing = length(HomotopyContinuation.singular(R))
    catch e
        if verbose
            printstyled("HomotopyContinuation.jl failed to solve the system. This face needs to be reexamined!\n"; color = :yellow)
        end
        regsols = []
        n_sing = 0
    end

    keep = []
    discard = []
    Jac = differentiate(expressions(incidence_param)[1:length(F)],variables(incidence_param))
    p = length(pars)

    printed_dominant = false
    for sol in regsols
        valJac = ComplexF64.(HomotopyContinuation.evaluate.(Jac,variables(incidence_param) => sol, parameters(incidence_param)=>[line[:];fibercut[:]]))
        k = nullspace(valJac)
        if (size(k,2) > 1 && rank(k[1:p,:], rtol = 1e-8) == p) || (size(k,2) == 1 && p == 1 && norm(k[1:p]) > 1e-6)
            push!(discard,sol)
            if verbose && !printed_dominant
                printstyled("This face has dominant components\n"; color = :magenta)
                printed_dominant = true
            end
        else
            push!(keep,sol)
        end
    end
    regsols = keep
    if verbose
            println("   no. regular solutions: $(length(regsols))")
            println("   no. singular solutions: $n_sing")
            println("   discarded $(length(discard)) regular solutions on a dominant component")
    end

    # Make groups of solutions using monodromy
    groups = []
    if !isempty(regsols)
        for sol in regsols
            loop = true
            for g in groups
                mindist = minimum(norm.([sol] .- g))
                if mindist < 1e-6
                    loop = false
                end
            end
            if loop
                try
                    MS = monodromy_solve(incidence_param, [sol], [line[:];fibercut[:]]; compile = false, threading = true, show_progress = verbose, target_solutions_count = 10000)
                    push!(groups,solutions(MS))
                catch e
                    if verbose
                        printstyled("HomotopyContinuation.jl failed to solve the system. This face needs to be reexamined!\n"; color = :yellow)
                    end
                end
            end
        end
    end
    degs = []
    for g in groups
        if length(g) == 0
        elseif length(g) == 1
            push!(degs, 1)
        else
            push!(degs,length(unique_points([sol[1:length(pars)] for sol in g], rtol=1e-6)))
        end
    end
    return groups, # groups of solutions on different components of the incidence variety
           incidence_param, # parameterized system used to compute these solutions
           line, # instance of the line in parameter space that was used
           fibercut, # instance of the linear space in variable space that was used
           degs, # estimated degrees of the projections of the components
           n_sing # number of singular solutions. If nonzero, we should look at higher fiber dimension
end

# Let x = vars, p = pars. This function computes the equations for the hypersurface V obtained as follows.
# Let Y = {LE(x;p) = 0} be the incidence variety in (x,p)-space. Let X be an irreducible component,
# such that the projection πₚ(X) onto parameter space is a hypersurface. Among those, keep only the components X
# which are generically smooth (no scheme structure). V is the union of the projections of all such X.
# By default, the code does not attempt to interpolate components for which it needs more than 10000 samples.
# This bound can be altered manually in the code.
# -------------  Input:
# LE                            a system of polynomial equations with variables vars and parameters pars
# pars                          parameters of F
# vars                          variables of F
# HP                            use high precision, default = false
# first_fd                      smallest fiber dimension to try, default = 0
# verbose                       print intermediate results, default = false
# homogeneous                   use homogeneous ansatz for polynomial reconstruction, default = true
# -------------  Output:
# interpolants_per_dim          a list of lists of discriminants, organized per dimension of a generic fiber of the projection map πₚ
# samplegroups_per_dimension    a list of lists of samples used to interpolate these discriminants
# gaps                          a list of lists of gaps, indicating how accurately the equations are expected to be computed (larger gap = more accurate)
function project_codim1(LE,pars,vars; HP = false, first_fd = 0, verbose = false, homogeneous = true)
    groups_per_dim = []
    degs_per_dim = []
    parametrized_systems = []
    lines = []
    fibercuts = []
    n_sing = 1
    if verbose
        println("Scanning over the possible fiber dimensions...")
    end
    fd = first_fd # fiber dimension
    while (n_sing > 0) && (fd - 1 < length(vars))
        groups, incidence_param, line, fibercut, degs, n_sing = assess(LE,pars,vars; fiberdim = fd, verbose = verbose)
        push!(groups_per_dim, groups)
        push!(parametrized_systems, incidence_param)
        push!(lines, line)
        push!(fibercuts, fibercut)
        push!(degs_per_dim,degs)
        if verbose
            println("   found $(length(groups)) group(s) of solutions with fiber dimension $fd")
        end
        fd += 1
    end
    if verbose
        println("Sampling all identified components...")
    end
    samplegroups_per_dimension = []
    foundpositive_per_dimension = []
    for i = 1:length(groups_per_dim)
        incidence_param = parametrized_systems[i]
        line = lines[i]
        fibercut = fibercuts[i]
        push!(samplegroups_per_dimension,[])
        push!(foundpositive_per_dimension,[])
        for j = 1:length(groups_per_dim[i])
            deg = degs_per_dim[i][j]
            foundpositive = false
            if verbose
                println("Component $j with fiber dimension $(i-1) has estimated degree $deg")
            end
            if homogeneous
                nsamples = Int64(ceil(binomial(length(pars) + deg - 1, deg)*1.1))
            else
                nsamples = Int64(ceil(binomial(length(pars) + deg, deg)*1.1))
            end
            start_samples_upstairs = groups_per_dim[i][j]

            if HP
                vrs = [pars;vars]
                eqns = subs(expressions(incidence_param),parameters(incidence_param)=>[line[:];fibercut[:]])
                if length(vrs) < length(eqns)
                    eqns = randn(ComplexF64,length(vrs),length(eqns))*eqns
                end

                specialized_system = System(eqns,variables = vrs)
                start_samples_upstairs = [Landau.refine(specialized_system, x₀) for x₀ in start_samples_upstairs]

                start_samples_upstairs = [
                    convert.(BigFloat, real(sol)) + im * convert.(BigFloat, imag(sol))
                    for sol ∈ start_samples_upstairs
                    ]
            end

            if verbose && !isempty(findall([all(a -> isapprox(imag(a), 0.0, atol=1e-6) && real(a) > 0, s[length(pars)+1:end-1]/s[length(pars)+1]) for s in start_samples_upstairs]))
                foundpositive = true
            end

            start_samples = unique_points([sol[1:length(pars)] for sol in start_samples_upstairs], rtol = 1e-6)
            samples = start_samples
            if verbose
                println("Number of samples needed: $(nsamples)")
            end
            if nsamples < 10000
                ctr = 0
                while length(samples) < nsamples && ctr < 2*nsamples
                    newline = randn(ComplexF64,size(line))
                    newfibercut = randn(ComplexF64,size(fibercut))

                    local sols
                    try
                        if isempty([line[:];fibercut[:]])
                            R = HomotopyContinuation.solve(incidence_param, start_samples_upstairs, compile = false, threading = true, show_progress = verbose)
                        else
                            R = HomotopyContinuation.solve(incidence_param, start_samples_upstairs, start_parameters = [line[:];fibercut[:]], target_parameters = [newline[:];newfibercut[:]], compile = false, threading = true, show_progress = verbose)
                        end
                        sols = solutions(R)
                    catch e
                        if verbose
                            printstyled("HomotopyContinuation.jl failed to solve the system. This face needs to be reexamined!\n"; color = :yellow)
                        end
                        sols = []
                    end

                    if length(sols) > 0
                        
                        if HP
                            vrs = [pars;vars]
                            eqns = subs(expressions(incidence_param),parameters(incidence_param)=>[newline[:];newfibercut[:]])
                            if length(vrs) < length(eqns)
                                eqns = randn(ComplexF64,length(vrs),length(eqns))*eqns
                            end

                            specialized_system = System(eqns,variables = vrs)

                            sols = [Landau.refine(specialized_system, x₀) for x₀ in sols]
                            sols = [
                                convert.(BigFloat, real(sol)) + im * convert.(BigFloat, imag(sol))
                                for sol ∈ sols
                                    ]
                            end

                        if verbose && !foundpositive && !isempty(findall([all(a -> isapprox(imag(a), 0.0, atol=1e-6) && real(a) > 0, s[length(pars)+1:end-1]/s[length(pars)+1]) for s in sols]))
                            foundpositive = true
                        end

                        newsamples = unique_points([sol[1:length(pars)] for sol in sols], rtol = 1e-6)
                        samples = vcat(samples,newsamples)
                    end
                    ctr +=1
                end
                if ctr == 2*nsamples
                    if verbose
                        println(" couldn't collect sufficiently many samples")
                    end
                    push!(samplegroups_per_dimension[i],[])
                    push!(found_positive[i],[])
                else
                    push!(samplegroups_per_dimension[i],samples)
                    push!(foundpositive_per_dimension[i],foundpositive)
                end
            else
                if verbose
                    printstyled("Skipping due to too many samples needed\n"; color = :magenta)
                end
            end
        end
    end
    interpolants_per_dim = []
    gaps = []
    for i = 1:length(groups_per_dim)
        push!(interpolants_per_dim,[])
        push!(gaps,[])
        if length(samplegroups_per_dimension[i]) > 0
            for j = 1:length(samplegroups_per_dimension[i])
                samples = samplegroups_per_dimension[i][j]
                if length(samples) > 0
                    pol, N, gap = Landau.interpolate_deg(samples, degs_per_dim[i][j], pars; homogeneous = homogeneous)
                else
                    pol = 1.0 + 0*pars[1]
                    gap = 1
                end
                push!(interpolants_per_dim[i], Landau.rat(pol, tolerance = HP ? 1e-16 : 1e-8))
                push!(gaps[i],gap)
                if verbose
                    @printf("Component %i with fiber dimension %i has interpolation gap = %.2e\n", j, i-1, gap)
                    if foundpositive_per_dimension[i][j]
                        printstyled("Found α-positive samples in component $(j) with fiber dimension $(i-1)\n"; color = :cyan)
                    end
                end
            end
        end
    end
    return interpolants_per_dim, samplegroups_per_dimension, gaps
end

# Computes a specialized principal A-determinant for the polynomial f(x;a) with parameters a and variables x.
# For each face of the newton polytope of f, it computes the specialized discriminant of the face equation.
# See getSpecializedDiscriminant above.
# -------------  Input:
# f                     the polynomial we want the principal determinant of
# pars                  parameters of f
# vars                  variables of f. f is assumed to be a polynomial in both pars and vars
# method                use symbolic or numerical computation (:sym or :num), default = :sym
# high_prec             use high precision, default = false,
# codim_start           we compute discriminants for faces with codimension ≦ codim_start, by default, it is set to the dimension of the Newton polytope of f
# face_start            start computing from the face_start'th face, default = 1
# single_face           terminate after computing the discriminant of the face specified by codim_start and face_start, default = false
# single_weight         terminate after computing the discriminant of the face specified by single_weight, default = nothing
# verbose               print intermediate output, default = false
# homogeneous           use homogeneous ansatz for polynomial reconstruction, default = true
# save_output           save output to a file or "" if output is not to be saved, default = ""
# -------------  Output:
# discriminants         a list of discriminants in Oscar format. The principal Landau determinant is the product of all these discriminants
function getSpecializedPAD(f, pars, vars; method = :sym, high_prec = false, codim_start = -1, face_start = 1, single_face = false, single_weight = nothing, verbose = false, homogeneous = true, save_output = "")

    P = newton_polytope(f)
    Σ = normal_fan(P)

    E2 = convert.(Array{Int64},collect(Oscar.AbstractAlgebra.exponent_vectors(f)))
    mons2 = collect(Oscar.monomials(f))
    coeff2 = collect(Oscar.coefficients(f))

    discriminants = []
    unique_discriminants = []

    min_codim = Oscar.ambient_dim(P) - Oscar.dim(P)
    fvector = Oscar.f_vector(P)

    if verbose
        println("f = $(f)");
        println("pars = [$(split(string(pars),'[')[2])");
        println("vars = [$(split(string(vars),'[')[2])");
        println("method = $(method)");
        println("high_prec = $(high_prec)");
        println("codim_start = $(codim_start)");
        println("face_start = $(face_start)");
        println("single_face = $(single_face)");
        println("single_weight = $(single_weight)");
        println("verbose = $(verbose)");
        println("homogeneous = $(homogeneous)");
        println("save_output = $(save_output)");
        println("f_fector = [$(split(string(fvector),'[')[2])\n");
    end

    if codim_start < 0
        codim_start = Oscar.ambient_dim(P)
    end

    if single_weight !== nothing
        single_face = true
        weights = getWeights(f)
        indices = [(i, j) for i in 1:length(weights) for j in findall(x -> x == single_weight, weights[i])]
        if length(indices) > 0
            codim_start = indices[1][1] - 1
            face_start = indices[1][2]
        else
            printstyled("The specified weight $(single_weight) was not found\n"; color = :yellow)
            return
        end
    end

    for codim = codim_start : -1 : (single_face ? codim_start : min_codim)

        # Cones in this codimension
        codim_cones = cones(Σ,codim)
        if codim > min_codim && isnothing(codim_cones)
            continue
        end
        
        # Number of faces in this codimension
        noFaces = codim > min_codim ? length(codim_cones) : 1

        if !single_face && isnothing(single_weight)
            printstyled("------- codim = $(codim), $(noFaces) faces\n"; color = :blue)
            if verbose println("") end
        end

        # Loop over all the faces in this codimension
        for i in (codim == codim_start ? face_start : 1) : (single_face ? face_start : noFaces)
            if codim > min_codim
                weight = sum(Oscar.rays(codim_cones[i]))
                vals = convert.(Rational{Int64}, [transpose(weight)*e for e in E2] )
                mininds = findall(weight -> weight == minimum(vals), vals)
                IF = transpose(coeff2[mininds])*mons2[mininds]
            else
                weight = zeros(Int,Oscar.ambient_dim(P))
                IF = f
            end

            if verbose
                @show IF
            end

            disc = getSpecializedDiscriminant(IF,vars,pars; method = method, high_prec = high_prec, codim = codim, verbose = verbose, homogeneous = homogeneous);
            len_before = length(unique_discriminants)
            push!(discriminants, disc)

            disc_string = [replace(s, "//" => "/") for s in string.(vcat(disc...))]

            append!(unique_discriminants, disc_string)
            unique_discriminants = sort(unique(unique_discriminants))

            printstyled("codim: $(codim), face: $(i)/$(noFaces), weights: [$(join(string.(vcat(weight...)),", "))], discriminant: $(join(disc_string,", "))\n"; color = :red)
            if !single_face && isnothing(single_weight) && verbose println("") end
            if save_output != ""
                open(save_output, "a") do file
                    print(file, "codim: $(codim), face: $(i)/$(noFaces), weights: [$(join(string.(vcat(weight...)),", "))], discriminant: $(join(disc_string,", "))\n")
                end
            end

            len_after = length(unique_discriminants)

            if !single_face && isnothing(single_weight) && len_after > len_before
                printstyled("New discriminants after codim $(codim), face $(i)/$(noFaces). The list is: $(join(unique_discriminants, ", "))\n"; color = :yellow)
                if verbose println("") end
            end

            flush(stdout)
        end

        if !single_face && isnothing(single_weight)
            printstyled("Unique discriminants after codim $(codim): $(join(unique_discriminants,", "))\n"; color = :green)
            if verbose println("") end
        end

        flush(stdout)

    end

    return(discriminants)
end

# For a diagram specified with edges, nodes, and list of internal and external masses, getUF returns
# the first and second Symanzik polynomials, U and F, in the Oscar format. Kinematic parameters are
# automatically named for n ≤ 7 according to the substitutions below. User can also specify a custom
# list of substitutions. The masses can be specified as either of the four options:
# 1) :zero means zero masses
# 2) :equal means all equal and non-zero masses
# 3) :generic means all non-equal and non-zero masses
# 4) a custom list of masses in the HomotopyContinuation.jl format
# -------------  Input:
# edges                 list of internal edges specified as pairs of vertices
# nodes                 list of vertices to which external legs are attached
# internal_masses       list of internal masses-squared in the order of appearance in edges (:zero, :generic, or a list), default = :zero
# external_masses       list of external masses-squared in the order of appearance in nodes (:zero, :generic, or a list), default = :zero
# substitute            optional list of substitutions for kinematic parameters, default = [] 
# -------------  Output:
# U_oscar               first Symanzik polynomial in the Oscar format
# F_oscar               second Symanzik polynomial in the Oscar format
# pars_oscar            the list of kinematic parameters in the Oscar format
# vars_oscar            the list of variables (Schwinger parameters) in the Oscar format
function getUF(edges, nodes; internal_masses = :zero, external_masses = :zero, substitute = [])

    F, U, αHC, pp, mm = Landau.getF(edges,nodes)
    n = length(nodes)
    E = length(edges)
    @var M[1:n]

    if n == 2
        subslist_pp = [[M[1], -M[1]]
                       [-M[1], M[1]]];
    elseif n == 3
        subslist_pp = [[M[1], (-M[1] - M[2] + M[3])/2, (-M[1] + M[2] - M[3])/2]
                       [(-M[1] - M[2] + M[3])/2, M[2], (M[1] - M[2] - M[3])/2]
                       [(-M[1] + M[2] - M[3])/2, (M[1] - M[2] - M[3])/2, M[3]]];
    elseif n == 4
        @var s, t;
        subslist_pp = [[M[1], (s - M[1] - M[2])/2, (-s - t + M[2] + M[4])/2, (t - M[1] - M[4])/2]
                       [(s - M[1] - M[2])/2, M[2], (t - M[2] - M[3])/2, (-s - t + M[1] + M[3])/2]
                       [(-s - t + M[2] + M[4])/2, (t - M[2] - M[3])/2, M[3], (s - M[3] - M[4])/2]
                       [(t - M[1] - M[4])/2, (-s - t + M[1] + M[3])/2, (s - M[3] - M[4])/2, M[4]]];
    elseif n == 5
        @var s12, s23, s34, s45, s51;
        subslist_pp = [[M[1], (s12 - M[1] - M[2])/2, (-s12 - s23 + s45 + M[2])/2, (s23 - s45 - s51 + M[5])/2, (s51 - M[1] - M[5])/2]
                       [(s12 - M[1] - M[2])/2, M[2], (s23 - M[2] - M[3])/2, (-s23 - s34 + s51 + M[3])/2, (-s12 + s34 - s51 + M[1])/2]
                       [(-s12 - s23 + s45 + M[2])/2, (s23 - M[2] - M[3])/2, M[3], (s34 - M[3] - M[4])/2, (s12 - s34 - s45 + M[4])/2]
                       [(s23 - s45 - s51 + M[5])/2, (-s23 - s34 + s51 + M[3])/2, (s34 - M[3] - M[4])/2, M[4], (s45 - M[4] - M[5])/2]
                       [(s51 - M[1] - M[5])/2, (-s12 + s34 - s51 + M[1])/2, (s12 - s34 - s45 + M[4])/2, (s45 - M[4] - M[5])/2, M[5]]];
    elseif n == 6
        @var s12, s23, s34, s45, s56, s61, s123, s234, s345;
        subslist_pp = [[M[1], (s12 - M[1] - M[2])/2, (-s12 + s123 - s23 + M[2])/2, (-s123 + s23 - s234 + s56)/2, (s234 - s56 - s61 + M[6])/2, (s61 - M[1] - M[6])/2]
                       [(s12 - M[1] - M[2])/2, M[2], (s23 - M[2] - M[3])/2, (-s23 + s234 - s34 + M[3])/2, (-s234 + s34 - s345 + s61)/2, (-s12 + s345 - s61 + M[1])/2]
                       [(-s12 + s123 - s23 + M[2])/2, (s23 - M[2] - M[3])/2, M[3], (s34 - M[3] - M[4])/2, (-s34 + s345 - s45 + M[4])/2, (s12 - s123 - s345 + s45)/2]
                       [(-s123 + s23 - s234 + s56)/2, (-s23 + s234 - s34 + M[3])/2, (s34 - M[3] - M[4])/2, M[4], (s45 - M[4] - M[5])/2, (s123 - s45 - s56 + M[5])/2]
                       [(s234 - s56 - s61 + M[6])/2, (-s234 + s34 - s345 + s61)/2, (-s34 + s345 - s45 + M[4])/2, (s45 - M[4] - M[5])/2, M[5], (s56 - M[5] - M[6])/2]
                       [(s61 - M[1] - M[6])/2, (-s12 + s345 - s61 + M[1])/2, (s12 - s123 - s345 + s45)/2, (s123 - s45 - s56 + M[5])/2, (s56 - M[5] - M[6])/2, M[6]]];
    elseif n == 7
        @var s12, s23, s34, s45, s56, s67, s71, s123, s234, s345, s456, s567, s671, s712;
        subslist_pp = [[M[1], (s12 - M[1] - M[2])/2, (-s12 + s123 - s23 + M[2])/2, (-s123 + s23 - s234 + s567)/2, (s234 - s567 + s67 - s671)/2, (-s67 + s671 - s71 + M[7])/2, (s71 - M[1] - M[7])/2]
                       [(s12 - M[1] - M[2])/2, M[2], (s23 - M[2] - M[3])/2, (-s23 + s234 - s34 + M[3])/2, (-s234 + s34 - s345 + s671)/2, (s345 - s671 + s71 - s712)/2, (-s12 - s71 + s712 + M[1])/2]
                       [(-s12 + s123 - s23 + M[2])/2, (s23 - M[2] - M[3])/2, M[3], (s34 - M[3] - M[4])/2, (-s34 + s345 - s45 + M[4])/2, (-s345 + s45 - s456 + s712)/2, (s12 - s123 + s456 - s712)/2]
                       [(-s123 + s23 - s234 + s567)/2, (-s23 + s234 - s34 + M[3])/2, (s34 - M[3] - M[4])/2, M[4], (s45 - M[4] - M[5])/2, (-s45 + s456 - s56 + M[5])/2, (s123 - s456 + s56 - s567)/2]
                       [(s234 - s567 + s67 - s671)/2, (-s234 + s34 - s345 + s671)/2, (-s34 + s345 - s45 + M[4])/2, (s45 - M[4] - M[5])/2, M[5], (s56 - M[5] - M[6])/2, (-s56 + s567 - s67 + M[6])/2]
                       [(-s67 + s671 - s71 + M[7])/2, (s345 - s671 + s71 - s712)/2, (-s345 + s45 - s456 + s712)/2, (-s45 + s456 - s56 + M[5])/2, (s56 - M[5] - M[6])/2, M[6], (s67 - M[6] - M[7])/2]
                       [(s71 - M[1] - M[7])/2, (-s12 - s71 + s712 + M[1])/2, (s12 - s123 + s456 - s712)/2, (s123 - s456 + s56 - s567)/2, (-s56 + s567 - s67 + M[6])/2, (s67 - M[6] - M[7])/2, M[7]]];
    else
        println("Unsupported number of external legs in the automatic mode. Try running getSpecializedPAD manually.")
        return [], []
    end

    Fsubs = subs(F, pp[:] => subslist_pp)

    if internal_masses == :zero
        Fsubs = subs(Fsubs, mm => zeros(Int64, E))
    elseif internal_masses == :generic
        @var m[1:E]
        Fsubs = subs(Fsubs, mm => m)
    elseif internal_masses == :equal
        @var m2
        Fsubs = subs(Fsubs, mm => [m2 for i in 1:length(mm)]) 
    else
        Fsubs = subs(Fsubs, mm => internal_masses)
    end

    if external_masses == :zero
        Fsubs = subs(Fsubs, M => zeros(Int64, n))
    elseif external_masses == :generic
    elseif external_masses == :equal
        @var M2
        Fsubs = subs(Fsubs, M => [M2 for i in 1:length(M)])
    else
        Fsubs = subs(Fsubs, M => external_masses)
    end

    if !isempty(substitute)
        for subrule in substitute
            Fsubs = subs(Fsubs, subrule)
        end
    end

    Fsubs = HomotopyContinuation.ModelKit.expand(Fsubs)

    vars = αHC
    pars = setdiff(variables(Fsubs),vars)

    U_oscar, S, pars_oscar, vars_oscar = HC_to_oscar_S_mynames(U, pars, vars; parnames = string.(pars), varnames = string.(vars))
    F_oscar, S, pars_oscar, vars_oscar = HC_to_oscar_S_mynames(Fsubs, pars, vars; parnames = string.(pars), varnames = string.(vars))

    return U_oscar, F_oscar, pars_oscar, vars_oscar

end

# Converts numerical subscripts in a string to arguments in square brackets, for example "m₁" => "m[1]" or "M₁₂" => "M[12]"
# -------------  Input:
# str                 input string
# -------------  Output:
# output string
function subscript_to_bracket(str)
    for (sub, num) in zip('₀':'₉', '0':'9')
        str = replace(str, string(sub) => "["*string(num)*"]")
    end
    return replace(str, "][" => "")
end

# A simplified interface for the function getSpecializedPAD.
# As an input, it takes a Feynman diagram specified by the lists of internal edges, external legs,
# and the internal and external masses in the following way. We label the vertices with consecutive
# integers 1,2,…,V. The list of edges is specified as a length-E list of pairs of vertices that
# are connected (multi-edges are allowed). The list of nodes is a length-n list of vertices to
# which external legs are attached (multiple legs can be attached to the same vertex). They are
# assigned the external momenta pᵢ for i=1,2,…,n ≤ 7 in the order of appearance.
# The lists of internal and external masses have lengths E and n respectively and take either
# symbolic or numerical values (if method = :num and numerical values are provided, don't forget to
# turn on homogeneous = false). The function getPLD automatically assigns the Mandelstam invariants
# sᵢⱼ… = (pᵢ + pⱼ + ⋯)² in the cyclic basis, see below. In the special case when n ≤ 3, there are
# no Mandelstam invariants, and for n = 4, it uses the convention s = s₁₂, t = s₂₃.
# The remaining options are the same as for getSpecializedPAD, see above.
#
# The option load_output allows the user to load the text output from a previous run instead of
# computing the discriminants from scratch. It scrapes all lines of the form:
#
# codim: 6, face: 22/153, weights: [-1, 0, 0, -1, 0, 0, -2], discriminant: t, m2 - 1/4*t
#
# from the file load_output and adds the discriminats to the running list in Oscar format.
#
# Example usage for a double-box diagram with massless external legs and a massive outer loop
# (all internal edges have masses-squared equal to m2, except for the middle rung which is massless):
#
#   p₂ ---2----3----4--- p₃
#         |    |    |
#         |    |    |
#   p₁ ---1----6----5--- p₄
#
# @var m2
# edges = [[1,2],[2,3],[3,4],[4,5],[5,6],[6,1],[3,6]];
# nodes = [1,2,4,5];
# internal_masses = [m2,m2,m2,m2,m2,m2,0];
# external_masses = [0,0,0,0];
#
# getPLD(edges, nodes, internal_masses = internal_masses, external_masses = external_masses, method = :num)
#
# -------------  Input:
# edges                 list of internal edges specified as pairs of vertices
# nodes                 list of vertices to which external legs are attached
# internal_masses       list of internal masses-squared in the order of appearance in edges (:zero, :generic, or a list), default = :zero
# external_masses       list of external masses-squared in the order of appearance in nodes (:zero, :generic, or a list), default = :zero
# method                use symbolic or numerical computation (:sym or :num), default = :sym
# high_prec             use high precision, default = false,
# codim_start           we compute discriminants for faces with codimension ≦ codim_start, by default, it is set to the dimension of the Newton polytope of f
# face_start            start computing from the face_start'th face, default = 1
# single_face           terminate after computing the discriminant of the face specified by codim_start and face_start, default = false
# single_weight         terminate after computing the discriminant of the face with the specified weight, default = nothing
# verbose               print intermediate output, default = false
# homogeneous           use homogeneous ansatz for polynomial reconstruction, default = true
# save_output           save output to a file, default = ""
# load_output           load output from a previous run instead of computing discriminants, default = ""
# substitutions         extra substitutions for the external kinematics (use this to study certain kinematic limits)
# -------------  Output:
# discriminants         a list of discriminants in Oscar format. The principal Landau determinant is the product of all these discriminants
# pars_oscar            a list of kinematic variables in Oscar format: Mandelstam invariants set internally and the possible masses provided by the user
# vars_oscar            a list of Schwinger parameters in Oscar format
# U_oscar               the first Symanzik polynomial in Oscar format
# F_oscar               the second Symanzik polynomial in Oscar format
function getPLD(edges, nodes; internal_masses = :zero, external_masses = :zero, method = :sym, high_prec = false, codim_start = -1, face_start = 1, single_weight = nothing, single_face = false, verbose = false, homogeneous = true, save_output = "", load_output = "", substitutions = [])

    U_oscar, F_oscar, pars_oscar, vars_oscar = getUF(edges, nodes; internal_masses = internal_masses, external_masses = external_masses, substitute = substitutions)


    if load_output == ""
        discriminants = getSpecializedPAD(U_oscar + F_oscar, pars_oscar, vars_oscar; method = method, high_prec = high_prec, codim_start = codim_start, face_start = face_start, single_face = single_face, single_weight = single_weight, verbose = verbose, homogeneous = homogeneous, save_output = save_output)
    else
        pars_string = [subscript_to_bracket(str) for str in string.(pars_oscar)]

        if internal_masses == :generic
            global m = Vector{fmpq_mpoly}(undef, length(edges))
        end
        if external_masses == :generic
            global M = Vector{fmpq_mpoly}(undef, length(nodes))
        end

        global R, pars = polynomial_ring(QQ, pars_string)
        eval(Meta.parse("global ("*join(pars_string, ", ")*") = pars"))

        discriminants = []
        open(load_output, "r") do f
            for line in eachline(f)
                # Extract discriminants using regex
                disc = match(r"discriminant: (.+)$", line)
                if disc != nothing
                    push!(discriminants, eval(Meta.parse(replace("["*subscript_to_bracket(disc[1])*"]", "/" => "//"))))
                end
            end
        end
    end

    return discriminants, pars_oscar, vars_oscar, U_oscar, F_oscar

end

# For a given polynomial f, computes the signed Euler characteristic |χ(X)| of X = (C^*)^{|vars|} \ {f(pᵣ)=0}
# for a random kinematic point pᵣ on {candidate=0}. Here, candidate is a candidate component of the Euler discriminant written
# as a polynomial in pars. The function uses homotopy continuation methods to compute the Euler characteristic
# as a number of solutions of the log-likelihood equations. Alternatively, setting candidate = :random, the
# function computes the Euler characteristic for generic point in the kinematic space. 
# -------------  Input:
# f                     the polynomial f
# pars                  parameters of f
# vars                  variables of f. f is assumed to be a polynomial in both pars and vars
# candidate             candidate component of the Euler discriminant, default = :random
# homogeneous           is f a homogeneous polynomial in vars?, default = false
# cert                  whether to certify the result of homotopy continuation, default = true
# verbose               print intermediate output, default = false
# -------------  Output:
# signed Euler characteristic
function getGenericEuler(f, pars, vars, candidate = :random; homogeneous = false, cert = true, verbose = false, allowed_retries = 5)

    parsstring = prod([string(p)*" " for p in pars])
    HCpars = [eval(Meta.parse("@var "*parsstring))...]
    varsstring = prod([string(p)*" " for p in vars])
    HCvars = [eval(Meta.parse("@var "*varsstring))...]

    fHC = oscar_to_HC_S(f,HCvars,HCpars)

    newcoeffs = []
    retry_count = 0

    while is_empty(newcoeffs)
        if candidate == :random # Random point
            newcoeffs = randn(ComplexF64,length(pars))
        else # Random point on candidate = 0
            candidateHC = oscar_to_HC_Q(candidate,HCpars)
            if verbose
                println("Checking candidate ", candidateHC)
            end
            if length(pars) == 1
                # Due to a bug in HomotopyContinuation.jl, we cannot solve the univariate equation x = 0
                # Here, we implement a quick and dirty fix
                @var z
                R = HomotopyContinuation.solve(System([candidateHC;z-1], variables = [HCpars;z]), compile = false, threading = true, show_progress = verbose)
                newcoeffs = [solutions(R)[1][1]]
            else
                A = randn(ComplexF64,length(pars)-1,length(pars))
                b = randn(ComplexF64,length(pars)-1)
                lineq = A*HCpars - b

                S = System([lineq;candidateHC], variables = HCpars)
                R = HomotopyContinuation.solve(S; compile = false, threading = true, show_progress = verbose)
                
                if is_empty(solutions(R))
                    println("No solutions found for the Euler characteristic. Retrying...")
                    retry_count += 1
                    if retry_count > allowed_retries
                        println("PLD.jl couldn't get the Euler characteristic!")
                        break
                    end
                    continue
                else
                    newcoeffs = solutions(R)[1]
                end
            end
        end
    end

    if homogeneous # Dehomogenize if necessary
        fHC = subs(fHC, HCvars[1]=>1)
        HCvars = HCvars[2:end]
    end

    # Likelihood equations
    @var μ, ν[1:length(vars)]
    L =  μ*log(fHC) + sum([ν[i]*log(HCvars[i]) for i = 1:length(vars)])

    S = System(subs(differentiate(L,[HCvars...]), HCpars => newcoeffs), variables = HCvars, parameters = [μ;ν])
    sol = monodromy_solve(S; compile = false, threading = true, show_progress = false)

    # Possibly certify
    if (cert == true) && (length(pars) > 1)
        cert_result = certify(S, sol)
        return ndistinct_certified(cert_result)
    else
        sol_distinct = unique_points(solutions(sol), rtol=1e-6)
        return length(sol_distinct)
    end

end

# A wrapper for getEuler that takes as input a Feynman diagram specified by the lists of edges, nodes, internal_masses,
# external_masses, and returns the Euler characteristic.
# -------------  Input:
# edges                 list of internal edges specified as pairs of vertices
# nodes                 list of vertices to which external legs are attached
# internal_masses       list of internal masses-squared in the order of appearance in edges (:zero, :generic, or a list), default = :zero
# external_masses       list of external masses-squared in the order of appearance in nodes (:zero, :generic, or a list), default = :zero
# substitute            optional list of substitutions for kinematic parameters, default = [] 
# -------------  Output:
# signed Euler characteristic
function getEuler(edges, nodes; internal_masses = :zero, external_masses = :zero, substitute = [])
    U, F, pars, vars = getUF(edges,nodes; internal_masses = internal_masses, external_masses = external_masses, substitute = substitute)
    return getGenericEuler(U+F,pars,vars)
end

# Checks whether the provided candidate components belong to the Euler discriminant by computing their signed Euler
# characteristics χ and comparing them to the generic signed Euler characteristic χ∗. The function returns a list
# of those candidates for which χ < χ∗. To increase reliability of the results, each Euler characteristic is computed
# multiple times and we select the maximum of the computed values.
# -------------  Input:
# f                     the polynomial f
# pars                  parameters of f
# vars                  variables of f. f is assumed to be a polynomial in both pars and vars
# candidates            list of polynomials specifying candidate components of the Euler discriminant
# repeats               number of times to compute the Euler characteristic for each candidate, default = 10
# skip_1                do we skip the empty discriminants "1"?, default = true
# -------------  Output:
# discs                 list of candidate components of the Euler discriminant for which χ < χ∗
# eulerchars            list of signed Euler characteristics χ for the candidate components in discs
function EulerDiscriminantQ(f, pars, vars, candidates; repeats = 10, skip_1 = true)

    genericEuler = maximum([getGenericEuler(f, pars, vars, :random) for i in 1:repeats])
    println("Generic |Euler characteristic|, χ∗ = $(genericEuler)")

    if skip_1
        idx = findall(d -> string(d) != "1" , candidates)
        candidates = candidates[idx]
    end

    @show candidates

    discs = []
    eulerchars = []
    for cand in candidates
        candEuler = maximum([getGenericEuler(f, pars, vars, cand) for i in 1:repeats])
        if candEuler < genericEuler
            printstyled("Subspace $(cand) has χ = $(candEuler) < χ∗\n"; color = :green)
            push!(discs, cand)
            push!(eulerchars,candEuler)
        else
            printstyled("Subspace $(cand) has χ = $(candEuler) = χ∗\n"; color = :red)
        end
    end
    return discs, eulerchars
end

# Takes the list of discriminants from the output of getSpecializedPAD or getPLD for a given polynomial f
# and returns the list of unique discriminants with their correcsponding face weights.
# -------------  Input:
# f                     the polynomial f
# discs                 list of discriminants from the output of getSpecializedPAD or getPLD
# skip_1                do we skip the empty discriminants "1"?, default = true
# verbose               print intermediate output, default = true
# -------------  Output:
# unique_discs          list of unique discriminants
# weight_list           list of weights for each unique discriminant
function discriminants_with_weights(f, discs; faces_done = [], skip_1 = true, verbose = true)

    all_wghts = getWeights(f)

    if isempty(faces_done)
        wghts = vcat(reverse(all_wghts)...)
    else
        wghts = []
        for i in 1:length(faces_done)
            push!(wghts, all_wghts[faces_done[i][1]+1][faces_done[i][2]])
        end
    end


    flat_discs = reduce(vcat, discs)
    unique_discs = unique(flat_discs)
    unique_discs_string = string.(unique_discs)

    sort = sortperm(unique_discs_string)
    unique_discs = unique_discs[sort]
    unique_discs_string = unique_discs_string[sort]

    if skip_1
        idx = findall(d -> d != "1" , unique_discs_string)
        unique_discs = unique_discs[idx]
        unique_discs_string = unique_discs_string[idx]
    end

    Σ = Dict(unique_discs.=>1:length(unique_discs))
    weight_list = [[] for i = 1:length(unique_discs)]

    for i = 1:length(discs)
        for j = 1:length(discs[i])
            if (string(discs[i][j]) in unique_discs_string)
                push!(weight_list[Σ[discs[i][j]]], wghts[i])
            end
        end
    end

    if verbose
        println("Unique discriminants with weights:")
        for i = 1:length(unique_discs)
            println(string(unique_discs[i])*" => "*string(weight_list[i]))
        end
    end
    
    return unique_discs, weight_list

end

end