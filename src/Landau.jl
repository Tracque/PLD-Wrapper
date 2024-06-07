module Landau

# Authors: Sebastian Mizera and Simon Telen
# Date: October 16, 2023
# Short description: This code accompanies the paper
# [LD] Sebastian Mizera and Simon Telen, "Landau Discriminants", arXiv:2109.08036

using HomotopyContinuation
using LinearAlgebra
using GenericSVD
using Arpack
include("refine.jl")

export getU,
    getF,
    LandauEquations,
    affineLandauEquations,
    sampleProjection,
    degreeProjection,
    getVdm,
    interpolate_deg,
    interpolate_mons,
    rat
export substitute4legs,
    substitute5legs,
    sampleProjection_HP,
    decompose,
    getInterpolantSparse,
    getLandauDiscriminant,
    distinct_inds

# Computes the first Symanzik polynomial of a diagram.
# -------------  Input:
# edges     a list of edges of the diagram
# -------------  Output:
# U         the first Symanzik polynomial
# M         an auxiliary matrix M used in computing U
# α         an array of Schwinger parameters
# vtcs      a list of vertex labels for the diagram
function getU(edges)
    vtcs = unique(vcat(edges...))
    unique_edges = unique(edges)
    VG = length(unique(vcat(edges...)))
    I = Dict(vtcs .=> 1:VG)
    @var α[1:length(edges)]
    M = fill(zero(typeof((sum(α)))), VG, VG)
    J = Dict(edges .=> 1:length(edges))
    for i = 1:length(edges)6
        e = edges[i]
        M[I[e[1]], I[e[2]]] += -α[i]
        M[I[e[2]], I[e[1]]] += -α[i]
        M[I[e[1]], I[e[1]]] += α[i]
        M[I[e[2]], I[e[2]]] += α[i]
    end
    U0 = det(M[1:VG-1, 1:VG-1])
    E, C = exponents_coefficients(U0, α)
    E = -(E .- 1)
    U = sum([prod(α .^ E[:, i]) for i = 1:length(C)])
    return U, M, α, vtcs
end

# Computes the second Symanzik polynomial of a diagram.
# -------------  Input:
# edges     a list of edges of the diagram
# nodes     a list of nodes of the diagram to which external legs are attached
# -------------  Output:
# F         the second Symanzik polynomial
# U         the first Symanzik polynomial
# α         an array of Schwinger parameters
# p         an array of momentum vectors
# mm        an array of (squared) internal mass parameters
function getF(edges, nodes)
    U, M, α, vtcs = getU(edges)
    VG = length(vtcs)
    NG = length(nodes)
    I = Dict(vtcs .=> 1:VG)
    EG = length(edges)
    @var p[1:NG, 1:NG] mm[1:EG]
    F = 0
    for i = 1:NG
        for j = i+1:NG
            if nodes[i] == nodes[j]
                continue
            end
            if size(M, 1) == 2
                g = 1 + 0.0 * sum(α)
            else
                g = det(M[
                    setdiff(1:VG, I[nodes[i]], I[nodes[j]]),
                    setdiff(1:VG, I[nodes[i]], I[nodes[j]]),
                ])
            end
            E, C = exponents_coefficients(g, α)
            E = -(E .- 1)
            g = sum([prod(α .^ E[:, i]) for i = 1:length(C)])
            F = F - p[i, j] * g
        end
    end
    F = F - U * sum([α[i] * mm[i] for i = 1:EG])
    return F, U, α, p, mm
end

# Computes the Landau equations of a diagram.
# -------------  Input:
# edges     a list of edges of the diagram
# nodes     a list of nodes of the diagram to which external legs are attached
# -------------  Output:
# LE        an array of polynomials, LE = 0 are the Landau equations (without the condition U ≠ 0)
# α         an array of Schwinger parameters
# p         an array of momentum vectors
# mm        an array of (squared) internal mass parameters
function LandauEquations(edges, nodes)
    F, U, α, p, mm = getF(edges, nodes)
    LE = differentiate(F, α)
    LE, α, p, mm
end

# Computes the affine Landau equations of a diagram.
# -------------  Input:
# edges     a list of edges of the diagram
# nodes     a list of nodes of the diagram to which external legs are attached
# -------------  Output:
# ALE       an array of polynomials, ALE = 0 are the affine Landau equations from Section 3.2 in [LD]
# y         an extra variable to incorporate U ≂̸ 0 in the equations
# α         an array of Schwinger parameters
# p         an array of momentum vectors
# mm        an array of (squared) internal mass parameters
function affineLandauEquations(edges, nodes)
    F, U, α, p, mm = getF(edges, nodes)
    LE = differentiate(F, α)
    @var y
    ALE =
        [evaluate(LE, α[end] => 1.0); prod(α[1:end-1]) * evaluate(U, α[end] => 1.0) * y - 1]
    return ALE, y, α, p, mm
end

# Sample the projection of a variety given by polynomial equations on ℂᵐ × ℂᵉ onto ℂᵉ.
# -------------  Input:
# F         an array of polynomials giving the polynomial system in variables x and p.
# x and p   two disjoint groups of coordinates, x on ℂᵐ, p on ℂᵉ.
# -------------  Optional input:
# The system F(x,p) is augmented with linear equations A*p + B = 0. If a solution (x,p) is known for some A, B,
# it can be passed to this function via seedA, seedB, seedsol. Other solutions will then be computed via monodromy on (A,B).
# codimen           expected codimension of the projection
# distinct_tol      tolerance for deciding whether two projected points are the same
# npoints           a lower bound on the number of samples that need to be computed
# printprogress     boolean variable, set to true if you want information to be printed throughout the execution
# findSingular      boolean variable, set to true if you want to allow singular points on the incidence variety to produce samples
# denseStart        boolean variable, set to true if you want to use a total degree start system instead of the default polyhedral start system
# -------------  Output:
# samples       sample points on the projection
# monres        result of the initial solving step using HomotopyContinuation.jl
# AA, BB        parameters used for the linear slice AA*p + BB = 0 in the initial solving step
# H             augmented system of equations (F(x,p);AA*p+BB) = 0
# singsamples   sample points produced by singular points on the incidence variety (empty if findSingular = false)
function sampleProjection(
    F,
    x,
    p;
    seedA = nothing,
    seedB = nothing,
    seedsol = nothing,
    codimen = 1,
    distinct_tol = 1e-6,
    npoints = 100,
    printprogress = false,
    findSingular = false,
    denseStart = false,
)
    e = length(p)
    s = e - codimen
    @var A[1:s, 1:e] b[1:s]
    L = A * p + b
    if codimen > 1 # fibers are positive dimensional, add generic linear equations in all variables
        A2 = randn(ComplexF64, codimen - 1, length(x))
        b2 = randn(ComplexF64, codimen - 1)
        L2 = A2 * x + b2
        H = System([F; L; L2], variables = [p; x], parameters = [A[:]; b])
    else
        H = System([F; L], variables = [p; x], parameters = [A[:]; b])
    end
    if isnothing(seedsol)
        targetpars = randn(ComplexF64, length(parameters(H)))
        if denseStart
            monres = solve(H, target_parameters = targetpars; start_system = :total_degree)
        else
            monres = solve(H, target_parameters = targetpars)
        end
        if length(solutions(monres)) == 0
            println("found no starting solutions.")
        end
    else
        monres = monodromy_solve(H, seedsol, [seedA[:]; seedB])
        targetpars = parameters(monres)
    end
    psols = [sol[1:e] for sol ∈ solutions(monres)]
    if length(psols) > 0
        samples, inds = distinct_inds(psols, tol = distinct_tol)
    else
        samples = []
    end
    startsamples = samples
    if findSingular && isnothing(seedsol) && (length(singular(monres)) > 0)
        println("found some singular solutions")
        singsols = [ss.solution for ss ∈ singular(monres)]
        psols = [sol[1:e] for sol ∈ singsols]
        singsamples, inds = distinct_inds(psols, tol = distinct_tol)
        dsing = length(singsamples)
        println("singular component is estimated to have degree >= $dsing")
    else
        singsamples = []
    end
    if length(samples) + length(singsamples) == 0
        println("Found no initial samples, are you sure about the codimension?")
    else
        while length(samples) + length(singsamples) < npoints
            newparams = randn(ComplexF64, length(targetpars))
            if findSingular && length(singsamples) > 0
                R = solve(H, target_parameters = newparams)
                psols = [sol[1:e] for sol ∈ solutions(R)]
                if length(psols) > 0
                    newsamples, inds = distinct_inds(psols, tol = distinct_tol)
                else
                    newsamples = []
                end
                samples = [samples; newsamples]
                singsols = [ss.solution for ss ∈ singular(R)]
                psols = [sol[1:e] for sol ∈ singsols]
                newsingsamples, inds = distinct_inds(psols, tol = distinct_tol)
                singsamples = [singsamples; newsingsamples]
            else
                R = solve(
                    H,
                    solutions(monres);
                    start_parameters = targetpars,
                    target_parameters = newparams,
                )
                psols = [sol[1:e] for sol ∈ solutions(R)]
                newsamples, inds = distinct_inds(psols, tol = distinct_tol)
                samples = [samples; newsamples]
            end
            if printprogress
                println("Collected $(length(samples)) samples")
                if findSingular
                    println("singlar: $(length(singsamples))")
                end
            end
        end
    end
    AA = targetpars[1:length(A[:])]
    BB = targetpars[length(A[:])+1:end]
    return samples, monres, AA, BB, H, singsamples
end

# Compute the degree of the projection of a variety given by polynomial equations on ℂᵐ × ℂᵉ onto ℂᵉ.
# -------------  Input:
# F         an array of polynomials giving the polynomial system in variables x and p.
# x and p   two disjoint groups of coordinates, x on ℂᵐ, p on ℂᵉ.
# -------------  Optional input:
# The system F(x,p) is augmented with linear equations A*p + B = 0. If a solution (x,p) is known for some A, B,
# it can be passed to this function via seedA, seedB, seedsol. Other solutions will then be computed via monodromy on (A,B).
# codimen           expected codimension of the projection
# distinct_tol      tolerance for deciding whether two projected points are the same
# findSingular      boolean variable, set to true if you want to allow singular points on the incidence variety to produce samples
# -------------  Output:
# The degree of the projection
function degreeProjection(
    F,
    x,
    p;
    codimen = 1,
    distinct_tol = 1e-6,
    findSingular = false,
    seedA = nothing,
    seedB = nothing,
    seedsol = nothing,
)
    samples, monres, AA, bb, H, singsamples = sampleProjection(
        F,
        x,
        p;
        npoints = 1,
        codimen = codimen,
        distinct_tol = distinct_tol,
        findSingular = findSingular,
        seedA = seedA,
        seedB = seedB,
        seedsol = seedsol,
    )
    length(samples) + length(singsamples)
end

# Computes all monomials of (or up to) a given degree
# -------------  Input:
# varlist   a list of variables for the monomials
# degrees   degree of the required monomials
# -------------  Optional input:
# homogeneous       boolean variable, set to true if you want all monomials of the fixed degree,
#                   set to false if you want all monomials up to the fixed degree
# -------------  Output:
# mons      a list of the requested monomials
function getMonomials(varlist, degree; homogeneous = true)
    if homogeneous
        E, C = exponents_coefficients(sum(varlist)^degree, varlist)
        mons = [prod(varlist .^ E[:, i]) for i = 1:size(E, 2)]
    else
        E, C = exponents_coefficients(sum([1; varlist])^degree, varlist)
        mons = [prod(varlist .^ E[:, i]) for i = 1:size(E, 2)]
    end
    return mons
end

# Computes a Vandermonde matrix of a set of monomials evaluated in a set of points
# -------------  Input:
# samples   a list of sample points that you want to interpolate
# mons      a set of monomials that are allowed to occur in the interpolant
# varlist   a list of variables
# -------------  Optional input:
# scaling   boolean variable, set to true if you want the rows of the Vandermonde matrix to be scaled to have norm 1
# -------------  Output:
# vdm      a Vandermonde matrix
# vdmvecs  the rows of vdm
function getVdm(samples, mons, varlist; scaling = true)
    exps, C = exponents_coefficients(sum(mons), varlist)
    I = Dict(mons .=> 1:length(mons))
    newmons = prod(varlist .^ (exps), dims = 1)
    sortinds = sortperm([I[m] for m ∈ newmons][:]) # make sure the monomials are used in the correct order.
    exps = exps[:, sortinds]
    vdmvecs = [prod(smpl .^ (exps), dims = 1) for smpl ∈ samples]
    if scaling
        vdmvecs = [vec / norm(convert.(ComplexF64, vec)) for vec ∈ vdmvecs]
    end
    vdm = vcat(vdmvecs...)
    return vdm, vdmvecs
end

# Interpolate a set of points by a polynomial of a fixed degree
# -------------  Input:
# samples   a list of sample points that you want to interpolate
# degree    degree of the interpolant
# varlist   a list of variables
# -------------  Optional input:
# homogeneous       boolean variable, set to true if you want all monomials of the fixed degree,
#                   set to false if you want all monomials up to the fixed degree
# -------------  Output:
# pol      an interpolant
# N        coefficients of the interpolant
# gap      ratio between the second to last and the last singular value of the interpolation problem
function interpolate_deg(samples, degree::Int, varlist; homogeneous = true)
    mons = getMonomials(varlist, degree; homogeneous = homogeneous)
    V, vdmvecs = getVdm(samples, mons, varlist)
    if size(V,2) == 1
        pol = 1.0*mons[1]
        N = [1.0]
        gap = Inf
    else
        svdobj = svd(V)
        N = svdobj.V[:, end]
        gap = svdobj.S[end-1] / svdobj.S[end]
        nonzind = argmax(abs.(N))
        N = N / N[nonzind]
        pol = sum([N[i] * mons[i] for i = 1:length(mons)])
    end
    return pol, N, gap
end

# Interpolate a set of points by a polynomial with predescribed support
# -------------  Input:
# samples   a list of sample points that you want to interpolate
# mons      monomials in the support of the interpolant
# varlist   a list of variables
# -------------  Output:
# pol      an interpolant
# N        coefficients of the interpolant
# gap      ratio between the second to last and the last singular value of the interpolation problem
function interpolate_mons(samples, mons::Array{Expression}, varlist)
    V, vdmvecs = getVdm(samples, mons, varlist)
    svdobj = svd(V)
    N = svdobj.V[:, end]
    gap = svdobj.S[end-1] / svdobj.S[end]
    nonzind = argmax(abs.(N))
    N = N / N[nonzind]
    pol = sum([N[i] * mons[i] for i = 1:length(mons)])
    return pol, N, gap
end

# Approximate a polynomial with floating point coefficients by one with rational coefficients
# -------------  Input:
# Δ         the floating point polynomial you want to approximate
# -------------  Optional input:
# tol       tolerance for the rationalization
# -------------  Output:
# the rational approximation
function rat(Δ; tolerance = 1e-8)
    vars = variables(Δ)
    E, C = exponents_coefficients(Δ, vars)
    mons = [prod(vars .^ E[:, i]) for i = 1:size(E, 2)]
    C = convert.(ComplexF64, C)
    C = real.(C)
    Crat = rationalize.(C; tol = tolerance)
    sum([Crat[i] * mons[i] for i = 1:length(mons)])
end

# Substitute inner products of momentum vectors by expressions in the Mandelstam invariants, internal and external masses
# in the second Symanzik polynomial for a diagram with 4 legs.
# -------------  Input:
# F         the second Symanzik polynomial of the 4-leg diagram
# p         the `Gram matrix' of the momentum vectors p[i,j] = pᵢ ⋅ pⱼ
# mm        a vector of internal mass parameters
# -------------  Optional input:
# equalM    boolean variable, set to true if you want to assume that all external particles have the same mass
# equalm    boolean variable, set to true if you want to assume that all internal particles have the same mass
# -------------  Output:
# Fsubs     the substituted Symanzik polynomial
# s         the mandelstam invariant (p₁ + p₂)²
# t         the mandelstam invariant (p₂ + p₃)²
# M         squared external masses
# m         squared internal masses
function substitute4legs(F, p, mm; equalM = false, equalm = false)
    @var s t
    if equalM
        @var M
        subslist = [
            p[1, 1] => M,
            p[2, 2] => M,
            p[3, 3] => M,
            p[4, 4] => M,
            p[1, 2] => (s / 2 - M),
            p[1, 3] => ((4 * M - s - t) / 2 - M),
            p[1, 4] => (t / 2 - M),
            p[2, 3] => (t / 2 - M),
            p[2, 4] => ((4 * M - s - t) / 2 - M),
            p[3, 4] => (s / 2 - M),
        ]
    else
        @var M[1:4]
        subslist = [
            p[1, 1] => M[1],
            p[2, 2] => M[2],
            p[3, 3] => M[3],
            p[4, 4] => M[4],
            p[1, 2] => (s - M[1] - M[2]) / 2,
            p[1, 3] => ((sum(M) - s - t) - M[1] - M[3]) / 2,
            p[1, 4] => (t - M[1] - M[4]) / 2,
            p[2, 3] => (t - M[2] - M[3]) / 2,
            p[2, 4] => ((sum(M) - s - t) - M[2] - M[4]) / 2,
            p[3, 4] => (s - M[3] - M[4]) / 2,
        ]
    end
    if equalm
        @var m
        subslist = [subslist; [ℓ => m for ℓ ∈ mm]]
    else
        @var m[1:length(mm)]
        subslist = [subslist; [mm[i] => m[i] for i = 1:length(mm)]]
    end
    Fsubs = subs(F, subslist...)
    return Fsubs, s, t, M, m
end

# Substitute inner products of momentum vectors by expressions in the Mandelstam invariants, internal and external masses
# in the second Symanzik polynomial for a diagram with 5 legs.
# -------------  Input:
# F         the second Symanzik polynomial of the 5-leg diagram
# p         the `Gram matrix' of the momentum vectors p[i,j] = pᵢ ⋅ pⱼ
# mm        a vector of internal mass parameters
# -------------  Optional input:
# equalM    boolean variable, set to true if you want to assume that all external particles have the same mass
# equalm    boolean variable, set to true if you want to assume that all internal particles have the same mass
# -------------  Output:
# Fsubs     the substituted Symanzik polynomial
# sij       the mandelstam invariant (pᵢ + pⱼ)²
# M         squared external masses
# m         squared internal masses
function substitute5legs(F, p, mm; equalM = false, equalm = false)
    @var s12, s23, s34, s45, s51
    if equalM
        @var M
        subslist = [
            p[1, 1] => M,
            p[2, 2] => M,
            p[3, 3] => M,
            p[4, 4] => M,
            p[5, 5] => M,
            p[1, 2] => (s12 - 2 * M) / 2,
            p[1, 3] => (-s12 - s23 + s45 + M) / 2,
            p[1, 4] => (s23 - s45 - s51 + M) / 2,
            p[1, 5] => (s51 - 2 * M) / 2,
            p[2, 3] => (s23 - 2 * M) / 2,
            p[2, 4] => (-s23 - s34 + s51 + M) / 2,
            p[2, 5] => (-s12 + s34 - s51 + M) / 2,
            p[3, 4] => (s34 - 2 * M) / 2,
            p[3, 5] => (s12 - s34 - s45 + M) / 2,
            p[4, 5] => (s45 - 2 * M) / 2,
        ]
    else
        @var M[1:5]
        subslist = [
            p[1, 1] => M[1],
            p[2, 2] => M[2],
            p[3, 3] => M[3],
            p[4, 4] => M[4],
            p[5, 5] => M[5],
            p[1, 2] => (s12 - M[1] - M[2]) / 2,
            p[1, 3] => (-s12 - s23 + s45 + M[2]) / 2,
            p[1, 4] => (s23 - s45 - s51 + M[5]) / 2,
            p[1, 5] => (s51 - M[1] - M[5]) / 2,
            p[2, 3] => (s23 - M[2] - M[3]) / 2,
            p[2, 4] => (-s23 - s34 + s51 + M[3]) / 2,
            p[2, 5] => (-s12 + s34 - s51 + M[1]) / 2,
            p[3, 4] => (s34 - M[3] - M[4]) / 2,
            p[3, 5] => (s12 - s34 - s45 + M[4]) / 2,
            p[4, 5] => (s45 - M[4] - M[5]) / 2,
        ]
    end
    if equalm
        @var m
        subslist = [subslist; [ℓ => m for ℓ ∈ mm]]
    else
        @var m[1:length(mm)]
        subslist = [subslist; [mm[i] => m[i] for i = 1:length(mm)]]
    end
    Fsubs = subs(F, subslist...)
    return Fsubs, s12, s23, s34, s45, s51, M, m
end

# Sample the projection of a variety given by polynomial equations on ℂᵐ × ℂᵉ onto ℂᵉ.
# This is the high precision version of `sampleProjection`. The coordinates of the samples are BigFloats.
# -------------  Input:
# F         an array of polynomials giving the polynomial system in variables x and p.
# x and p   two disjoint groups of coordinates, x on ℂᵐ, p on ℂᵉ.
# -------------  Optional input:
# The system F(x,p) is augmented with linear equations A*p + B = 0. If a solution (x,p) is known for some A, B,
# it can be passed to this function via seedA, seedB, seedsol. Other solutions will then be computed via monodromy on (A,B).
# codimen           expected codimension of the projection
# distinct_tol      tolerance for deciding whether two projected points are the same
# npoints           a lower bound on the number of samples that need to be computed
# printprogress     boolean variable, set to true if you want information to be printed throughout the execution
# -------------  Output:
# samples       sample points on the projection
# monres        result of the initial solving step using HomotopyContinuation.jl
# AA, BB        parameters used for the linear slice AA*p + BB = 0 in the initial solving step
# H             augmented system of equations (F(x,p);AA*p+BB) = 0
function sampleProjection_HP(
    F,
    x,
    p;
    seedA = nothing,
    seedB = nothing,
    seedsol = nothing,
    codimen = 1,
    distinct_tol = 1e-6,
    npoints = 100,
    printprogress = false,
)
    e = length(p)
    s = e - codimen
    @var A[1:s, 1:e] b[1:s]
    L = A * p + b
    H = System([F; L], variables = [p; x], parameters = [A[:]; b])
    if isnothing(seedsol)
        targetpars = randn(ComplexF64, length(parameters(H)))
        monres = solve(H, target_parameters = targetpars)
    else
        monres = monodromy_solve(H, seedsol, [seedA[:]; seedB])
        targetpars = parameters(monres)
    end
    refinedsols = [
        refine(
            System(
                subs(H.expressions, parameters(H) => targetpars),
                variables = variables(H),
            ),
            sol,
        ) for sol ∈ solutions(monres)
    ]
    refinedsols = [
        convert.(BigFloat, real(sol)) + im * convert.(BigFloat, imag(sol))
        for sol ∈ refinedsols
    ]
    psols = [sol[1:e] for sol ∈ refinedsols]
    samples, inds = distinct_inds(psols, tol = distinct_tol)
    startsamples = samples
    if length(samples) == 0
        println("Found no initial samples, are you sure about the codimension?")
    else
        while length(samples) < npoints
            newparams = randn(ComplexF64, length(targetpars))
            R = solve(
                H,
                solutions(monres);
                start_parameters = targetpars,
                target_parameters = newparams,
            )
            refinedsols = [
                refine(
                    System(
                        subs(H.expressions, parameters(H) => newparams),
                        variables = variables(H),
                    ),
                    sol,
                ) for sol ∈ solutions(R)
            ]
            refinedsols = [
                convert.(BigFloat, real(sol)) + im * convert.(BigFloat, imag(sol))
                for sol ∈ refinedsols
            ]
            psols = [sol[1:e] for sol ∈ refinedsols]
            newsamples, inds = distinct_inds(psols, tol = distinct_tol)
            samples = [samples; newsamples]
            if printprogress
                println("Collected $(length(samples)) samples")
            end
        end
    end
    return samples, monres, targetpars[1:length(A[:])], targetpars[length(A[:])+1:end], H
end

# decompose a set of solutions to F(x;p) = 0 into the groups lying on different irreducible components of the incidence variety.
# -------------  Input:
# F             an array of polynomials giving the polynomial system in variables x and p.
# seedsols      seed solutions: values of x for which F(x^*;p^*) = 0
# seedparams    seed parameters: values of p^* for which F(x^*;p^*) = 0
# -------------  Optional input:
# distinct_tol      tolerance for deciding whether two points are the same
# -------------  Output:
# representatives       a subset of solutions from seedsols for which each lies on a different irreducible component.
function decompose(F, seedsols, seedparams; distinct_tol = 1e-6)
    monres = monodromy_solve(F, seedsols[1], seedparams)
    representatives = [seedsols[1]]
    covered = solutions(monres)
    i = 2
    while length(covered) < length(seedsols) && i <= length(seedsols)
        if minimum(norm.(covered .- fill(seedsols[i], length(covered)))) > distinct_tol
            monres = monodromy_solve(F, seedsols[i], seedparams)
            representatives = [representatives; [seedsols[i]]]
            covered = [covered; solutions(monres)]
        end
        i += 1
    end
    return representatives
end

# Interpolate a set of points by a polynomial of a fixed degree by first learning its support and then using
# it for sparse interpolation.
# -------------  Input:
# samples   a list of sample points that you want to interpolate
# deg       degree of the interpolant
# parlist   a list of variables
# -------------  Optional input:
# homogeneous   boolean variable, set to true if you want all monomials of the fixed degree,
#               set to false if you want all monomials up to the fixed degree
# -------------  Output:
# pol           an interpolant
# sparsemons    monomials appearing in pol
# coeffs        coefficients of pol
function getInterpolantSparse(samples, parlist, deg; homogeneous = true, tol = 1e-9)
    mons = getMonomials(parlist, deg; homogeneous = homogeneous)
    V, vdmvecs = getVdm(samples, mons, parlist)
    VtV = V' * V
    EVobj = eigs(VtV; nev = 1, which = :SM)
    abseigenvec = abs.(EVobj[2][:, 1])
    nonzinds = findall(ℓ -> ℓ > tol, abseigenvec[:])
    sparsemons = mons[nonzinds]
    coeffs = EVobj[2][nonzinds, 1][:]
    maxind = findfirst(ℓ -> abs(ℓ) == maximum(abs.(coeffs)), coeffs)
    coeffs = coeffs / coeffs[maxind]
    pol = dot(sparsemons, coeffs)
    return pol, sparsemons, coeffs
end

# Compute the Landau discriminant of a given diagram.
# -------------  Input:
# edges     a list of edges of the diagram
# nodes     a list of nodes of the diagram to which external legs are attached
# -------------  Optional input:
# equalM    boolean variable, set to true if you want to assume that all external particles have the same mass
# equalm    boolean variable, set to true if you want to assume that all internal particles have the same mass
# HP        boolean variable, set to true if you want to use higher precision for the sampling
# mindegHP  a component is automatically sampled with higher precision if its degree is greater than or equal to mindegHP
# -------------  Output:
# ratpols   a list of polynomials with rational coefficients defining the irreducible components of the Landau discriminant
# parlist   a list of parameters in which the discriminant is expressed
function getLandauDiscriminant(
    edges,
    nodes;
    equalM = true,
    equalm = true,
    HP = false,
    mindegHP = 15,
)
    ALE, y, α, p, mm = affineLandauEquations(edges, nodes)
    if length(nodes) == 4
        ALE, s, t, M, m = substitute4legs(ALE, p, mm; equalM = equalM, equalm = equalm)
        parlist = [M; m; s; t]
    else
        ALE, s12, s23, s34, s45, s51, M, m =
            substitute5legs(ALE, p, mm; equalM = true, equalm = true)
        parlist = [s12, s23, s34, s45, s51, M, m]
    end
    println("-------------------------------------------")
    println("Slicing the incidence variety...")
    samples, monres, AA, bb, H, singsamples = sampleProjection(
        ALE,
        [α[1:end-1]; y],
        parlist;
        npoints = 1,
        printprogress = true,
        findSingular = true,
    )

    if length(singsamples) == 0
        println("Irreducible decomposition using monodromy...")
        reps = decompose(H, solutions(monres), [AA[:]; bb])
        println("Found $(length(reps)) components. ")
        println("Computing degrees of the components...")
        degs = [
            degreeProjection(
                ALE,
                [α[1:end-1]; y],
                parlist,
                seedsol = rep,
                seedA = AA,
                seedB = bb,
            ) for rep ∈ reps
        ]
        totaldeg = degreeProjection(ALE, [α[1:end-1]; y], parlist)
        println("Components have degrees $degs, total estimated degree is $totaldeg.")
        pols = []
        ratpols = []
        for i = 1:length(reps)
            println("-------------------------------------------")
            println("Sampling component $i ...")
            M = round(1.2 * binomial(degs[i] + length(parlist) - 1, length(parlist) - 1))
            if !HP && degs[i] < mindegHP
                samples, monres, AA, bb, H, singsamples = sampleProjection(
                    ALE,
                    [α[1:end-1]; y],
                    parlist;
                    npoints = M,
                    printprogress = false,
                    seedsol = reps[i],
                    seedA = AA,
                    seedB = bb,
                    findSingular = true,
                )
            else
                samples, monres, AA, bb, H = sampleProjection_HP(
                    ALE,
                    [α[1:end-1]; y],
                    parlist;
                    npoints = M,
                    printprogress = false,
                    seedsol = reps[i],
                    seedA = AA,
                    seedB = bb,
                )
            end
            println("Interpolating $(length(samples)) samples ...")
            if !HP && degs[i] < mindegHP
                Δ, sparsemons, coeff = getInterpolantSparse(samples, parlist, degs[i])
            else
                Δ, N, gap = interpolate_deg(samples, degs[i], parlist)
            end
            ratpol = rat(Δ)
            ratpols = [ratpols; ratpol]
            println(ratpol)
        end
    else
        deg = degreeProjection(ALE, [α[1:end-1]; y], parlist; findSingular = true)
        M = round(1.2 * binomial(deg + length(parlist) - 1, length(parlist) - 1))
        samples, monres, AA, bb, H, singsamples = sampleProjection(
            ALE,
            [α[1:end-1]; y],
            parlist;
            npoints = M,
            printprogress = false,
            findSingular = true,
        )
        samples = [samples; singsamples]
        println("Interpolating $(length(samples)) samples ...")
        Δ, sparsemons, coeff = getInterpolantSparse(samples, parlist, deg)
        ratpol = rat(Δ)
        ratpols = [ratpol]
        println(ratpol)
    end
    return ratpols, parlist
end

# Compute the distinct elements in a vector v, as well as an index in v for each of these elements.
# -------------  Input:
# v     a vector of scalars or arrays of (possibly floating point) numbers
# -------------  Optional input:
# tol   tolerance for deciding when two elements of v are the same
# -------------  Output:
# newv  a vector with the distinct elements of v
# inds  a vector of indices in v such that newv = v[inds]
function distinct_inds(v; tol = 1e-6)
    newv = [v[1]]
    inds = [1]
    counter = 2
    for a ∈ v[2:end]
        if minimum([norm(a - b) / norm(a) for b ∈ newv]) > tol
            push!(newv, a)
            push!(inds, counter)
        else
            #println(a)
        end
        counter += 1
    end
    return newv, inds
end

# Computes the number of master integrals of a Feynman diagram.
# -------------  Input:
# edges     a list of edges of the diagram
# nodes     a list of nodes of the diagram to which external legs are attached
# -------------  Output:
# the number of master integrals
function χ(edges, nodes)
    E = length(edges)
    F, U, α, p, mm = getF(edges, nodes)
    F, s, t, M, m = substitute4legs(F, p, mm)
    @var u[1:E+2]
    W = u[1] * log(U) + u[2] * log(F) + dot(u[3:E+2], log.(α))
    dW = System(differentiate(subs(W, α[E] => 1), α[1:E-1]), parameters = [s; t; M; m; u])
    Crit = monodromy_solve(dW)
    crt = certify(dW, Crit)
    return ndistinct_certified(crt)
end


end # module
