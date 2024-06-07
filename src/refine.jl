using HomotopyContinuation, Arblib, LinearAlgebra

# This is code for refining a solution to a system of polynomial equations. That is,
# to compute the coordinates in higher precision.
# The code was kindly shared with the authors by Sascha Timme.

#=== examples
@var x y
f = System([y^2 - 2x - 3, x - 1])

# real input
x0 = [1,√5+1e-12]
x1 = refine(f, x0, atol = 1e-12)
# as bigfloat
big_x1 = BigFloat.(x1)

# complex input
z0 = [1+1e-7im,√5+1e-12]
z1 = refine(f, z0, prec = 512)
=#

function refine(f::System, x0::AbstractVector; prec = 256, atol = 0.0, iters = 10)
    F = InterpretedSystem(f)
    n, = size(F)
    y = ArbRefMatrix(x0; prec = prec)
    Δy = ArbRefMatrix(n, 1; prec = prec)
    u = ArbRefMatrix(n,1; prec = prec)
    U = ArbRefMatrix(n, n; prec = prec)
    mag = Arblib.Mag()
    for k in 1:iters
        evaluate_and_jacobian!(u, U, F, y)
        Arblib.get_mid!(u, u)
        Arblib.get_mid!(U, U)
        ldiv!(Δy, U, u)
        Arblib.get_mid!(Δy, Δy)
        Arblib.sub!(y, y, Δy)
        Arblib.get_mid!(y, y)
        Arblib.bound_inf_norm!(mag, Δy)
        if mag < atol
            break
        end
    end
    ArbVector(y[:,1])
end

function refine(f::System, x0::AbstractVector{<:Complex}; prec = 256, atol = 0.0, iters = 10)
    F = InterpretedSystem(f)
    n, = size(F)
    y = AcbRefMatrix(x0; prec = prec)
    Δy = AcbRefMatrix(n, 1; prec = prec)
    u = AcbRefMatrix(n,1; prec = prec)
    U = AcbRefMatrix(n, n; prec = prec)
    mag = Arblib.Mag()
    for k in 1:iters
        evaluate_and_jacobian!(u, U, F, y)
        Arblib.get_mid!(u, u)
        Arblib.get_mid!(U, U)
        ldiv!(Δy, U, u)
        Arblib.get_mid!(Δy, Δy)
        Arblib.sub!(y, y, Δy)
        Arblib.get_mid!(y, y)
        Arblib.bound_inf_norm!(mag, Δy)
        if mag < atol
            break
        end
    end
    AcbVector(y[:,1])
end
