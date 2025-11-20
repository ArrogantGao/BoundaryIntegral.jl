solve_lu(A::AbstractMatrix, b::AbstractVector) = lu(A) \ b

function solve_gmres(A, b, atol::T = 1e-6, rtol::T = 1e-6, itmax::Int = 1000) where T
    x, _ = Krylov.gmres(A, b, atol = atol, rtol = rtol, itmax = itmax)
    return x
end