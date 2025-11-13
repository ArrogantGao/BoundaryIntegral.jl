solve_lu(A::AbstractMatrix, b::AbstractVector) = lu(A) \ b

function solve_gmres(A, b::Vector{T}, atol::T = 1e-6, itmax::Int = 1000) where T
    x = IterativeSolvers.gmres(A, b, abstol = atol, maxiter = itmax)
    return x
end

function solve_gmres(A, b::CuVector{T}, atol::T = 1e-4, itmax::Int = 1000) where T
    x, status = Krylov.gmres(A, b, atol = atol, itmax = itmax)
    return x
end