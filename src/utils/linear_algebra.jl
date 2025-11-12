solve_lu(A::AbstractMatrix, b::AbstractVector) = lu(A) \ b

function solve_gmres(A, b::AbstractVector, reltol::Float64 = 1e-10, maxiter::Int = 1000)
    x = gmres(A, b, reltol = reltol, maxiter = maxiter)
    return x
end