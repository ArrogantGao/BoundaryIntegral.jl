"""
    legendreP_matrix(x, M)

Return an M×N matrix P where P[k,j] = P_{k-1}(x[j]),
the (unnormalized) Legendre polynomial of degree k-1 evaluated at x[j].

Uses the standard recurrence:
P₀(x)=1, P₁(x)=x,
(n+1) P_{n+1}(x) = (2n+1)x P_n(x) - n P_{n-1}(x).
"""
function legendreP_matrix(x::AbstractVector{T}, M::Integer) where {T}
    N = length(x)
    P = Matrix{T}(undef, M, N)

    # P₀(x) = 1
    @inbounds for j in 1:N
        P[1, j] = one(T)
    end

    # P₁(x) = x
    if M > 1
        @inbounds for j in 1:N
            P[2, j] = x[j]
        end
    end

    # Recurrence for n = 1,2,...,M-2 giving P_{n+1}
    for n in 1:(M-2)
        np1 = n + 1
        @inbounds for j in 1:N
            P[np1 + 1, j] = ((2n + 1) * x[j] * P[np1, j] - n * P[n, j]) / (n + 1)
        end
    end

    return P
end

"""
    legendre_transform_matrix(x, w; M=length(x)) -> T

Build the M×N matrix T such that, for values f[i] = f(x[i]),

    c = T * f

gives Legendre coefficients c[k+1] ≈ ∫_{-1}^1 f(x) P_k(x) dx * (2k+1)/2.
"""
function legendre_transform_matrix(x::AbstractVector,
                                   w::AbstractVector;
                                   M::Integer = length(x))
    N = length(x)
    @assert length(w) == N "weights length must match nodes length"
    @assert M <= N "for exactness on polynomials, require M ≤ N"

    T = legendreP_matrix(x, M)  # M×N
    @inbounds for k in 0:(M-1)
        row = k + 1
        s = (2k + 1) / 2
        for j in 1:N
            T[row, j] *= s * w[j]
        end
    end
    return T
end

"""
    values_to_legendre_coeffs(x, w, f; M=length(x))

Matrix form: c = T * f.
"""
function values_to_legendre_coeffs(x::AbstractVector,
                                   w::AbstractVector,
                                   f::AbstractVector;
                                   M::Integer = length(x))
    T = legendre_transform_matrix(x, w; M=M)
    @assert length(f) == size(T, 2)
    return T * f
end

"""
    eval_legendre_series(x, c)

Evaluate the Legendre series

    S(x) = ∑_{k=0}^{M-1} c[k+1] * P_k(x)

at a scalar `x` ∈ [-1,1].

Uses the standard three-term recurrence for P_k.
"""
function eval_legendre_series(x::Real, c::AbstractVector)
    M = length(c)
    T = promote_type(eltype(c), typeof(x))
    xx = T(x)

    if M == 0
        return zero(T)
    end

    # P₀(x) = 1
    Pnm2 = one(T)
    s = c[1] * Pnm2

    if M == 1
        return s
    end

    # P₁(x) = x
    Pnm1 = xx
    s += c[2] * Pnm1

    # n runs 1,...,M-2 to generate P_{n+1}
    @inbounds for n in 1:(M-2)
        Pn = ((2n + 1) * xx * Pnm1 - n * Pnm2) / (n + 1)
        s += c[n + 2] * Pn
        Pnm2, Pnm1 = Pnm1, Pn
    end

    return s
end

"""
    eval_legendre_series(x::AbstractVector, c)

Vectorized version: evaluate the Legendre series at each x[j].
Returns a vector of the same length as `x`.
"""
function eval_legendre_series(x::AbstractVector, c::AbstractVector)
    T = promote_type(eltype(x), eltype(c))
    y = similar(x, T)
    @inbounds for j in eachindex(x)
        y[j] = eval_legendre_series(x[j], c)
    end
    return y
end

"""
    legendre_eval_matrix(x_target, M) -> E

Build the evaluation matrix E of size length(x_target) × M such that

    f_vals ≈ E * c

where c[k+1] are Legendre coefficients and
E[j,k] = P_{k-1}(x_target[j]).
"""
function legendre_eval_matrix(x_target::AbstractVector, M::Integer)
    # legendreP_matrix gives M×N; transpose to N×M.
    P = legendreP_matrix(x_target, M)  # M×N_target
    return permutedims(P)              # N_target×M
end


# -----------------------------------------------------------
# 2D transform matrix (tensor product)
# -----------------------------------------------------------

"""
    legendre_transform_matrix_2d(x, wx, y, wy; Mx=length(x), My=length(y)) -> T2D

Build the matrix T2D of size (Mx*My) × (Nx*Ny) such that for a vector `f`
of nodal values on the tensor-product grid (stored as `vec(F)` with F[i,j]),
the 2D Legendre coefficients `c` (for basis P_m(x)P_n(y)) satisfy

    c = T2D * f

where `c` is ordered as `vec(C)` with C[m+1,n+1].
"""
function legendre_transform_matrix_2d(x::AbstractVector, wx::AbstractVector,
                                      y::AbstractVector, wy::AbstractVector;
                                      Mx::Integer = length(x),
                                      My::Integer = length(y))

    Tx = legendre_transform_matrix(x, wx; M = Mx)  # Mx×Nx
    Ty = legendre_transform_matrix(y, wy; M = My)  # My×Ny

    # For F (Nx×Ny), vec(C) = (Ty ⊗ Tx) * vec(F)
    return kron(Ty, Tx)
end

# Convenience: apply 2D transform to a vector of nodal values
function values_to_legendre_coeffs_2d(x::AbstractVector, wx::AbstractVector,
                                      y::AbstractVector, wy::AbstractVector,
                                      f::AbstractVector;
                                      Mx::Integer = length(x),
                                      My::Integer = length(y))
    Nx = length(x)
    Ny = length(y)
    @assert length(wx) == Nx
    @assert length(wy) == Ny
    @assert length(f) == Nx * Ny

    T2D = legendre_transform_matrix_2d(x, wx, y, wy; Mx=Mx, My=My)
    return T2D * f          # length Mx*My
end

# Optional: matrix form if you prefer working with F as Nx×Ny
function values_to_legendre_coeffs_2d(x::AbstractVector, wx::AbstractVector,
                                      y::AbstractVector, wy::AbstractVector,
                                      F::AbstractMatrix;
                                      Mx::Integer = length(x),
                                      My::Integer = length(y))
    Nx, Ny = size(F)
    @assert length(x)  == Nx
    @assert length(y)  == Ny
    @assert length(wx) == Nx
    @assert length(wy) == Ny

    Tx = legendre_transform_matrix(x, wx; M = Mx)  # Mx×Nx
    Ty = legendre_transform_matrix(y, wy; M = My)  # My×Ny

    C = Tx * F * transpose(Ty)   # Mx×My
    return C                    # or vec(C) if you want a vector
end



# -----------------------------------------------------------
# 2D evaluation matrix
# -----------------------------------------------------------

"""
    legendre_eval_matrix_2d(x_target, y_target; Mx, My) -> E2D

Build the matrix E2D of size (Nx*Ny) × (Mx*My) such that

    fvec = E2D * cvec

where:
- `x_target` (length Nx), `y_target` (length Ny) define the tensor grid,
- `cvec = vec(C)` with C[m+1,n+1] the 2D Legendre coefficients
  in basis P_m(x)P_n(y),
- `fvec = vec(F)` with F[i,j] = f(x_target[i], y_target[j]).
"""
function legendre_eval_matrix_2d(x_target::AbstractVector,
                                 y_target::AbstractVector;
                                 Mx::Integer,
                                 My::Integer)
    Ex = legendre_eval_matrix(x_target, Mx)  # Nx×Mx
    Ey = legendre_eval_matrix(y_target, My)  # Ny×My

    # vec(F) = (Ey ⊗ Ex) * vec(C)
    return kron(Ey, Ex)
end