# --------------------------
# 1D: barycentric weights for Gauss–Legendre nodes
# λ_i = (-1)^(n-i) * sqrt((1-x_i^2) w_i / 2)  (up to a global scaling)
# x must be sorted ascending and wG in the same order.
# --------------------------
function gl_barycentric_weights(x::AbstractVector, wG::AbstractVector)
    n = length(x)
    λ = similar(x, Float64)
    for i in 1:n
        s = (-1)^(n - i)
        λ[i] = s * sqrt((1 - x[i]^2) * wG[i] / 2)
    end
    # optional normalization (does not change the interpolant)
    λ ./= maximum(abs.(λ))
    return λ
end

function barycentric_row(x::AbstractVector, λ::AbstractVector, xq::Real)
    n = length(x)
    r = zeros(Float64, n)

    hit = findfirst(i -> isapprox(xq, x[i]), 1:n)   # same as xq ≈ x[i]
    if hit !== nothing
        r[hit] = 1.0
        return r
    end

    denom = 0.0
    for i in 1:n
        t = λ[i] / (xq - x[i])
        r[i] = t
        denom += t
    end
    r ./= denom
    return r
end

function interp_matrix_1d_gl(x::AbstractVector, wG::AbstractVector, xq::AbstractVector)
    λ = gl_barycentric_weights(x, wG)
    m, n = length(xq), length(x)
    E = zeros(Float64, m, n)
    for j in 1:m
        E[j, :] .= barycentric_row(x, λ, xq[j])
    end
    return E
end

function interp_matrix_2d_gl_scattered(x::AbstractVector, wx::AbstractVector,
                                       y::AbstractVector, wy::AbstractVector,
                                       xt::AbstractVector, yt::AbstractVector)
    @assert length(xt) == length(yt)
    nx, ny = length(x), length(y)
    Nt = length(xt)

    λx = gl_barycentric_weights(x, wx)
    λy = gl_barycentric_weights(y, wy)

    A = zeros(Float64, Nt, nx * ny)
    for k in 1:Nt
        rx = barycentric_row(x, λx, xt[k])
        ry = barycentric_row(y, λy, yt[k])
        A[k, :] .= kron(ry, rx)   # vec(F) is column-major: i + (j-1)*nx
    end
    return A
end

function interp_matrix_2d_gl_tensor(x::AbstractVector, wx::AbstractVector,
                                    y::AbstractVector, wy::AbstractVector,
                                    Xout::AbstractVector, Yout::AbstractVector)
    Ex = interp_matrix_1d_gl(x, wx, Xout)
    Ey = interp_matrix_1d_gl(y, wy, Yout)
    return kron(Ey, Ex)
end