# defination of the laplace kernel in 2d, potential and gradient
function laplace2d_pot(src::NTuple{2, T}, trg::NTuple{2, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    return log(r) / 2π
end

function laplace2d_grad(src::NTuple{2, T}, trg::NTuple{2, T}, norm::NTuple{2, T}) where T
    r2 = sum((src .- trg).^2)
    inv_r2 = one(T) / r2

    return dot(norm, inv_r2 .* (trg .- src)) / 2π
end

# filling the matrix directly, including S, D and DT
function laplace2d_S(interface::DielectricInterface{P, T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    weights = all_weights(interface)
    S = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(interface))
        for (j, pointj) in enumerate(eachpoint(interface))
            i == j && continue
            S[i, j] = laplace2d_pot(pointj.panel_point.point, pointi.panel_point.point)
        end
    end
    return S * diagm(weights)
end

function laplace2d_D(interface::DielectricInterface{P, T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    weights = all_weights(interface)

    D = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(interface))
        for (j, pointj) in enumerate(eachpoint(interface))
            i == j && continue
            D[j, i] = laplace2d_grad(pointj.panel_point.point, pointi.panel_point.point, pointi.panel_point.normal)
        end
    end

    return D * diagm(weights)
end

function laplace2d_DT(interface::DielectricInterface{P, T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    weights = all_weights(interface)

    DT = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(interface))
        for (j, pointj) in enumerate(eachpoint(interface))
            i == j && continue
            DT[i, j] = laplace2d_grad(pointj.panel_point.point, pointi.panel_point.point, pointi.panel_point.normal)
        end
    end

    return DT * diagm(weights)
end

function laplace2d_pottrg(interface::DielectricInterface{P, T}, targets::Matrix{T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    n_targets = size(targets, 2)
    weights = all_weights(interface)
    pot = zeros(T, n_points, n_targets)
    for (i, pointi) in enumerate(eachpoint(interface))
        for j in 1:n_targets
            target = (targets[1, j], targets[2, j])
            pot[i, j] = laplace2d_pot(pointi.panel_point.point, target)
        end
    end
    return pot' * diagm(weights)
end

# fmm2d based fast evaluation
# matrix free via fmm2d, return a linear map
# potential: log(r), gradient: 1/r
function _laplace2d_S_fmm2d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, thresh::Float64)
    n = length(charges)
    @assert size(sources) == (2, n)
    @assert size(weights) == (n,)
    vals = rfmm2d(eps = thresh, sources = sources, charges = weights .* charges, pg = 1)
    return vals.pot ./ 2π
end

function laplace2d_S_fmm2d(interface::DielectricInterface{P, Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 2, n_points)
    weights = zeros(Float64, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
    end
    f = charges -> _laplace2d_S_fmm2d(charges, sources, weights, thresh)
    return LinearMap{Float64}(f, n_points, n_points)
end

function _laplace2d_DT_fmm2d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64)
    n = length(charges)
    @assert size(sources) == (2, n)
    @assert size(norms) == (2, n)
    @assert size(weights) == (n,)
    vals = rfmm2d(eps = thresh, sources = sources, charges = weights .* charges, pg = 2)
    grad = vals.grad
    gradn = zeros(Float64, n)

    for i in 1:n
        gradn[i] = dot(norms[:, i], grad[:, i])
    end

    return gradn ./ 2π
end


function laplace2d_DT_fmm2d(interface::DielectricInterface{P, Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 2, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 2, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
        norms[1, i] = point.panel_point.normal[1]
        norms[2, i] = point.panel_point.normal[2]
    end

    f = charges -> _laplace2d_DT_fmm2d(charges, sources, weights, norms, thresh)
    return LinearMap{Float64}(f, n_points, n_points)
end

function _laplace2d_D_fmm2d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64)
    n = length(weights)
    @assert size(sources) == (2, n)
    @assert size(norms) == (2, n)
    @assert size(weights) == (n,)

    dipvecs = zeros(Float64, 2, n)
    for i in 1:n
        dipvecs[1, i] = norms[1, i] * (charges[i] * weights[i])
        dipvecs[2, i] = norms[2, i] * (charges[i] * weights[i])
    end

    vals = rfmm2d(eps = thresh, sources = sources, dipvecs = dipvecs, pg = 1)
    return vals.pot ./ 2π
end

function laplace2d_D_fmm2d(interface::DielectricInterface{P, Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 2, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 2, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
        norms[1, i] = point.panel_point.normal[1]
        norms[2, i] = point.panel_point.normal[2]
    end

    f = charges -> _laplace2d_D_fmm2d(charges, sources, weights, norms, thresh)
    return LinearMap{Float64}(f, n_points, n_points)
end

function _laplace2d_pottrg_fmm2d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, targets::Matrix{Float64}, thresh::Float64)
    n = length(charges)
    m = size(targets, 2)

    @assert size(sources) == (2, n)
    @assert size(targets) == (2, m)
    @assert size(weights) == (n,)
    vals = rfmm2d(eps = thresh, sources = sources, charges = weights .* charges, targets = targets, pgt = 1)
    return vals.pottarg ./ 2π
end

function laplace2d_pottrg_fmm2d(interface::DielectricInterface{P, Float64}, targets::Matrix{Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 2, n_points)
    weights = zeros(Float64, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
    end

    f = charges -> _laplace2d_pottrg_fmm2d(charges, sources, weights, targets, thresh)
    return LinearMap{Float64}(f, size(targets, 2), n_points)
end

# I am not sure why we need to evaluate this
# function _laplace2d_gradtrg_fmm2d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, targets::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64)
#     n = length(charges)
#     m = size(targets, 2)
#     @assert size(sources) == (2, n)
#     @assert size(targets) == (2, m)
#     @assert size(weights) == (n,)
#     @assert size(norms) == (2, n)

#     dipvecs = zeros(Float64, 2, n)
#     for i in 1:n
#         dipvecs[1, i] = norms[1, i] * (charges[i] * weights[i])
#         dipvecs[2, i] = norms[2, i] * (charges[i] * weights[i])
#     end

#     vals = rfmm2d(eps = thresh, sources = sources, dipvecs = dipvecs, targets = targets, pgt = 1)
#     return - vals.pottarg ./ 2π
# end

# function laplace2d_gradtrg_fmm2d(interface::DielectricInterface{P, Float64}, targets::Matrix{Float64}, thresh::Float64) where {P <: AbstractPanel}
#     n_points = num_points(interface)
#     sources = zeros(Float64, 2, n_points)
#     weights = zeros(Float64, n_points)
#     norms = zeros(Float64, 2, n_points)
#     for (i, point) in enumerate(eachpoint(interface))
#         weights[i] = point.panel_point.weight
#         sources[1, i] = point.panel_point.point[1]
#         sources[2, i] = point.panel_point.point[2]
#         norms[1, i] = point.panel_point.normal[1]
#         norms[2, i] = point.panel_point.normal[2]
#     end

#     f = charges -> _laplace2d_gradtrg_fmm2d(charges, sources, targets, weights, norms, thresh)
#     return LinearMap{Float64}(f, size(targets, 2), n_points)
# end
