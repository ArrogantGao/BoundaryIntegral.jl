function laplace3d_singlelayer(src::NTuple{3, T}, trg::NTuple{3, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    inv_r = one(T) / r

    return inv_r / 4π
end

function laplace3d_doublelayer(src::NTuple{3, T}, trg::NTuple{3, T}, norm::NTuple{3, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    inv_r = one(T) / r

    return dot(norm, inv_r^3 .* (trg .- src))  / 4π
end

function laplace3d_DT(dielectric_interfaces::DielectricInterfaces{T, 3}) where T
    n_points = num_points(dielectric_interfaces)
    weights = all_weights(dielectric_interfaces)

    DT = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(dielectric_interfaces))
        for (j, pointj) in enumerate(eachpoint(dielectric_interfaces))
            i == j && continue
            DT[i, j] = - laplace3d_doublelayer(pointj.point, pointi.point, pointi.normal)
        end
    end
    return DT * diagm(weights)
end

function laplace3d_D(dielectric_interfaces::DielectricInterfaces{T, 3}) where T
    n_points = num_points(dielectric_interfaces)
    weights = all_weights(dielectric_interfaces)

    D = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(dielectric_interfaces))
        for (j, pointj) in enumerate(eachpoint(dielectric_interfaces))
            i == j && continue
            D[j, i] = laplace3d_doublelayer(pointj.point, pointi.point, pointi.normal)
        end
    end
    return D * diagm(weights)
end

function _laplace3d_DT_fmm3d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64) 
    n = length(charges)
    @assert size(sources) == (3, n)
    @assert size(norms) == (3, n)
    @assert size(weights) == (n,)
    eps = thresh
    vals = lfmm3d(eps, sources, charges = weights .* charges, pg = 2)
    grad = vals.grad
    gradn = zeros(Float64, n)

    for i in 1:n
        gradn[i] = dot(norms[:, i], grad[:, i])
    end

    return - gradn ./ 4π
end

function laplace3d_DT_fmm3d(dielectric_interfaces::DielectricInterfaces{Float64, 3}, thresh::Float64)
    n_points = num_points(dielectric_interfaces)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 3, n_points)
    for (i, point) in enumerate(eachpoint(dielectric_interfaces))
        weights[i] = point.weight
        sources[1, i] = point.point[1]
        sources[2, i] = point.point[2]
        sources[3, i] = point.point[3]
        norms[1, i] = point.normal[1]
        norms[2, i] = point.normal[2]
        norms[3, i] = point.normal[3]
    end

    f = charges -> _laplace3d_DT_fmm3d(charges, sources, weights, norms, thresh)
    return LinearMap{Float64}(f, n_points, n_points)
end

function _laplace3d_D_fmm3d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64) 
    n = length(weights)
    @assert length(charges) == n
    @assert size(sources) == (3, n)
    @assert size(norms) == (3, n)

    dipvecs = zeros(Float64, 3, n)
    for i in 1:n
        dipvecs[1, i] = norms[1, i] * (charges[i] * weights[i])
        dipvecs[2, i] = norms[2, i] * (charges[i] * weights[i])
        dipvecs[3, i] = norms[3, i] * (charges[i] * weights[i])
    end

    vals = lfmm3d(thresh, sources, dipvecs = dipvecs, pg = 1)

    return vals.pot ./ 4π
end

function laplace3d_D_fmm3d(dielectric_interfaces::DielectricInterfaces{Float64, 3}, thresh::Float64)
    n_points = num_points(dielectric_interfaces)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 3, n_points)
    for (i, point) in enumerate(eachpoint(dielectric_interfaces))
        weights[i] = point.weight
        sources[1, i] = point.point[1]
        sources[2, i] = point.point[2]
        sources[3, i] = point.point[3]
        norms[1, i] = point.normal[1]
        norms[2, i] = point.normal[2]
        norms[3, i] = point.normal[3]
    end

    f = charges -> _laplace3d_D_fmm3d(charges, sources, weights, norms, thresh)
    return LinearMap{Float64}(f, n_points, n_points)
end

function _laplace3d_D_trg_fmm3d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, targets::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64) 
    n = length(weights)
    m = size(targets, 2)
    @assert length(charges) == n
    @assert size(sources) == (3, n)
    @assert size(norms) == (3, n)

    dipvecs = zeros(Float64, 3, n)
    for i in 1:n
        dipvecs[1, i] = norms[1, i] * (charges[i] * weights[i])
        dipvecs[2, i] = norms[2, i] * (charges[i] * weights[i])
        dipvecs[3, i] = norms[3, i] * (charges[i] * weights[i])
    end

    vals = lfmm3d(thresh, sources, dipvecs = dipvecs, targets = targets, pgt = 1)

    return vals.pottarg ./ 4π
end

function laplace3d_D_trg_fmm3d(dielectric_interfaces::DielectricInterfaces{Float64, 3}, targets::Matrix{Float64}, thresh::Float64)
    n_points = num_points(dielectric_interfaces)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 3, n_points)
    for (i, point) in enumerate(eachpoint(dielectric_interfaces))
        weights[i] = point.weight
        sources[1, i] = point.point[1]
        sources[2, i] = point.point[2]
        sources[3, i] = point.point[3]
        norms[1, i] = point.normal[1]
        norms[2, i] = point.normal[2]
        norms[3, i] = point.normal[3]
    end

    f = charges -> _laplace3d_D_trg_fmm3d(charges, sources, targets, weights, norms, thresh)
    return LinearMap{Float64}(f, size(targets, 2), n_points)
end

function _laplace3d_pottarg_fmm3d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, targets::Matrix{Float64}, thresh::Float64)
    n = length(charges)
    m = size(targets, 2)
    @assert size(sources) == (3, n)
    @assert size(targets) == (3, m)
    @assert size(weights) == (n,)
    vals = lfmm3d(thresh, sources, charges = weights .* charges, targets = targets, pgt = 1)
    return vals.pottarg ./ 4π
end

function laplace3d_pottarg_fmm3d(dielectric_interfaces::DielectricInterfaces{Float64, 3}, targets::Matrix{Float64}, thresh::Float64)
    n_points = num_points(dielectric_interfaces)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    for (i, point) in enumerate(eachpoint(dielectric_interfaces))
        weights[i] = point.weight
        sources[1, i] = point.point[1]
        sources[2, i] = point.point[2]
        sources[3, i] = point.point[3]
    end

    f = charges -> _laplace3d_pottarg_fmm3d(charges, sources, weights, targets, thresh)
    return LinearMap{Float64}(f, size(targets, 2), n_points)
end