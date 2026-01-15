function laplace3d_pot(src::NTuple{3, T}, trg::NTuple{3, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    inv_r = one(T) / r

    return inv_r / 4π
end

function laplace3d_grad(src::NTuple{3, T}, trg::NTuple{3, T}, norm::NTuple{3, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    inv_r = one(T) / r

    return dot(norm, inv_r^3 .* (trg .- src))  / 4π
end

# direct evaluation of S, D and DT
function laplace3d_S(interface::DielectricInterface{P, T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    weights = all_weights(interface)
    S = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(interface))
        for (j, pointj) in enumerate(eachpoint(interface))
            i == j && continue
            S[i, j] = laplace3d_pot(pointj.panel_point.point, pointi.panel_point.point)
        end
    end
    return S * diagm(weights)
end

function laplace3d_DT(interface::DielectricInterface{P, T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    weights = all_weights(interface)

    DT = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(interface))
        for (j, pointj) in enumerate(eachpoint(interface))
            i == j && continue
            DT[i, j] = laplace3d_grad(pointj.panel_point.point, pointi.panel_point.point, pointi.panel_point.normal)
        end
    end
    return DT * diagm(weights)
end

function laplace3d_DT_panel(panel_src::FlatPanel{T, 3}, panel_trg::FlatPanel{T, 3}) where T

    np_src = num_points(panel_src)
    np_trg = num_points(panel_trg)

    DT = zeros(T, np_trg, np_src)
    for (i, pointi) in enumerate(eachpoint(panel_src))
        for (j, pointj) in enumerate(eachpoint(panel_trg))
            DT[j, i] = laplace3d_grad(pointi.point, pointj.point, pointj.normal)
        end
    end
    return DT * diagm(panel_src.weights)
end

# this function generate a block of the correction matrix
function laplace3d_DT_panel_upsampled(panel_src::FlatPanel{T, 3}, panel_trg::FlatPanel{T, 3}, n_up::Int) where T
    ns_up, ws_up = gausslegendre(n_up)
    ns_up = T.(ns_up)
    ws_up = T.(ws_up)

    ns0 = panel_src.gl_xs
    ws0 = panel_src.gl_ws

    M_up = interp_matrix_2d_gl_tensor(ns0, ws0, ns0, ws0, ns_up, ns_up)

    a, b, c, d = panel_src.corners
    cc = (a .+ b .+ c .+ d) ./ 4
    Lx = norm(b .- a)
    Ly = norm(d .- a)

    weights_up = Vector{T}(undef, n_up * n_up)
    points_up = Vector{NTuple{3, T}}(undef, n_up * n_up)
    idx = 1
    for j in 1:n_up
        for i in 1:n_up
            points_up[idx] = cc .+ (b .- a) .* (ns_up[i] / 2) .+ (d .- a) .* (ns_up[j] / 2)
            weights_up[idx] = ws_up[i] * ws_up[j] * Lx * Ly / 4
            idx += 1
        end
    end

    np_trg = num_points(panel_trg)
    D_up_trg = zeros(T, np_trg, n_up * n_up)
    for (ti, point_trg) in enumerate(eachpoint(panel_trg))
        for ui in 1:length(points_up)
            D_up_trg[ti, ui] = laplace3d_grad(points_up[ui], point_trg.point, point_trg.normal)
        end
    end

    DT_up = D_up_trg * diagm(weights_up) * M_up

    return DT_up
end

function laplace3d_D(interface::DielectricInterface{P, T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    weights = all_weights(interface)

    D = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(interface))
        for (j, pointj) in enumerate(eachpoint(interface))
            i == j && continue
            D[j, i] = laplace3d_grad(pointj.panel_point.point, pointi.panel_point.point, pointi.panel_point.normal)
        end
    end
    return D * diagm(weights)
end

function laplace3d_D_trg(interface::DielectricInterface{P, T}, targets::Matrix{T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    n_targets = size(targets, 2)
    weights = all_weights(interface)
    D = zeros(T, n_points, n_targets)
    for (i, pointi) in enumerate(eachpoint(interface))
        for j in 1:n_targets
            target = (targets[1, j], targets[2, j], targets[3, j])
            D[i, j] = laplace3d_grad(pointi.panel_point.point, target, pointi.panel_point.normal)
        end
    end
    return D' * diagm(weights)
end

function laplace3d_pottrg(interface::DielectricInterface{P, T}, targets::Matrix{T}) where {P <: AbstractPanel, T}
    n_points = num_points(interface)
    n_targets = size(targets, 2)
    weights = all_weights(interface)
    pot = zeros(T, n_points, n_targets)
    for (i, pointi) in enumerate(eachpoint(interface))
        for j in 1:n_targets
            target = (targets[1, j], targets[2, j], targets[3, j])
            pot[i, j] = laplace3d_pot(pointi.panel_point.point, target)
        end
    end
    return pot' * diagm(weights)
end

function _laplace3d_DT_fmm3d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64)
    n = length(charges)
    @assert size(sources) == (3, n)
    @assert size(norms) == (3, n)
    @assert size(weights) == (n,)
    vals = lfmm3d(thresh, sources, charges = weights .* charges, pg = 2)
    grad = vals.grad
    gradn = zeros(Float64, n)

    for i in 1:n
        gradn[i] = dot(norms[:, i], grad[:, i])
    end

    return - gradn ./ 4π
end

function laplace3d_DT_fmm3d(interface::DielectricInterface{P, Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 3, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
        sources[3, i] = point.panel_point.point[3]
        norms[1, i] = point.panel_point.normal[1]
        norms[2, i] = point.panel_point.normal[2]
        norms[3, i] = point.panel_point.normal[3]
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
    return - vals.pot ./ 4π
end

function laplace3d_D_fmm3d(interface::DielectricInterface{P, Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 3, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
        sources[3, i] = point.panel_point.point[3]
        norms[1, i] = point.panel_point.normal[1]
        norms[2, i] = point.panel_point.normal[2]
        norms[3, i] = point.panel_point.normal[3]
    end

    f = charges -> _laplace3d_D_fmm3d(charges, sources, weights, norms, thresh)
    return LinearMap{Float64}(f, n_points, n_points)
end

function _laplace3d_D_trg_fmm3d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, targets::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64)
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

    vals = lfmm3d(thresh, sources, dipvecs = dipvecs, targets = targets, pgt = 1)
    return vals.pottarg ./ 4π
end

function laplace3d_D_trg_fmm3d(interface::DielectricInterface{P, Float64}, targets::Matrix{Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 3, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
        sources[3, i] = point.panel_point.point[3]
        norms[1, i] = point.panel_point.normal[1]
        norms[2, i] = point.panel_point.normal[2]
        norms[3, i] = point.panel_point.normal[3]
    end

    f = charges -> _laplace3d_D_trg_fmm3d(charges, sources, targets, weights, norms, thresh)
    return LinearMap{Float64}(f, size(targets, 2), n_points)
end

function _laplace3d_pottrg_fmm3d(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, targets::Matrix{Float64}, thresh::Float64)
    n = length(charges)
    @assert size(sources) == (3, n)
    @assert size(targets, 1) == 3
    @assert size(weights) == (n,)
    vals = lfmm3d(thresh, sources, charges = weights .* charges, targets = targets, pgt = 1)
    return vals.pottarg ./ 4π
end

function laplace3d_pottrg_fmm3d(interface::DielectricInterface{P, Float64}, targets::Matrix{Float64}, thresh::Float64) where {P <: AbstractPanel}
    n_points = num_points(interface)
    sources = zeros(Float64, 3, n_points)
    weights = zeros(Float64, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        weights[i] = point.panel_point.weight
        sources[1, i] = point.panel_point.point[1]
        sources[2, i] = point.panel_point.point[2]
        sources[3, i] = point.panel_point.point[3]
    end

    f = charges -> _laplace3d_pottrg_fmm3d(charges, sources, weights, targets, thresh)
    return LinearMap{Float64}(f, size(targets, 2), n_points)
end
