function laplace3d_singlelayer(src::NTuple{3, T}, trg::NTuple{3, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    inv_r = one(T) / r

    return inv_r / 4π
end

function laplace2d_singlelayer(src::NTuple{2, T}, trg::NTuple{2, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    return - log(r) / 2π
end

function laplace2d_singlelayer_interface(interface::Interface{T, 2}, sigma::Vector{T}, trg::NTuple{2, T}) where T
    t = 0.0
    i = 0
    for panel in interface.panels
        for (point, weight) in zip(panel.points, panel.weights)
            i += 1
            t += laplace2d_singlelayer(point, trg) * weight * sigma[i]
        end
    end
    return t
end

function laplace2d_singlelayer_interface(dielectric_interfaces::DielectricInterfaces{T, 2}, sigma::Vector{T}, trg::NTuple{2, T}) where T
    t = 0.0
    for (i, pointi) in enumerate(eachpoint(dielectric_interfaces))
        t += laplace2d_singlelayer(pointi.point, trg) * pointi.weight * sigma[i]
    end
    return t
end


function laplace3d_doublelayer(src::NTuple{3, T}, trg::NTuple{3, T}, norm::NTuple{3, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    inv_r = one(T) / r

    return dot(norm, inv_r^3 .* (trg .- src))  / 4π
end

function laplace2d_doublelayer(src::NTuple{2, T}, trg::NTuple{2, T}, norm::NTuple{2, T}) where T
    r2 = sum((src .- trg).^2)
    inv_r2 = one(T) / r2

    return dot(norm, inv_r2 .* (trg .- src)) / 2π
end

function laplace2d_D(interface::Interface{T, 2}) where{T}
    n_points = num_points(interface)
    weights = all_weights(interface)

    D = zeros(T, n_points, n_points)
    i_offset = 0;
    for paneli in interface.panels
        for (i, pointi) in enumerate(paneli.points)
            i_offset += 1
            j_offset = 0;
            for panelj in interface.panels
                for (j, pointj) in enumerate(panelj.points)
                    j_offset += 1
                    i_offset == j_offset && continue
                    D[j_offset, i_offset] = laplace2d_doublelayer(pointj, pointi, paneli.normal)
                end
            end
        end
    end

    return D * diagm(weights)
end

function laplace2d_DT(interface::Interface{T, 2}) where{T}
    n_points = num_points(interface)
    weights = all_weights(interface)

    DT = zeros(T, n_points, n_points)
    i_offset = 0;
    for paneli in interface.panels
        for (i, pointi) in enumerate(paneli.points)
            i_offset += 1
            j_offset = 0;
            for panelj in interface.panels
                for (j, pointj) in enumerate(panelj.points)
                    j_offset += 1
                    i_offset == j_offset && continue
                    DT[i_offset, j_offset] = laplace2d_doublelayer(pointj, pointi, paneli.normal)
                end
            end
        end
    end

    return DT * diagm(weights)
end

function laplace2d_DT(dielectric_interfaces::DielectricInterfaces{T, 2}) where T
    n_points = num_points(dielectric_interfaces)
    weights = all_weights(dielectric_interfaces)

    DT = zeros(T, n_points, n_points)
    for (i, pointi) in enumerate(eachpoint(dielectric_interfaces))
        for (j, pointj) in enumerate(eachpoint(dielectric_interfaces))
            i == j && continue
            DT[i, j] = laplace2d_doublelayer(pointj.point, pointi.point, pointi.normal)
        end
    end

    return DT * diagm(weights)
end

# matrix free via fmm2d, return a linear map
# potential: log(r), gradient: 1/r
function laplace2d_DT_fmm2d_kernel(charges::AbstractVector{Float64}, sources::Matrix{Float64}, weights::Vector{Float64}, norms::Matrix{Float64}, thresh::Float64)
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


function laplace2d_DT_fmm2d(dielectric_interfaces::DielectricInterfaces{T, 2}, thresh::Float64) where T
    n_points = num_points(dielectric_interfaces)
    sources = zeros(Float64, 2, n_points)
    weights = zeros(Float64, n_points)
    norms = zeros(Float64, 2, n_points)
    for (i, point) in enumerate(eachpoint(dielectric_interfaces))
        weights[i] = point.weight
        sources[1, i] = point.point[1]
        sources[2, i] = point.point[2]
        norms[1, i] = point.normal[1]
        norms[2, i] = point.normal[2]
    end

    f = charges -> laplace2d_DT_fmm2d_kernel(charges, sources, weights, norms, thresh)
    return LinearMap(f, n_points, n_points)
end