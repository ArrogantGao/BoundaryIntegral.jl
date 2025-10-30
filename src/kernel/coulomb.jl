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

function laplace2d_DT(interface::Interface{T, 2}) where{T}
    n_points = num_points(interface)

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
    return DT
end