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
            DT[i, j] = laplace3d_doublelayer(pointj.point, pointi.point, pointi.normal)
        end
    end
    return DT * diagm(weights)
end