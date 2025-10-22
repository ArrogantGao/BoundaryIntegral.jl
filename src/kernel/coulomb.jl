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