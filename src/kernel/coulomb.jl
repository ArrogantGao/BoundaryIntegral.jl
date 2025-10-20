function coulomb_En(src::NTuple{3, T}, trg::NTuple{3, T}, norm::NTuple{3, T}) where T
    r2 = sum((src .- trg).^2)
    r = sqrt(r2)
    inv_r = one(T) / r

    En = dot(norm, inv_r^3 .* (trg .- src))  / 4π

    return En
end