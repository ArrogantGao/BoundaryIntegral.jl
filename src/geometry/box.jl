function uniform_box(n::Int, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n)
    t1 = one(T)
    t0 = zero(T)

    locations = Vector{NTuple{3, T}}()
    norms = Vector{NTuple{3, T}}()
    weights = Vector{T}()

    # generate locations, norms and weights on each surfaces
    for norm in [
        (t1, t0, t0), (-t1, t0, t0),
        (t0, t1, t0), (t0, -t1, t0),
        (t0, t0, t1), (t0, t0, -t1)]
        for i in 1:n, j in 1:n
            w = ws[i] * ws[j]
            n1, n2 = ns[i], ns[j]

            if norm[1] != t0
                loc = (norm[1], n1, n2)
            elseif norm[2] != t0
                loc = (n1, norm[2], n2)
            else
                loc = (n1, n2, norm[3])
            end

            push!(locations, loc)
            push!(norms, norm)
            push!(weights, w)
        end
    end

    return Surface(6 * n^2, locations, norms, weights)
end