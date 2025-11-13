@kernel function laplace2d_DT_kernel!(@Const(sources), @Const(charges), @Const(weights), @Const(norms), gradn)

    i = @index(Global, Linear)

    # creating a temporary sum variable for matrix multiplication
    # i is the target point, j is the source point
    tmp_sum = zero(eltype(gradn))
    @inbounds x_i = sources[1, i]
    @inbounds y_i = sources[2, i]
    @inbounds normx_i = norms[1, i]
    @inbounds normy_i = norms[2, i]
    for j in 1:size(sources, 2)
        i == j && continue

        @inbounds x_j = sources[1, j]
        @inbounds y_j = sources[2, j]

        @inbounds tmp_sum += charges[j] * weights[j] * (normx_i * (x_i - x_j) + normy_i * (y_i - y_j)) / (2π * ((x_j - x_i)^2 + (y_j - y_i)^2))
    end
   
    @inbounds gradn[i] = tmp_sum
end

function _laplace2d_DT_ka(sources, charges, weights, norms)
    n_points = size(sources, 2)
    gradn = KernelAbstractions.zeros(backend, eltype(sources), n_points)
    # KernelAbstractions.@launch(backend, KernelAbstractions.Threads(), n_points, laplace2d_DT_kernel!(sources, charges, weights, norms, gradn))

    laplace2d_DT_kernel!(backend, 256, n_points)(sources, charges, weights, norms, gradn)
    KernelAbstractions.synchronize(backend)
    return gradn
end

function laplace2d_DT_ka(dielectric_interfaces::DielectricInterfaces{T, 2}) where T
    n_points = num_points(dielectric_interfaces)
    sources = zeros(T, 2, n_points)
    weights = zeros(T, n_points)
    norms = zeros(T, 2, n_points)

    for (i, point) in enumerate(eachpoint(dielectric_interfaces))
        weights[i] = point.weight
        sources[1, i] = point.point[1]
        sources[2, i] = point.point[2]
        norms[1, i] = point.normal[1]
        norms[2, i] = point.normal[2]
    end

    sources = adapt(backend, sources)
    weights = adapt(backend, weights)
    norms = adapt(backend, norms)

    f = charges -> _laplace2d_DT_ka(sources, charges, weights, norms)

    return LinearMap{T}(f, num_points(dielectric_interfaces), num_points(dielectric_interfaces))
end