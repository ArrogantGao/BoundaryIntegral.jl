function viz_2d_surface(surface::Surface{T, 2}) where T
    p = plot()
    for panel in surface.panels
        xs = T[]
        ys = T[]
        weights = T[]
        for (point, weight) in zip(panel.points, panel.weights)
            push!(xs, point[1])
            push!(ys, point[2])
            push!(weights, weight)
        end
        scatter!(xs, ys, label=nothing, markersize=1)
    end
    plot!(p, aspect_ratio=:equal)
    return p
end