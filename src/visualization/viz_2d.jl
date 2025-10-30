function viz_2d_surface(surface::Surface{T, 2}) where T
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1])

    xs = T[]
    ys = T[]
    weights = T[]
    colors = T[]

    for (i, panel) in enumerate(surface.panels)
        for (point, weight) in zip(panel.points, panel.weights)
            push!(xs, point[1])
            push!(ys, point[2])
            push!(weights, weight)
            push!(colors, i * 10)
        end
    end

    points = Point2f.(xs, ys)
    scatter!(ax, points, label=nothing, markersize=6, strokewidth = 0, color = colors)

    return fig
end