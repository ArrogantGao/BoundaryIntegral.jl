function l2d_point_source(surface::Surface{T, 2}, src::NTuple{2, T}) where T
    t = 0.0
    for panel in surface.panels
        for (point, weight) in zip(panel.points, panel.weights)
            t += laplace2d_doublelayer(src, point, panel.normal) * weight
        end
    end
    return t
end