module GLMakieExt

using GLMakie
using BoundaryIntegral
using BoundaryIntegral: AbstractPanel, DielectricInterface

import BoundaryIntegral: viz_2d, viz_3d

function viz_2d!(ax::Axis, interface::DielectricInterface{P, T}; show_normals::Bool = true) where {P <: AbstractPanel, T}
    t = 0.05

    for panel in interface.panels
        a, b = panel.corners
        lines!(ax, [a[1], b[1]], [a[2], b[2]], color = :blue)

        show_normals || continue
        norm = panel.normal
        for (x1, x2) in panel.points
            lines!(ax, [x1, x1 + t * norm[1]], [x2, x2 + t * norm[2]], color = :black, linewidth = 0.2)
        end
    end
end

function viz_2d(interface::DielectricInterface{P, T}; show_normals::Bool = true, size = (600, 600)) where {P <: AbstractPanel, T}
    fig = Figure(size = size)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    viz_2d!(ax, interface; show_normals = show_normals)

    return fig
end

function viz_3d!(ax::Axis3, interface::DielectricInterface{P, T}; show_normals::Bool = false, show_points::Bool = true) where {P <: AbstractPanel, T}
    t = 0.2
    for panel in interface.panels
        a, b, c, d = panel.corners
        lines!(ax, [a[1], b[1], c[1], d[1], a[1]], [a[2], b[2], c[2], d[2], a[2]], [a[3], b[3], c[3], d[3], a[3]], color = :blue, linewidth = 0.6)

        if show_points
            xs = [p[1] for p in panel.points]
            ys = [p[2] for p in panel.points]
            zs = [p[3] for p in panel.points]
            scatter!(ax, xs, ys, zs, color = :red, markersize = 3)
        end

        show_normals || continue
        nx, ny, nz = t .* panel.normal
        xs = [p[1] for p in panel.points]
        ys = [p[2] for p in panel.points]
        zs = [p[3] for p in panel.points]
        lines!(ax, xs .+ nx, ys .+ ny, zs .+ nz, color = :black, linewidth = 0.4)
    end
end

function viz_3d(interface::DielectricInterface{P, T}; show_normals::Bool = false, show_points::Bool = true, size = (700, 600)) where {P <: AbstractPanel, T}
    fig = Figure(size = size)
    ax = Axis3(fig[1, 1], aspect = :data)
    viz_3d!(ax, interface; show_normals = show_normals, show_points = show_points)
    return fig
end

function viz_3d(interfaces::Vector{<:DielectricInterface{P, T}}; show_normals::Bool = false, show_points::Bool = true, size = (700, 600)) where {P <: AbstractPanel, T}
    fig = Figure(size = size)
    ax = Axis3(fig[1, 1], aspect = :data)
    for interface in interfaces
        viz_3d!(ax, interface; show_normals = show_normals, show_points = show_points)
    end
    return fig
end

end
