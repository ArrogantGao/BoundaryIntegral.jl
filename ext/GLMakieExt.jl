module GLMakieExt

using GLMakie
using BoundaryIntegral
using BoundaryIntegral: Interface, DielectricInterfaces, Panel

import BoundaryIntegral: viz_2d_interfaces, viz_2d_dielectric_interfaces
import BoundaryIntegral: viz_3d_squares, viz_3d_interface

function viz_2d_interface!(ax::Axis, interface::Interface{T, 2}) where T
    t = 0.05

    for (i, panel) in enumerate(interface.panels)
        lines!(ax, [panel.points[1][1], panel.points[end][1]], [panel.points[1][2], panel.points[end][2]], color = :blue)

        norm = panel.normal
        for (x1, x2) in panel.points
            lines!(ax, [x1, x1 + t * norm[1]], [x2, x2 + t * norm[2]], color = :black, linewidth = 0.2)
        end
    end

end

function viz_2d_interface(interface::Interface{T, 2}) where T
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1])

    viz_2d_interface!(ax, interface)

    return fig
end

function viz_2d_interfaces(interfaces::Vector{Interface{T, 2}}) where T
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1])

    for interface in interfaces
        viz_2d_interface!(ax, interface)
    end

    return fig
end

function viz_2d_interfaces(interfaces_dict::Dict{Tuple{Int, Int}, Interface{T, 2}}) where T
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1])

    for (id1, id2) in keys(interfaces_dict)
        viz_2d_interface!(ax, interfaces_dict[(id1, id2)])
    end

    return fig
end

function viz_2d_dielectric_interfaces(dbox::DielectricInterfaces{T, 2}) where T
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1])

    for (interface, eps_in, eps_out) in dbox.interfaces
        viz_2d_interface!(ax, interface)
    end

    return fig
end

function viz_3d_squares(squares::Vector{NTuple{4, NTuple{3, T}}}) where T
    fig = Figure()
    ax = Axis3(fig[1, 1])
    viz_3d_squares!(ax, squares)
    return fig
end

function viz_3d_squares!(ax::Axis3, squares::Vector{NTuple{4, NTuple{3, T}}}) where T
    for square in squares
        lines!(ax, [square[1][1], square[2][1], square[3][1], square[4][1], square[1][1]], [square[1][2], square[2][2], square[3][2], square[4][2], square[1][2]], [square[1][3], square[2][3], square[3][3], square[4][3], square[1][3]], color = :blue)
    end
end

function viz_3d_interface(interface::Interface{T, 3}; show_normal::Bool = false) where T
    fig = Figure()
    ax = Axis3(fig[1, 1])
    viz_3d_panels!(ax, interface.panels, show_normal = show_normal)
    return fig
end

function viz_3d_panels!(ax::Axis3, panels::Vector{Panel{T, 3}}; show_normal::Bool = false) where T
    for panel in panels
        scatter!(ax, [p[1] for p in panel.points], [p[2] for p in panel.points], [p[3] for p in panel.points], color = :blue, markersize = 1.5)
        if show_normal
            nx, ny, nz = 0.2 .* panel.normal
            lines!(ax, [panel.points[1][1], panel.points[1][1] + nx], [panel.points[1][2], panel.points[1][2] + ny], [panel.points[1][3], panel.points[1][3] + nz], color = :black, linewidth = 0.4)
        end
    end
end

end