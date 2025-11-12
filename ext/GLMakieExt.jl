module GLMakieExt

using GLMakie
using BoundaryIntegral

import BoundaryIntegral: viz_2d_interface!, viz_2d_interfaces, viz_2d_interfaces, viz_2d_dielectric_interfaces

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

end