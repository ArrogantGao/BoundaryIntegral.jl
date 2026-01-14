module GLMakieExt

using GLMakie
using BoundaryIntegral
using BoundaryIntegral: AbstractPanel, DielectricInterface

import BoundaryIntegral: viz_2d

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

# function viz_3d_squares(squares::Vector{NTuple{4, NTuple{3, T}}}) where T
#     fig = Figure()
#     ax = Axis3(fig[1, 1])
#     viz_3d_squares!(ax, squares)
#     return fig
# end

# function viz_3d_squares!(ax::Axis3, squares::Vector{NTuple{4, NTuple{3, T}}}) where T
#     for square in squares
#         lines!(ax, [square[1][1], square[2][1], square[3][1], square[4][1], square[1][1]], [square[1][2], square[2][2], square[3][2], square[4][2], square[1][2]], [square[1][3], square[2][3], square[3][3], square[4][3], square[1][3]], color = :blue)
#     end
# end

# function viz_3d_interface(interface::Interface{T, 3}; show_normal::Bool = false) where T
#     fig = Figure()
#     ax = Axis3(fig[1, 1])
#     viz_3d_panels!(ax, interface.panels, show_normal = show_normal)
#     return fig
# end

# function viz_3d_panels!(ax::Axis3, panels::Vector{Panel{T, 3}}; show_normal::Bool = false) where T
#     for panel in panels
#         xs = [p[1] for p in panel.points]
#         ys = [p[2] for p in panel.points]
#         zs = [p[3] for p in panel.points]
#         scatter!(ax, xs, ys, zs, color = :blue, markersize = 1.5)
#         x_max = maximum(xs)
#         x_min = minimum(xs)
#         y_max = maximum(ys)
#         y_min = minimum(ys)
#         z_max = maximum(zs)
#         z_min = minimum(zs)

#         # link the corners with lines
#         a, b, c, d = panel.corners
#         lines!(ax, [a[1], b[1], c[1], d[1], a[1]], [a[2], b[2], c[2], d[2], a[2]], [a[3], b[3], c[3], d[3], a[3]], color = :red, linewidth = 0.4)

#         if show_normal
#             nx, ny, nz = 0.2 .* panel.normal
#             lines!(ax, xs .+ nx, ys .+ ny, zs .+ nz, color = :black, linewidth = 0.4)
#         end
#     end
# end

# function viz_3d_dielectric_interfaces(dbox::DielectricInterfaces{T, 3}; show_normal::Bool = false) where T
#     fig = Figure()
#     ax = Axis3(fig[1, 1], aspect = :data)
#     for (interface, eps_in, eps_out) in dbox.interfaces
#         viz_3d_panels!(ax, interface.panels, show_normal = show_normal)
#     end
#     return fig
# end

end
