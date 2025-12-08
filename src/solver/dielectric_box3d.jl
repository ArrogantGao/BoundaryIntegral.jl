function dielectric_box3d(eps_box::T, eps_out::T, n_boxes::Int, n_quad::Int, reduce_quad::Bool, n_adapt_edge::Int, n_adapt_corner::Int, ::Type{T} = Float64) where T
    box = single_box3d(n_boxes, n_quad, reduce_quad, n_adapt_edge, n_adapt_corner, T)
    return DielectricInterfaces(1, [(box, eps_box, eps_out)])
end

function dielectric_double_box3d(eps_box1::T, eps_box2::T, eps_out::T, n_quad::Int, n_adapt_edge::Int, n_adapt_corner::Int, ::Type{T} = Float64) where T

    ns, ws = gausslegendre(n_quad)

    interfaces = Vector{Tuple{Interface{T, 3}, T, T}}()

    # box1 at left, (0, -1, 0) -> (1, 0, 1), box2 at right, (0, 0, 0) -> (1, 1, 1)
    square_1 = [
        ((1.0, -1.0, 1.0), (0.0, -1.0, 1.0), (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), 
        ((1.0, -1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, -1.0, 0.0)), 
        ((1.0, 0.0, 0.0), (1.0, -1.0, 0.0), (1.0, -1.0, 1.0), (1.0, 0.0, 1.0)), 
        ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, -1.0, 1.0), (0.0, -1.0, 0.0)), 
        ((0.0, -1.0, 0.0), (0.0, -1.0, 1.0), (1.0, -1.0, 1.0), (1.0, -1.0, 0.0)), 
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 0.0, 1.0), (0.0, 0.0, 1.0))]
    square_2 = [
        ((1.0, 0.0, 1.0), (0.0, 0.0, 1.0), (0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), 
        ((1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), 
        ((1.0, 1.0, 0.0), (1.0, 0.0, 0.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)), 
        ((0.0, 1.0, 0.0), (0.0, 1.0, 1.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0)), 
        ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 0.0, 0.0)), 
        ((0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0))]
    
    normals = [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, - 1.0, 0.0), (0.0, 1.0, 0.0)]

    surf_1 = [(square_1[i], normals[i]) for i in (1:5)]
    surf_2 = [(square_2[i], normals[i]) for i in [1, 2, 3, 4, 6]]
    surf_3 = (square_1[6], normals[5])

    panels_1 = Vector{Panel{T, 3}}()
    for surf in surf_1
        vertices, normal = surf
        new_panels = square_surface_adaptive_panels(vertices..., ns, ws, normal, (true, true, true, true), (true, true, true, true), n_adapt_edge, n_adapt_corner)
        append!(panels_1, new_panels)
    end
    interface_1 = Interface(length(panels_1), panels_1)

    panels_2 = Vector{Panel{T, 3}}()
    for surf in surf_2
        vertices, normal = surf
        new_panels = square_surface_adaptive_panels(vertices..., ns, ws, normal, (true, true, true, true), (true, true, true, true), n_adapt_edge, n_adapt_corner)
        append!(panels_2, new_panels)
    end
    interface_2 = Interface(length(panels_2), panels_2)

    panels_3 = square_surface_adaptive_panels(surf_3[1][1], surf_3[1][2], surf_3[1][3], surf_3[1][4], ns, ws, surf_3[2], (true, true, true, true), (true, true, true, true), n_adapt_edge, n_adapt_corner)
    interface_3 = Interface(length(panels_3), panels_3)

    interfaces = [(interface_1, eps_box1, eps_out), (interface_2, eps_box2, eps_out), (interface_3, eps_box2, eps_box1)]

    return DielectricInterfaces(length(interfaces), interfaces)
end

function Lhs_dielectric_mbox3d_direct(dbox::DielectricInterfaces{T, 3}) where T
    D_transpose = laplace3d_DT(dbox)
    Lhs = D_transpose
    offset = 0
    for (interface, eps_in, eps_out) in dbox.interfaces
        t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
        for i in 1:num_points(interface)
            Lhs[offset + i, offset + i] += t
        end
        offset += num_points(interface)
    end
    return Lhs
end

function Lhs_dielectric_mbox3d_fmm3d(dbox::DielectricInterfaces{T, 3}, tol::Float64 = 1e-12) where T
    D_transpose = laplace3d_DT_fmm3d(dbox, tol)

    function g(x)
        Dx = D_transpose * x

        offset = 0
        for (interface, eps_in, eps_out) in dbox.interfaces
            t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
            for i in 1:num_points(interface)
                Dx[offset + i] += t * x[offset + i]
            end
            offset += num_points(interface)
        end

        return Dx
    end

    Lhs = LinearMap{T}(g, num_points(dbox), num_points(dbox))

    return Lhs
end

function Rhs_dielectric_mbox3d(dbox::DielectricInterfaces{T, 3}, src::NTuple{3, T}, eps_src::T) where T
    n_points = num_points(dbox)
    Rhs = zeros(T, n_points)
    for point in eachpoint(dbox)
        Rhs[point.global_idx] = - laplace3d_doublelayer(src, point.point, point.normal) / eps_src
    end
    return Rhs
end