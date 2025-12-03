function dielectric_box3d(eps_box::T, n_quad::Int, n_adapt_edge::Int, n_adapt_corner::Int, ::Type{T} = Float64) where T
    box = single_box3d(n_quad, n_adapt_edge, n_adapt_corner, T)
    return DielectricInterfaces(1, [(box, eps_box, one(T))])
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