function dielectric_box2d(n_panels::Int, n_quad::Int,::Type{T} = Float64; adapt::Bool = false, n_adapt::Int = 3) where T
    box = adapt ? box2d_adaptive_panels(n_panels, n_quad, n_adapt, T) : box2d_uniform_panels(n_panels, n_quad, T)
    return box
end

function Lhs_dielectric_box2d(eps_box::T, box::Interface{T, 2}) where T
    D_transpose = laplace2d_DT(box)
    Lhs = I(num_points(box)) * 0.5 * (1.0 + eps_box) / (1.0 - eps_box) + ein"ij,j->ij"(D_transpose, all_weights(box))
    return Lhs
end

function Rhs_dielectric_box2d(eps_box::T, box::Interface{T, 2}, src::NTuple{2, T}) where T
    Rhs = - [laplace2d_doublelayer(src, point, panel.normal) for panel in box.panels for point in panel.points] ./ eps_box
    return Rhs
end

function dielectric_dbox2d(eps_box1::T, eps_box2::T, n_panels::Int, n_quad::Int, n_adapt::Int) where T
    i1, i2, i12 = dbox2d_adaptive_panels(n_panels, n_quad, n_adapt, T)
    di = [(i1, eps_box1, one(T)), (i2, eps_box2, one(T)), (i12, eps_box1, eps_box2)]
    return DielectricInterfaces(length(di), di)
end

function Lhs_dielectric_mbox2d(dbox::DielectricInterfaces{T, 2}) where T
    D_transpose = laplace2d_DT(dbox)
    Lhs = ein"ij,j->ij"(D_transpose, all_weights(dbox))
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

function Rhs_dielectric_mbox2d(dbox::DielectricInterfaces{T, 2}, src::NTuple{2, T}, eps_src::T) where T
    n_points = num_points(dbox)
    Rhs = zeros(T, n_points)
    for point in eachpoint(dbox)
        Rhs[point.global_idx] = - laplace2d_doublelayer(src, point.point, point.normal) / eps_src
    end
    return Rhs
end