function dielectric_box2d(n_panels::Int, n_quad::Int, adapt::Bool = false, ::Type{T} = Float64; n_adapt::Int = 3) where T
    box = adapt ? box2d_adaptive_panels(n_panels, n_quad, n_adapt, T) : box2d_uniform_panels(n_panels, n_quad, T)
    return box
end

function Lhs_dielectric_box2d(eps_box::T, box::Surface{T, 2}) where T
    D_transpose = laplace2d_DT(box)
    Lhs = I(num_points(box)) * 0.5 * (1.0 + eps_box) / (1.0 - eps_box) + ein"ij,j->ij"(D_transpose, all_weights(box))
    return Lhs
end

function Rhs_dielectric_box2d(eps_box::T, box::Surface{T, 2}, src::NTuple{2, T}) where T
    Rhs = - [laplace2d_doublelayer(src, point, panel.normal) for panel in box.panels for point in panel.points] ./ eps_box
    return Rhs
end