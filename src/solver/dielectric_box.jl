function solve_dielectric_box(eps_box::T, src::NTuple{2, T}, n_panels::Int, n_quad::Int, radius::T, n_sphere::Int, adapt::Bool = false; n_adapt::Int = 3) where T

    b = adapt ? box2d_adaptive_panels(n_panels, n_quad, n_adapt) : box2d_uniform_panels(n_panels, n_quad)
    DT = laplace2d_DT(b);

    L = I(num_points(b)) * 0.5 * (1.0 + eps_box) / (1.0 - eps_box) + ein"ij,j->ij"(DT, all_weights(b));
    R = - [laplace2d_doublelayer(src, point, panel.normal) for panel in b.panels for point in panel.points] ./ eps_box;

    sigma = inv(L) * R;
    V = l2d_singlelayer_gi(b, sigma, radius, n_sphere) + 1.0 / eps_box

    return V
end