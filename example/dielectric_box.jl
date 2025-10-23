using BoundaryIntegral
import BoundaryIntegral as BI
using OMEinsum, LinearAlgebra

eps_box = 2.0
src = (0.0, 0.0)

function solve_dielectric_box(eps_box::Float64, src::NTuple{2, Float64}, n_panels::Int, n_quad::Int)

    b = BI.box2d_uniform_panels(n_panels, n_quad)
    DT = BI.laplace2d_DT(b);

    L = ein"ij,j->ij"(I(BI.num_points(b)) * 0.5 * (1.0 + eps_box) / (eps_box - 1.0) + DT, BI.all_weights(b));
    R = - [BI.laplace2d_doublelayer(src, point, panel.normal) for panel in b.panels for point in panel.points] ./ eps_box;

    sigma = inv(L) * R;
    V = BI.l2d_singlelayer_gi(b, sigma, 2.0, 32) + 1.0 / eps_box

    return V
end
