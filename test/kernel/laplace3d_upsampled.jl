using BoundaryIntegral
using FastGaussQuadrature
import BoundaryIntegral as BI
using LinearAlgebra
using Test

@testset "laplace3d DT panel upsampled" begin
    ns, ws = gausslegendre(6)
    ns = Float64.(ns)
    ws = Float64.(ws)

    a = (-0.5, -0.5, 0.0)
    b = (0.5, -0.5, 0.0)
    c = (0.5, 0.5, 0.0)
    d = (-0.5, 0.5, 0.0)
    normal = (0.0, 0.0, 1.0)
    panel_src = BI.rect_panel3d_discretize(a, b, c, d, ns, ws, normal)

    a2 = (-0.6, -0.6, 0.4)
    b2 = (0.6, -0.6, 0.4)
    c2 = (0.6, 0.6, 0.4)
    d2 = (-0.6, 0.6, 0.4)
    panel_trg = BI.rect_panel3d_discretize(a2, b2, c2, d2, ns, ws, normal)

    function gaussian_density(points)
        rho = Vector{Float64}(undef, length(points))
        alpha = 0.5
        for i in 1:length(points)
            p = points[i]
            rho[i] = exp(-alpha * (p[1]^2 + p[2]^2 + p[3]^2))
        end
        return rho
    end

    n_ref = 48
    ns_ref, ws_ref = gausslegendre(n_ref)
    ns_ref = Float64.(ns_ref)
    ws_ref = Float64.(ws_ref)
    panel_ref = BI.rect_panel3d_discretize(a, b, c, d, ns_ref, ws_ref, normal)
    DT_ref = BI.laplace3d_DT_panel(panel_ref, panel_trg)
    rho_ref = gaussian_density(panel_ref.points)
    ref_val = DT_ref * rho_ref

    rho_src = gaussian_density(panel_src.points)
    trg_point = panel_trg.points[1]

    for atol in (1e-4, 1e-6, 1e-8)
        n_up = BI.check_quad_order3d(panel_src, trg_point, atol, 24)
        DT_up = BI.laplace3d_DT_panel_upsampled(panel_src, panel_trg, n_up)

        ns_up, ws_up = gausslegendre(n_up)
        ns_up = Float64.(ns_up)
        ws_up = Float64.(ws_up)
        panel_up = BI.rect_panel3d_discretize(a, b, c, d, ns_up, ws_up, normal)
        DT_direct = BI.laplace3d_DT_panel(panel_up, panel_trg)
        rho_up = gaussian_density(panel_up.points)

        err_up = norm(DT_up * rho_src - ref_val, Inf)
        err_direct = norm(DT_direct * rho_up - ref_val, Inf)

        @test err_up <= 40 * atol
        @test err_direct <= 40 * atol
    end
end
