using BoundaryIntegral
using FastGaussQuadrature
using LinearAlgebra
import BoundaryIntegral as BI
using Test

@testset "quad_order3d" begin
    ns, ws = gausslegendre(4)
    ns = Float64.(ns)
    ws = Float64.(ws)
    a = (-0.5, -0.5, 0.0)
    b = (0.5, -0.5, 0.0)
    c = (0.5, 0.5, 0.0)
    d = (-0.5, 0.5, 0.0)
    normal = (0.0, 0.0, 1.0)
    panel = BI.rect_panel3d_discretize(a, b, c, d, ns, ws, normal)
    trg = (0.1, -0.2, 0.25)

    function gaussian_integral_direct(nq::Int)
        xs, ws = gausslegendre(nq)
        xs = Float64.(xs)
        ws = Float64.(ws)
        cc = (a .+ b .+ c .+ d) ./ 4
        Lx = norm(b .- a)
        Ly = norm(d .- a)
        alpha = 2.0
        acc = 0.0
        for i in 1:nq
            for j in 1:nq
                p = cc .+ (b .- a) .* (xs[i] / 2) .+ (d .- a) .* (xs[j] / 2)
                rho = exp(-alpha * (p[1]^2 + p[2]^2 + p[3]^2))
                acc += ws[i] * ws[j] * rho * BI.laplace3d_pot(p, trg) * Lx * Ly / 4
            end
        end
        return acc
    end

    function gaussian_integral_upsampled(nq::Int)
        xs, ws = gausslegendre(nq)
        xs = Float64.(xs)
        ws = Float64.(ws)
        cc = (a .+ b .+ c .+ d) ./ 4
        Lx = norm(b .- a)
        Ly = norm(d .- a)
        alpha = 2.0

        ns0 = panel.gl_xs
        ws0 = panel.gl_ws
        n0 = length(ns0)
        rho0 = zeros(Float64, n0 * n0)
        idx = 1
        for j in 1:n0
            for i in 1:n0
                p = cc .+ (b .- a) .* (ns0[i] / 2) .+ (d .- a) .* (ns0[j] / 2)
                rho0[idx] = exp(-alpha * (p[1]^2 + p[2]^2 + p[3]^2))
                idx += 1
            end
        end

        E = BI.interp_matrix_2d_gl_tensor(ns0, ws0, ns0, ws0, xs, xs)
        rho_up = E * rho0

        acc = 0.0
        idx = 1
        for j in 1:nq
            for i in 1:nq
                p = cc .+ (b .- a) .* (xs[i] / 2) .+ (d .- a) .* (xs[j] / 2)
                acc += ws[i] * ws[j] * rho_up[idx] * BI.laplace3d_pot(p, trg) * Lx * Ly / 4
                idx += 1
            end
        end
        return acc
    end

    n_quad = panel.n_quad
    atol = 1e-6
    max_order = 16
    order = BI.check_quad_order3d(panel, trg, atol, max_order)
    @test order >= n_quad
    ref_val = gaussian_integral_direct(max_order)
    test_val = gaussian_integral_upsampled(order)
    @test abs(test_val - ref_val) <= 1e-4
end
