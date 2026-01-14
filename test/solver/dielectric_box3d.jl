using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using Random
using Test

@testset "dielectric_box3d" begin
    eps_box = 4.0
    interface = BI.single_dielectric_box3d(1.2, 0.8, 0.6, 4, 0.4, 0.2, eps_box, 1.0, Float64)

    lhs = BI.Lhs_dielectric_box3d(interface)
    lhs_fmm3d = BI.Lhs_dielectric_box3d_fmm3d(interface, 1e-12)
    rhs = BI.Rhs_dielectric_box3d(interface, BI.PointSource((0.1, 0.1, 0.1), 1.0), eps_box)
    ws = BI.all_weights(interface)

    x = BI.solve_lu(lhs, rhs)
    @test norm(lhs * x - rhs) < 1e-10

    x_trial = randn(BI.num_points(interface))
    @test norm(lhs * x_trial - lhs_fmm3d * x_trial) < 1e-8

    x_gmres = BI.solve_gmres(lhs_fmm3d, rhs, 1e-12, 1e-12)
    @test norm(lhs_fmm3d * x_gmres - rhs) < 1e-10

    total_flux = dot(ws, x)
    @test isapprox(total_flux + 1.0 / eps_box, 1.0, atol = 1e-1)

    total_flux_gmres = dot(ws, x_gmres)
    @test isapprox(total_flux_gmres + 1.0 / eps_box, 1.0, atol = 1e-1)
end
