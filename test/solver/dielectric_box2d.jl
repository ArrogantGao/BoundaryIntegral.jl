using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using Test

@testset "dielectric_box" begin
    eps_box = 5.0

    for Lx in [1.0, 2.05]
        for Ly in [1.0, 3.03]
            box = BI.single_dielectric_box2d(8, 0.2, 0.05, 5.0, 1.0, Float64; Lx = Lx, Ly = Ly)
            lhs = BI.Lhs_dielectric_box2d(box)
            lhs_fmm2d = BI.Lhs_dielectric_box2d_fmm2d(box, 1e-12)
            rhs = BI.Rhs_dielectric_box2d(box, BI.PointSource((0.1, 0.1), 1.0), eps_box)
            ws = BI.all_weights(box)

            x = BI.solve_lu(lhs, rhs)
            @test norm(lhs * x - rhs) < 1e-10

            total_flux = dot(ws, x)
            @test isapprox(total_flux + 1.0 / eps_box, 1.0, atol = 1e-3)

            x_gmres = BI.solve_gmres(lhs_fmm2d, rhs,1e-12, 1e-12)
            @test norm(lhs_fmm2d * x_gmres - rhs) < 1e-10

            total_flux_gmres = dot(ws, x_gmres)
            @test isapprox(total_flux_gmres + 1.0 / eps_box, 1.0, atol = 1e-3)
        end
    end
end