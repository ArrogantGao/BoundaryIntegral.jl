using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, OMEinsum, Krylov
using Test

@testset "dielectric_box" begin
    eps_box = 2.0
    box1 = BI.dielectric_box2d(8,16, Float64, adapt = false)
    box2 = BI.dielectric_box2d(8,16, Float64, adapt = true)

    lhs1 = BI.Lhs_dielectric_box2d(eps_box, box1)
    lhs2 = BI.Lhs_dielectric_box2d(eps_box, box2)
    rhs1 = BI.Rhs_dielectric_box2d(eps_box, box1, (0.1, 0.1))
    rhs2 = BI.Rhs_dielectric_box2d(eps_box, box2, (0.1, 0.1))

    x1 = BI.solve_lu(lhs1, rhs1)
    x2 = BI.solve_lu(lhs2, rhs2)

    @test norm(lhs1 * x1 - rhs1) < 1e-10
    @test norm(lhs2 * x2 - rhs2) < 1e-10

    g1 = BI.l2d_singlelayer_gi(box1, x1, 2.0, 32) + 1.0 / eps_box
    g2 = BI.l2d_singlelayer_gi(box2, x2, 2.0, 32) + 1.0 / eps_box

    @test isapprox(g1, 1.0, atol = 1e-4)
    @test isapprox(g2, 1.0, atol = 1e-4)
end

@testset "2 dielectric box" begin
    eps_box1 = 2.0
    eps_box2 = 3.0
    dbox = BI.dielectric_dbox2d(eps_box1, eps_box2, 8, 16, 2)

    lhs_direct = BI.Lhs_dielectric_mbox2d(dbox)
    lhs_fmm2d = BI.Lhs_dielectric_mbox2d_fmm2d(dbox, 1e-12)
    rhs = BI.Rhs_dielectric_mbox2d(dbox, (0.1, 0.1), eps_box2)

    x_direct = BI.solve_lu(lhs_direct, rhs)
    x_fmm2d = BI.solve_gmres(lhs_fmm2d, rhs)

    @test norm(lhs_direct * x_direct - rhs) < 1e-8
    @test norm(lhs_direct * x_fmm2d - rhs) < 1e-6

    g_direct = BI.l2d_singlelayer_gi(dbox, x_direct, 2.0, 32) + 1.0 / eps_box2
    g_fmm2d = BI.l2d_singlelayer_gi(dbox, x_fmm2d, 2.0, 32) + 1.0 / eps_box2

    @test isapprox(g_direct, 1.0, atol = 1e-4)
    @test isapprox(g_fmm2d, 1.0, atol = 1e-4)
end

@testset "3 dielectric box" begin
    eps_box1 = 2.0
    eps_box2 = 3.0
    eps_box3 = 4.0
    
    rects = [BI.square(-1.0, -1.0), BI.square(0.0, -1.0), BI.square(-0.5, 0.0)]
    mbox = BI.dielectric_mbox2d([eps_box1, eps_box2, eps_box3], rects, 8, 16, 5)

    lhs_direct = BI.Lhs_dielectric_mbox2d(mbox)
    lhs_fmm2d = BI.Lhs_dielectric_mbox2d_fmm2d(mbox, 1e-12)
    rhs = BI.Rhs_dielectric_mbox2d(mbox, (0.0, 0.5), eps_box3)

    x_direct = BI.solve_lu(lhs_direct, rhs)
    x_fmm2d = BI.solve_gmres(lhs_fmm2d, rhs)

    @test norm(lhs_direct * x_direct - rhs) < 1e-8
    @test norm(lhs_direct * x_fmm2d - rhs) < 1e-6

    g_direct = BI.l2d_singlelayer_gi(mbox, x_direct, 2.0, 32) + 1.0 / eps_box3
    g_fmm2d = BI.l2d_singlelayer_gi(mbox, x_fmm2d, 2.0, 32) + 1.0 / eps_box3

    @test isapprox(g_direct, 1.0, atol = 1e-4)
    @test isapprox(g_fmm2d, 1.0, atol = 1e-4)
end

@testset "3 dielectric box float" begin
    eps_box1 = 2.0f0
    eps_box2 = 3.0f0
    eps_box3 = 4.0f0
    
    rects = [BI.square(-1.0f0, -1.0f0), BI.square(0.0f0, -1.0f0), BI.square(-0.5f0, 0.0f0)]
    mbox = BI.dielectric_mbox2d([eps_box1, eps_box2, eps_box3], rects, 8, 16, 5)

    lhs_direct = BI.Lhs_dielectric_mbox2d(mbox)
    rhs = BI.Rhs_dielectric_mbox2d(mbox, (0.0f0, 0.5f0), eps_box3)

    lhs_ka = BI.Lhs_dielectric_mbox2d_ka(mbox)
    rhs_ka = BI.Rhs_dielectric_mbox2d_ka(mbox, (0.0f0, 0.5f0), eps_box3)

    x_direct = BI.solve_lu(lhs_direct, rhs)
    x_ka = Array(Krylov.gmres(lhs_ka, rhs_ka)[1])

    @test norm(lhs_direct * x_direct - rhs) < 1e-5
    @test norm(lhs_direct * x_ka - rhs) < 1e-3

    g_direct = BI.l2d_singlelayer_gi(mbox, x_direct, 2.0f0, 32) + 1.0 / eps_box3
    g_ka = BI.l2d_singlelayer_gi(mbox, x_ka, 2.0f0, 32) + 1.0 / eps_box3

    @test isapprox(g_direct, 1.0, atol = 1e-4)
    @test isapprox(g_ka, 1.0, atol = 1e-3)
end