using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, OMEinsum
using Test

@testset "coulomb" begin
    box = BI.box2d_adaptive_panels(8,16,2)
    D = BI.laplace2d_D(box)

    i = ones(BI.num_points(box))
    phi = D * i

    @test norm(phi .- 0.5) < 0.02
end

@testset "fmm2d laplace2d_doublelayer" begin
    eps_box1 = 2.0
    eps_box2 = 3.0
    eps_box3 = 4.0
    
    rects = [BI.square(-1.0, -1.0), BI.square(0.0, -1.0), BI.square(-0.5, 0.0)]
    mbox = BI.dielectric_mbox2d([eps_box1, eps_box2, eps_box3], rects, 8, 16, 5)
    
    DT_direct = BI.laplace2d_DT(mbox)
    DT_fmm2d = BI.laplace2d_DT_fmm2d(mbox, 1e-14)

    x = rand(BI.num_points(mbox))
    @test norm(DT_direct * x - DT_fmm2d * x) < 1e-13
end