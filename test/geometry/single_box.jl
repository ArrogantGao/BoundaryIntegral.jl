using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, OMEinsum
using Test

@testset "single_box" begin
    box1 = BI.box2d_uniform_panels(4,16)
    @test BI.num_points(box1) == 4 * 4 * 16
    @test box1.n == 4 * 4

    box2 = BI.box2d_adaptive_panels(4,16,2)
    @test box2.n == 6 * 4
    @test BI.num_points(box2) == 6 * 4 * 16

    box3 = BI.box2d_adaptive_panels(4,16,3)
    @test box3.n == 8 * 4
    @test BI.num_points(box3) == 8 * 4 * 16
end