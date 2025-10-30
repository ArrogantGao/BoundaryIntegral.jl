using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, OMEinsum
using Test

@testset "coulomb" begin
    box = BI.box2d_adaptive_panels(8,16,2)
    DT = BI.laplace2d_DT(box)

    i = ones(BI.num_points(box))
    phi = ein"(ij, j), j -> i"(transpose(DT), BI.all_weights(box), i)

    @test norm(phi .- 0.5) < 0.02
end