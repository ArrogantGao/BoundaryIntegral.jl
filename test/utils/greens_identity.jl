using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, OMEinsum
using Test

@testset "greens_identity" begin
    gi_sphere = BI.l2d_sphere_gi(1.0, 32, (0.1, 0.1))
    @test isapprox(gi_sphere, 1.0, atol = 1e-10)
end