using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Test

@testset "multi_box2d" begin
    rects = [BI.square(0.0, 0.0), BI.square(1.0, 0.0)]
    eps_by_id = Dict(0 => 1.0, 1 => 2.0, 2 => 3.0)

    interfaces = BI.multi_box2d(8, 0.6, 0.3, rects, Float64; eps_by_id = eps_by_id)

    @test sort(collect(keys(interfaces))) == sort([(1, 0), (2, 0), (2, 1)])

    shared = interfaces[(2, 1)]
    @test !isempty(shared.panels)
    @test all(==(3.0), shared.eps_in)
    @test all(==(2.0), shared.eps_out)

    for panel in shared.panels
        @test isapprox(panel.normal[1], -1.0; atol = 1e-12)
        @test isapprox(panel.normal[2], 0.0; atol = 1e-12)
    end
end
