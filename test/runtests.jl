using BoundaryIntegral
using LinearAlgebra
using Test

@testset "BoundaryIntegral.jl" begin

    # kernel functions
    include("kernel/laplace2d_fmm2d.jl")
end
