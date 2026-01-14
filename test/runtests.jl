using BoundaryIntegral
using LinearAlgebra
using Test

@testset "BoundaryIntegral.jl" begin


    # utilities
    include("utils/linear_algebra.jl")
    include("utils/barycentric.jl")

    # kernel functions
    include("kernel/laplace2d_fmm2d.jl")
end
