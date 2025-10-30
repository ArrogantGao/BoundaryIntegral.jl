using BoundaryIntegral
using LinearAlgebra, OMEinsum
using Test

@testset "BoundaryIntegral.jl" begin

    # kernel functions
    include("kernel/coulomb.jl")

    # geometries
    include("geometry/single_box.jl")

    # utilities
    include("utils/greens_identity.jl")

    # solvers
    include("solver/dielectric_box.jl")
end
