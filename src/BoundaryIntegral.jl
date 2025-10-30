module BoundaryIntegral

using LinearAlgebra, OMEinsum
using FastGaussQuadrature

using GLMakie

#core types
export Interface

# kernel functions
export laplace3d_doublelayer, laplace2d_doublelayer, laplace3d_singlelayer, laplace2d_singlelayer

# geometries
export uniform_box3d, uniform_box2d

# linear algebra
export solve_lu

include("types.jl")

# kernel functions
include("kernel/coulomb.jl")

# geometries
include("geometry/single_box.jl")
include("geometry/multi_box.jl")

# utilities
include("utils/greens_identity.jl")
include("utils/linear_algebra.jl")

# solvers
include("solver/dielectric_box.jl")

# visualization
include("visualization/viz_2d.jl")

end
