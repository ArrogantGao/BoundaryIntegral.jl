module BoundaryIntegral

using LinearAlgebra
using FastGaussQuadrature
using Krylov, LinearMaps, Roots

using FMM2D, FMM3D

#core types
export FlatPanel
export DielectricInterface

export PointSource

# kernel functions
export laplace2d_doublelayer, laplace2d_singlelayer

# geometries
export uniform_box3d, uniform_box2d

# linear algebra
export solve_lu, solve_gmres

# visualization
export viz_2d

# core types
include("core/panels.jl")
include("core/sources.jl")

# kernel functions
include("kernel/laplace2d.jl")
include("kernel/laplace3d.jl")

# # geometries
include("shape/box2d.jl")
# include("shape/box3d.jl")

# # utilities
include("utils/linear_algebra.jl")
# include("utils/corner_singularity.jl")
include("utils/barycentric.jl")
include("utils/bernstein.jl")

# # solvers
include("solver/dielectric_box2d.jl")
# include("solver/dielectric_box3d.jl")

# # visualization
include("visualization/viz_2d.jl")
# include("visualization/viz_3d.jl")

end
