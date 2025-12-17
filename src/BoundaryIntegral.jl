module BoundaryIntegral

using LinearAlgebra, OMEinsum
using FastGaussQuadrature, Lebedev
using Krylov, LinearMaps, Roots

using FMM2D, FMM3D

using KernelAbstractions, Adapt
include(joinpath(dirname(pathof(KernelAbstractions)), "../examples/utils.jl")) # Load backend

#core types
export Interface

# kernel functions
export laplace3d_doublelayer, laplace2d_doublelayer, laplace3d_singlelayer, laplace2d_singlelayer

# geometries
export uniform_box3d, uniform_box2d

# linear algebra
export solve_lu, solve_gmres

# visualization
export viz_2d_interfaces, viz_2d_dielectric_interfaces
export viz_3d_squares, viz_3d_interface, viz_3d_dielectric_interfaces

include("types.jl")

# kernel functions
include("kernel/laplace2d.jl")
include("kernel/laplace3d.jl")
include("kernel/kernelabstractions.jl")

# geometries
include("geometry/box2d.jl")
include("geometry/box3d.jl")

# utilities
include("utils/greens_identity.jl")
include("utils/linear_algebra.jl")
include("utils/corner_singularity.jl")
include("utils/legendre.jl")

# solvers
include("solver/dielectric_box2d.jl")
include("solver/dielectric_box3d.jl")

# visualization
include("visualization/viz_2d.jl")
include("visualization/viz_3d.jl")

end
