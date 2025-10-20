module BoundaryIntegral

using LinearAlgebra
using FastGaussQuadrature

#core types
export Surface

# kernel functions
export coulomb_En

# geometries
export uniform_box

include("types.jl")

# kernel functions
include("kernel/coulomb.jl")

# geometries
include("geometry/box.jl")

end
