# each panel is a surface with quadrature
struct Panel{T, N}
    n::Int
    points::Vector{NTuple{N, T}}
    normal::NTuple{N, T}
    weights::Vector{T}
end
Base.show(io::IO, p::Panel{T, N}) where {T, N} = print(io, "Panel in $N-dimensional space, with $(length(p.points)) quadrature points in $T")

# 1d line panel with Guass-Legendre quadrature
function straight_line_panel(a::NTuple{2, T}, b::NTuple{2, T}, ns::Vector{T}, ws::Vector{T}, normal::NTuple{2, T}) where T

    points = [(b .+ a) ./ 2 .+ ns[i] .* (b .- a) ./ 2 for i in 1:length(ns)]
    L = norm(b .- a)
    weights = ws .* L ./ 2

    @assert norm(normal) ≈ 1 "Normal is not a unit vector"
    @assert dot(normal, b .- a) ≈ 0 "Normal is not perpendicular to the line segment"

    return Panel(length(ns), points, normal, weights)
end

struct Surface{T, N}
    n::Int
    panels::Vector{Panel{T, N}}
end

Base.show(io::IO, s::Surface{T, N}) where {T, N} = print(io, "Surface in $N-dimensional space, with $(length(s.panels)) panels in $T")