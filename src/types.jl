# each panel is a line segment with quadrature
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

struct Interface{T, N}
    n::Int
    panels::Vector{Panel{T, N}}
end

Base.show(io::IO, s::Interface{T, N}) where {T, N} = print(io, "Interface in $N-dimensional space, with $(length(s.panels)) panels in $T")

function num_points(interface::Interface{T, N}) where {T, N}
    return sum(panel.n for panel in interface.panels)
end

function all_weights(interface::Interface{T, N}) where {T, N}
    weights = Vector{T}()
    for panel in interface.panels
        for weight in panel.weights
            push!(weights, weight)
        end
    end
    return weights
end

struct DielectricInterfaces{T, N}
    n::Int
    interfaces::Vector{Tuple{Interface{T, N}, T, T}} # (interface, eps_in, eps_out)
end

function num_points(d::DielectricInterfaces{T, N}) where {T, N}
    return sum(num_points(interface) for (interface, eps_in, eps_out) in d.interfaces)
end

function all_weights(d::DielectricInterfaces{T, N}) where {T, N}
    weights = Vector{T}()
    for (interface, eps_in, eps_out) in d.interfaces
        append!(weights, all_weights(interface))
    end
    return weights
end