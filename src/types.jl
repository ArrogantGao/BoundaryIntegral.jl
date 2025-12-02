# each panel is a line segment with quadrature
struct Panel{T, N}
    n::Int
    points::Vector{NTuple{N, T}}
    normal::NTuple{N, T}
    weights::Vector{T}
end
Base.show(io::IO, p::Panel{T, N}) where {T, N} = print(io, "Panel in $N-dimensional space, with $(length(p.points)) quadrature points in $T")

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

# The default iterator for DielectricInterfaces iterates over the interfaces
Base.iterate(d::DielectricInterfaces{T, N}, state::Int = 1) where {T, N} = iterate(d.interfaces, state)

# To iterate over each point in all panels of all interfaces, use `eachpoint` helper
struct PointInfo{T, N}
    point::NTuple{N, T}
    normal::NTuple{N, T}
    weight::T
    global_idx::Int
    interface_idx::Int
    eps_in::T
    eps_out::T
end

struct AllPointsIterator{T, N}
    d::DielectricInterfaces{T, N}
end

eachpoint(d::DielectricInterfaces) = AllPointsIterator(d)

function Base.length(it::AllPointsIterator)
    return num_points(it.d)
end

function Base.eltype(::AllPointsIterator{T, N}) where {T, N}
    return PointInfo{T, N}
end

function Base.iterate(it::AllPointsIterator{T, N}, state = (1, 1, 1, 1)) where {T, N}
    interface_idx, panel_idx, point_idx, global_idx = state

    if interface_idx > length(it.d.interfaces)
        return nothing
    end

    interface, eps_in, eps_out = it.d.interfaces[interface_idx]

    if panel_idx > length(interface.panels)
        return Base.iterate(it, (interface_idx + 1, 1, 1, global_idx))
    end

    panel = interface.panels[panel_idx]

    if point_idx > length(panel.points)
        return Base.iterate(it, (interface_idx, panel_idx + 1, 1, global_idx))
    end

    point = panel.points[point_idx]
    weight = panel.weights[point_idx]
    normal = panel.normal

    info = PointInfo(point, normal, weight, global_idx, interface_idx, eps_in, eps_out)

    return (info, (interface_idx, panel_idx, point_idx + 1, global_idx + 1))
end