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
    @assert dot(normal, b .- a) < 1e-10 "Normal is not perpendicular to the line segment"

    return Panel(length(ns), points, normal, weights)
end

function straight_line_adaptive_panels(sp::NTuple{2, T}, ep::NTuple{2, T}, ns::Vector{T}, ws::Vector{T}, normal::NTuple{2, T}, n_panels::Int, n_adapt::Int) where T
    panels = Vector{Panel{T, 2}}()
    dt = (ep .- sp) ./ n_panels
    # adaptively refine the [-t1, -t1 + dt]
    t_start = sp .+ dt
    for i in 1:n_adapt-1
        t_end = t_start
        t_start = t_end .- dt ./ T(2^i)
        push!(panels, straight_line_panel(t_start, t_end, ns, ws, normal))
    end
    push!(panels, straight_line_panel(sp, t_start, ns, ws, normal))

    for i in 2:n_panels-1
        p_start = sp .+ (ep .- sp) .* (i - 1) ./ n_panels
        p_end = sp .+ (ep .- sp) .* i ./ n_panels
        push!(panels, straight_line_panel(p_start, p_end, ns, ws, normal))
    end

    # adaptively refine the [t1 - dt, t1]
    t_end = ep .- dt
    for i in 1:n_adapt-1
        t_start = t_end
        t_end = t_start .+ dt ./ T(2^i)
        push!(panels, straight_line_panel(t_start, t_end, ns, ws, normal))
    end
    push!(panels, straight_line_panel(t_end, ep, ns, ws, normal))

    return panels
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