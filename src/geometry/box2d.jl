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

# constructing boxes
function uniform_box2d(n_quad::Int, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n_quad)
    t1 = one(T)
    t0 = zero(T)

    panels = Vector{Panel{T, 2}}()

    for (sp, ep, normal) in zip([(-t1, t1), (t1, t1), (t1, -t1), (-t1, -t1)], [(t1, t1), (t1, -t1), (-t1, -t1), (-t1, t1)], [(t0, t1), (t1, t0), (t0, -t1), (-t1, t0)])
        push!(panels, straight_line_panel(sp, ep, ns, ws, normal))
    end

    return Interface(4, panels)
end


# generate uniform sized panels on boundary of a 2d box, each panel has n_quad quadrature points
function box2d_uniform_panels(n_panels::Int, n_quad::Int, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n_quad)
    t1 = one(T)
    t0 = zero(T)
    panels = Vector{Panel{T, 2}}()

    for (sp, ep, normal) in zip([(-t1, t1), (t1, t1), (t1, -t1), (-t1, -t1)], [(t1, t1), (t1, -t1), (-t1, -t1), (-t1, t1)], [(t0, t1), (t1, t0), (t0, -t1), (-t1, t0)])
        for i in 1:n_panels
            p_start = sp .+ (ep .- sp) .* (i - 1) ./ n_panels
            p_end = sp .+ (ep .- sp) .* i ./ n_panels
            push!(panels, straight_line_panel(p_start, p_end, ns, ws, normal))
        end
    end

    return Interface(4 * n_panels, panels)
end

function box2d_adaptive_panels(n_panels::Int, n_quad::Int, n_adapt::Int, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n_quad)
    panels = Vector{Panel{T, 2}}()

    t1 = one(T)
    t0 = zero(T)

    for (sp, ep, normal) in zip([(-t1, t1), (t1, t1), (t1, -t1), (-t1, -t1)], [(t1, t1), (t1, -t1), (-t1, -t1), (-t1, t1)], [(t0, t1), (t1, t0), (t0, -t1), (-t1, t0)])

        new_panels = straight_line_adaptive_panels(sp, ep, ns, ws, normal, n_panels, n_adapt)
        append!(panels, new_panels)
    end

    return Interface(length(panels), panels)
end

# first try a special case: 2 boxes
function dbox2d_adaptive_panels(n_panels::Int, n_quad::Int, n_adapt::Int, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n_quad)
    panels_1 = Vector{Panel{T, 2}}()
    panels_2 = Vector{Panel{T, 2}}()
    panels_12 = Vector{Panel{T, 2}}()

    t1 = one(T)
    t0 = zero(T)

    lines_1 = [(t0, -t1), (-t1, -t1), (-t1, t1), (t0, t1)]
    norms_1 = [(t0, -t1), (-t1, t0), (t0, t1)]
    lines_2 = [(t0, t1), (t1, t1), (t1, -t1), (t0, -t1)]
    norms_2 = [(t0, t1), (t1, t0), (t0, -t1)]

    for i in 1:3
        append!(panels_1, straight_line_adaptive_panels(lines_1[i], lines_1[i+1], ns, ws, norms_1[i], n_panels, n_adapt))
        append!(panels_2, straight_line_adaptive_panels(lines_2[i], lines_2[i+1], ns, ws, norms_2[i], n_panels, n_adapt))
    end

    append!(panels_12, straight_line_adaptive_panels((t0, t1), (t0, -t1), ns, ws, (-t1, t0), n_panels, n_adapt))

    return [Interface(length(panels_1), panels_1), Interface(length(panels_2), panels_2), Interface(length(panels_12), panels_12)]
end

square = (x, y) -> [(x, y), (x + 1, y), (x + 1, y + 1), (x, y+1)]
rect = (x, y, w, h) -> [(x, y), (x + w, y), (x + w, y + h), (x, y + h)]

# vertices in the same box are arranged counter-clockwise, for example: [(0,0), (1,0), (1,1), (0,1)]
# id of the box is the index of the box in the vector, id of vacuum is 0
# norm vector points from box with higher id to box with lower id
function multi_box2d(n_panels::Int, n_quad::Int, n_adapt::Int, vec_boxes::Vector{Vector{NTuple{2, T}}}, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n_quad)

    ns = T.(ns)
    ws = T.(ws)

    edges = build_edges(vec_boxes)
    interfaces_dict = Dict{Tuple{Int, Int}, Interface{T, 2}}()
    for (p1, p2, id1, id2, nvec) in edges
        if haskey(interfaces_dict, (id1, id2))
            append!(interfaces_dict[(id1, id2)].panels, straight_line_adaptive_panels(p1, p2, ns, ws, nvec, n_panels, n_adapt))
        else
            interfaces_dict[(id1, id2)] = Interface(length(straight_line_adaptive_panels(p1, p2, ns, ws, nvec, n_panels, n_adapt)), straight_line_adaptive_panels(p1, p2, ns, ws, nvec, n_panels, n_adapt))
        end
    end

    return interfaces_dict
end

function edge_normal(p1, p2)
    v = [p2[1]-p1[1], p2[2]-p1[2]]
    n = [v[2], -v[1]]
    n ./= norm(n)
    return tuple(n...)
end

function normalize_edge(p1, p2)
    return (p1 < p2) ? (p1, p2) : (p2, p1)
end

function overlap_segment(p1, p2, p3, p4; tol=1e-10)
    v1 = [p2[1]-p1[1], p2[2]-p1[2]]
    v2 = [p4[1]-p3[1], p4[2]-p3[2]]
    cross1 = abs(det([v1 [p3[1]-p1[1]; p3[2]-p1[2]]]))
    cross2 = abs(det([v1 [p4[1]-p1[1]; p4[2]-p1[2]]]))
    if cross1 > tol || cross2 > tol
        return false, nothing
    end

    if abs(v1[1]) >= abs(v1[2])
        pts = sort!([p1, p2, p3, p4], by = x -> x[1])
        lo, hi = pts[2], pts[3]
    else
        pts = sort!([p1, p2, p3, p4], by = x -> x[2])
        lo, hi = pts[2], pts[3]
    end

    if (hi[1]-lo[1])^2 + (hi[2]-lo[2])^2 > tol^2
        return true, (lo, hi)
    else
        return false, nothing
    end
end

function intersect_edges(rect1::Vector{NTuple{2, T}}, rect2::Vector{NTuple{2, T}}) where T
    for (v11, v12) in collect(zip(rect1, circshift(rect1, -1)))
        for (v21, v22) in collect(zip(rect2, circshift(rect2, -1)))
            # check if the edge (v11, v12) intersects with the edge (v21, v22)
            overlap, lohi = overlap_segment(v11, v12, v21, v22)
            if overlap
                # determine the normal vector of the edge, the same as outward normal vector of the rect2
                nvec = edge_normal(v21, v22)
                return true, (lohi[1], lohi[2], nvec)
            end
        end
    end
    return false, nothing
end

function split_edge_by_overlaps(p1, p2, overlaps; tol=1e-10)
    if abs(p2[1]-p1[1]) >= abs(p2[2]-p1[2])
        # x dominates
        smin, smax = min(p1[1],p2[1]), max(p1[1],p2[1])
        proj = x->x[1]
    else
        # y dominates
        smin, smax = min(p1[2],p2[2]), max(p1[2],p2[2])
        proj = x->x[2]
    end

    # collect the overlapping intervals
    intervals = []
    for (q1,q2) in overlaps
        lo, hi = sort([proj(q1), proj(q2)])
        push!(intervals, (lo,hi))
    end

    # merge the overlapping intervals
    sort!(intervals, by=i->i[1])
    merged = []
    for iv in intervals
        if isempty(merged) || iv[1] > merged[end][2] + tol
            push!(merged, iv)
        else
            merged[end] = (merged[end][1], max(merged[end][2], iv[2]))
        end
    end

    # find the non-overlapping parts
    segments = []
    cursor = smin
    for (lo,hi) in merged
        if lo - cursor > tol
            push!(segments,(cursor,lo))
        end
        cursor = hi
    end
    if smax - cursor > tol
        push!(segments,(cursor,smax))
    end

    # convert back to 2d endpoints
    res = []
    for (a,b) in segments
        if abs(p2[1]-p1[1]) >= abs(p2[2]-p1[2])
            # x dominates -> interpolate y
            yA = p1[2] + (a-p1[1])/(p2[1]-p1[1])*(p2[2]-p1[2])
            yB = p1[2] + (b-p1[1])/(p2[1]-p1[1])*(p2[2]-p1[2])
            push!(res, ((a,yA),(b,yB)))
        else
            # y dominates -> interpolate x
            xA = p1[1] + (a-p1[2])/(p2[2]-p1[2])*(p2[1]-p1[1])
            xB = p1[1] + (b-p1[2])/(p2[2]-p1[2])*(p2[1]-p1[1])
            push!(res, ((xA,a),(xB,b)))
        end
    end
    return res
end


function build_edges(rects::Vector{Vector{NTuple{2,T}}}) where T
    edges = []
    shared_edges = Dict{Tuple{NTuple{2,T},NTuple{2,T}}, Vector{Int}}()

    # find all shared edges and record them
    for i in 1:length(rects)-1
        for j in i+1:length(rects)
            overlap, lohinvec = intersect_edges(rects[i], rects[j])
            if overlap
                p1,p2 = lohinvec[1], lohinvec[2]
                nvec = lohinvec[3]
                key = normalize_edge(p1,p2)
                shared_edges[key] = get(shared_edges,key,Int[]) ∪ [i,j]
                push!(edges, (p1,p2,j,i,nvec))
            end
        end
    end

    # check the rest edges of each rectangle
    for (rid, rect) in enumerate(rects)
        for (p1,p2) in zip(rect, circshift(rect,1))
            key = normalize_edge(p1,p2)
            nvec = edge_normal(p1,p2)

            # collect the overlapping segments
            overlaps = []
            for (k,regs) in shared_edges
                if rid in regs
                    # match the overlapping of the current edge
                    overlap, seg = overlap_segment(p1,p2,k[1],k[2])
                    if overlap
                        push!(overlaps, seg)
                    end
                end
            end

            # split and record the non-overlapping parts
            remain = split_edge_by_overlaps(p1,p2,overlaps)
            for (q1,q2) in remain
                push!(edges,(q1,q2,rid,0,(-nvec[1], -nvec[2])))
            end
        end
    end
    return edges
end