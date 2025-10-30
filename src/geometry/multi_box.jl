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

    append!(panels_12, straight_line_adaptive_panels((t0, t1), (t0, -t1), ns, ws, (t1, t0), n_panels, n_adapt))

    return [Interface(length(panels_1), panels_1), Interface(length(panels_2), panels_2), Interface(length(panels_12), panels_12)]
end

# function multi_box2d(n_boxes::Int, vertices::Vector{NTuple{4, NTuple{2, T}}}) where T
#     @assert length(vertices) == n_boxes

#     # each length 4 tuple describes a box

# end