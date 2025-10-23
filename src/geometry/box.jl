# function uniform_box3d(n::Int, ::Type{T} = Float64) where T
#     ns, ws = gausslegendre(n)
#     t1 = one(T)
#     t0 = zero(T)

#     locations = Vector{NTuple{3, T}}()
#     norms = Vector{NTuple{3, T}}()
#     weights = Vector{T}()

#     # generate locations, norms and weights on each surfaces
#     for norm in [
#         (t1, t0, t0), (-t1, t0, t0),
#         (t0, t1, t0), (t0, -t1, t0),
#         (t0, t0, t1), (t0, t0, -t1)]
#         for i in 1:n, j in 1:n
#             w = ws[i] * ws[j]
#             n1, n2 = ns[i], ns[j]

#             if norm[1] != t0
#                 loc = (norm[1], n1, n2)
#             elseif norm[2] != t0
#                 loc = (n1, norm[2], n2)
#             else
#                 loc = (n1, n2, norm[3])
#             end

#             push!(locations, loc)
#             push!(norms, norm)
#             push!(weights, w)
#         end
#     end

#     return Surface(6 * n^2, locations, norms, weights)
# end

function uniform_box2d(n_quad::Int, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n_quad)
    t1 = one(T)
    t0 = zero(T)

    panels = Vector{Panel{T, 2}}()

    for (sp, ep, normal) in zip([(-t1, t1), (t1, t1), (t1, -t1), (-t1, -t1)], [(t1, t1), (t1, -t1), (-t1, -t1), (-t1, t1)], [(t0, t1), (t1, t0), (t0, -t1), (-t1, t0)])
        push!(panels, straight_line_panel(sp, ep, ns, ws, normal))
    end

    return Surface(4, panels)
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

    return Surface(4 * n_panels, panels)
end

function box2d_adaptive_panels(n_panels::Int, n_quad::Int, n_adapt::Int, ::Type{T} = Float64) where T
    ns, ws = gausslegendre(n_quad)
    panels = Vector{Panel{T, 2}}()

    t1 = one(T)
    t0 = zero(T)

    for (sp, ep, normal) in zip([(-t1, t1), (t1, t1), (t1, -t1), (-t1, -t1)], [(t1, t1), (t1, -t1), (-t1, -t1), (-t1, t1)], [(t0, t1), (t1, t0), (t0, -t1), (-t1, t0)])

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
    end

    return Surface(length(panels), panels)
end