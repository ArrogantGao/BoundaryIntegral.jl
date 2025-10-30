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