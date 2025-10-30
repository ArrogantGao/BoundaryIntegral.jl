function l2d_point_gi(interface::Interface{T, 2}, src::NTuple{2, T}) where T
    t = 0.0
    for panel in interface.panels
        for (point, weight) in zip(panel.points, panel.weights)
            t += laplace2d_doublelayer(src, point, panel.normal) * weight
        end
    end
    return t
end

function l2d_sphere_gi(radius::T, ns::Int, src::NTuple{2, T}) where T
    t = 0.0
    dl = 2π * radius / ns
    for theta in 0:2π/ns:(2π-2π/ns)
        cs = (cos(theta), sin(theta))
        rs = (radius * cs[1], radius * cs[2])
        t += laplace2d_doublelayer(src, rs, cs) * dl
    end
    return t
end

# integrate the single layer potential over the interface on a sphere
function l2d_singlelayer_gi(interface::Interface{T, 2}, sigma::Vector{T}, radius::T, ns::Int) where T
    s = zero(T)
    dl = 2π * radius / ns

    for theta in 0:2π/ns:(2π-2π/ns)
        cs = (cos(theta), sin(theta))
        rs = (radius * cs[1], radius * cs[2])

        t = zero(T)

        i_offset = 0
        for panel in interface.panels
            for (point, weight) in zip(panel.points, panel.weights)
                i_offset += 1
                t += laplace2d_doublelayer(point, rs, cs) * weight * sigma[i_offset]
            end
        end

        s += t * dl
    end

    return s
end

function l2d_singlelayer_gi(dbox::DielectricInterfaces{T, 2}, sigma::Vector{T}, radius::T, ns::Int) where T
    s = zero(T)
    dl = 2π * radius / ns

    for theta in 0:2π/ns:(2π-2π/ns)
        cs = (cos(theta), sin(theta))
        rs = (radius * cs[1], radius * cs[2])

        t = zero(T)
        for point in eachpoint(dbox)
            t += laplace2d_doublelayer(point.point, rs, cs) * point.weight * sigma[point.global_idx]
        end
        s += t * dl
    end
    return s
end