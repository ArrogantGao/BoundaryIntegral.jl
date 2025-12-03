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

function l3d_singlelayer_gi(dbox::DielectricInterfaces{T, 3}, sigma::Vector{T}, radius::T, order::Int) where T
    order_new = order
    while true
        if isavailable(order_new)
            order_new != order && @warn "lebedev order $order_new is the nearest available order"
            break
        else
            order_new += 1
        end
    end

    x, y, z, w = lebedev_by_order(order_new)
    x .*= radius
    y .*= radius
    z .*= radius
    w .*= 4π * radius^2

    t = zero(T)
    for (xi, yi, zi, wi) in zip(x, y, z, w)
        s = zero(T)
        for point in eachpoint(dbox)
            s += laplace3d_doublelayer(point.point, (xi, yi, zi), (xi, yi, zi) ./ radius) * point.weight * sigma[point.global_idx]
        end
        t += s * wi
    end
    return t
end

function l3d_singlelayer_charge(src::NTuple{3, T}, radius::T, order::Int) where T

    while true
        if isavailable(order)
            @warn "lebedev order $order is the nearest available order"
            break
        else
            order += 1
        end
    end

    x, y, z, w = lebedev_by_order(order)
    x .*= radius
    y .*= radius
    z .*= radius
    w .*= 4π * radius^2

    t = zero(T)
    for (xi, yi, zi, wi) in zip(x, y, z, w)
        t += laplace3d_doublelayer(src, (xi, yi, zi), (xi, yi, zi) ./ radius) * wi
    end
    return t
end

function l3d_box_gi(dbox::DielectricInterfaces{T, 3}, src::NTuple{3, T}) where T
    t = zero(T)
    for point in eachpoint(dbox)
        t += laplace3d_doublelayer(src, point.point, point.normal) * point.weight
    end
    return t
end