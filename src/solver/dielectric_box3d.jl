function Lhs_dielectric_mbox3d_direct(dbox::DielectricInterfaces{T, 3}) where T
    D_transpose = laplace3d_DT(dbox)
    Lhs = D_transpose
    offset = 0
    for (interface, eps_in, eps_out) in dbox.interfaces
        t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
        for i in 1:num_points(interface)
            Lhs[offset + i, offset + i] += t
        end
        offset += num_points(interface)
    end
    return Lhs
end

function Lhs_dielectric_mbox3d_fmm3d(dbox::DielectricInterfaces{T, 3}, tol::Float64 = 1e-6) where T
    D_transpose = laplace3d_DT_fmm3d(dbox, tol)

    function g(x)
        Dx = D_transpose * x

        offset = 0
        for (interface, eps_in, eps_out) in dbox.interfaces
            t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
            for i in 1:num_points(interface)
                Dx[offset + i] += t * x[offset + i]
            end
            offset += num_points(interface)
        end

        return Dx
    end

    Lhs = LinearMap{T}(g, num_points(dbox), num_points(dbox))

    return Lhs
end

function Lhs_dielectric_mbox3d_fmm3d_weighted(dbox::DielectricInterfaces{T, 3}, tol::Float64 = 1e-6) where T
    D_transpose = laplace3d_DT_fmm3d(dbox, tol)
    sqrt_weights = sqrt.(all_weights(dbox))

    function g(x)
        Dx = D_transpose * x

        offset = 0
        for (interface, eps_in, eps_out) in dbox.interfaces
            t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
            for i in 1:num_points(interface)
                Dx[offset + i] += t * x[offset + i]
            end
            offset += num_points(interface)
        end
        
        Dx .*= sqrt_weights

        return Dx
    end

    Lhs = LinearMap{T}(g, num_points(dbox), num_points(dbox))

    return Lhs
end

function Rhs_dielectric_mbox3d(dbox::DielectricInterfaces{T, 3}, src::NTuple{3, T}, eps_src::T) where T
    n_points = num_points(dbox)
    Rhs = zeros(T, n_points)
    for point in eachpoint(dbox)
        Rhs[point.global_idx] = - laplace3d_doublelayer(src, point.point, point.normal) / eps_src
    end
    return Rhs
end

function Rhs_dielectric_mbox3d_weighted(dbox::DielectricInterfaces{T, 3}, src::NTuple{3, T}, eps_src::T) where T
    n_points = num_points(dbox)
    sqrt_weights = sqrt.(all_weights(dbox))
    Rhs = zeros(T, n_points)
    for point in eachpoint(dbox)
        Rhs[point.global_idx] = - laplace3d_doublelayer(src, point.point, point.normal) / eps_src
    end
    Rhs .*= sqrt_weights
    return Rhs
end

function nystrom_interpolation_dielectric_box3d(dbox::DielectricInterfaces{T, 3}, target_interface::DielectricInterfaces{T, 3}, src::NTuple{3, T}, eps_src::T, sigma::Vector{T}, tol::Float64 = 1e-6) where T
    f = Rhs_dielectric_mbox3d(target_interface, src, eps_src)
    K = laplace3d_DT_trg_fmm3d(dbox, target_interface, tol)

    sigma_f = f - K * sigma

    offset = 0
    for (interface, eps_in, eps_out) in target_interface.interfaces
        t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
        (@view sigma_f[offset + 1:offset + num_points(interface)]) ./= t
        offset += num_points(interface)
    end

    return sigma_f
end