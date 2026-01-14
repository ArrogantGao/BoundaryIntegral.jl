function Lhs_dielectric_box2d(interface::DielectricInterface{P, T}) where {P <: AbstractPanel, T}
    D_transpose = laplace2d_DT(interface)
    Lhs = D_transpose
    offset = 0
    for i in 1:length(interface.panels)
        panel = interface.panels[i]
        eps_in = interface.eps_in[i]
        eps_out = interface.eps_out[i]
        n_pts = num_points(panel)
        t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
        for j in 1:n_pts
            Lhs[offset + j, offset + j] += t
        end
        offset += n_pts
    end
    return Lhs
end

function Lhs_dielectric_box2d_fmm2d(interface::DielectricInterface{P, T}, tol::Float64 = 1e-12) where {P <: AbstractPanel, T}
    D_transpose = laplace2d_DT_fmm2d(interface, tol)

    function g(x)
        Dx = D_transpose * x

        offset = 0
        for i in 1:length(interface.panels)
            panel = interface.panels[i]
            eps_in = interface.eps_in[i]
            eps_out = interface.eps_out[i]
            n_pts = num_points(panel)
            t = 0.5 * (eps_out + eps_in) / (eps_out - eps_in)
            for j in 1:n_pts
                Dx[offset + j] += t * x[offset + j]
            end
            offset += n_pts
        end

        return Dx
    end

    Lhs = LinearMap{T}(g, num_points(interface), num_points(interface))

    return Lhs
end

function Rhs_dielectric_box2d(interface::DielectricInterface{P, T}, ps::PointSource{T}, eps_src::T) where {P <: AbstractPanel, T}
    src = ps.point
    q = ps.charge
    n_points = num_points(interface)
    Rhs = zeros(T, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        Rhs[i] = - q * laplace2d_grad(src, point.panel_point.point, point.panel_point.normal) / eps_src
    end
    return Rhs
end