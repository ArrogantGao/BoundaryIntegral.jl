# mesh a triangle surface panel with tensor produce Gauss-Legendre quadrature points
function square_surface_uniform_panel(a::NTuple{3, T}, b::NTuple{3, T}, c::NTuple{3, T}, d::NTuple{3, T}, ns::Vector{T}, ws::Vector{T}, normal::NTuple{3, T}) where T

    # check edge lengths
    Lab = norm(b .- a)
    Lbc = norm(c .- b)
    Lcd = norm(d .- c)
    Lda = norm(a .- d)
    @assert (Lab ≈ Lcd) && (Lbc ≈ Lda) "Edges of the square are not equal"

    # check perpendicularity
    @assert (abs(dot(normal, b .- a)) < 1e-10) && (abs(dot(normal, c .- b)) < 1e-10) && (abs(dot(normal, d .- c)) < 1e-10) && (abs(dot(normal, a .- d)) < 1e-10) "Normal is not perpendicular to the edges"
    @assert (abs(dot(b .- a, c .- b)) < 1e-10) && (abs(dot(c .- b, d .- c)) < 1e-10) && (abs(dot(d .- c, a .- d)) < 1e-10) && (abs(dot(a .- d, b .- a)) < 1e-10) "Edges are not perpendicular"


    e_ab = (b .- a) ./ Lab
    e_bc = (c .- b) ./ Lbc

    points = Vector{NTuple{3, T}}()
    for i in 1:length(ns)
        for j in 1:length(ns)
            p = a .+ e_ab .* (ns[i] .+ 1) ./ 2 .* Lab .+ e_bc .* (ns[j] .+ 1) ./ 2 .* Lbc
            push!(points, p)
        end
    end    

    weights = Vector{T}()
    for i in 1:length(ns)
        for j in 1:length(ns)
            push!(weights, ws[i] * ws[j] * Lab * Lbc / 4)
        end
    end
    
    return Panel(length(ns) * length(ns), points, normal, weights)
end

function square_surface_adaptive_panels(a::NTuple{3, T}, b::NTuple{3, T}, c::NTuple{3, T}, d::NTuple{3, T}, ns::Vector{T}, ws::Vector{T}, normal::NTuple{3, T}, is_edge::NTuple{4, Bool}, is_corner::NTuple{4, Bool}, n_adapt_edge::Int, n_adapt_corner::Int) where T

    # no edges or corners
    all(isfalse.(is_edge)) && all(isfalse.(is_corner)) && return [square_surface_uniform_panel(a, b, c, d, ns, ws, normal)]

    # check the correctness of the is_edge and is_corner, a corner always leads to two edges, inv does not need to hold
    for i in 1:4
        is_corner[i] && (@assert is_edge[i] && is_edge[mod1(i + 1, 4)] "A corner must lead to two edges")
    end

    @assert n_adapt_edge <= n_adapt_corner "n_adapt_edge must be less than or equal to n_adapt_corner"

    squares = Vector{NTuple{4, NTuple{3, T}}}()

    # first handle the edges
    _square_surface_edge!(squares, (a, b, c, d), is_edge, n_adapt_edge)

    # then for the corners
    # do not handle corners differently now

    panels = Vector{Panel{T, 3}}()
    for square in squares
        push!(panels, square_surface_uniform_panel(square..., ns, ws, normal))
    end

    return panels
end

# after one step of adaptation, there are at most two edges with singularity
function _square_surface!(squares::Vector{NTuple{4, NTuple{3, T}}}, abcd::NTuple{4, NTuple{3, T}}, is_edge::NTuple{4, Bool}, is_corner::NTuple{4, Bool}, edge_depth::Int, corner_depth::Int) where T

    if (!any(is_edge) && !any(is_corner)) || (edge_depth == 0 && corner_depth == 0)
        push!(squares, abcd)
        return squares
    end

    a, b, c, d = abcd
    h_ab = (a .+ b) ./ 2
    h_bc = (b .+ c) ./ 2
    h_cd = (c .+ d) ./ 2
    h_da = (d .+ a) ./ 2
    c_abcd = (a .+ b .+ c .+ d) ./ 4

    lb = (a, h_ab, c_abcd, h_da)
    rb = (h_ab, b, h_bc, c_abcd)
    rt = (c_abcd, h_bc, c, h_cd)
    lt = (h_da, c_abcd, h_cd, d)

    sub_squares = (lb, rb, rt, lt)

    if (edge_depth == 0) 
        if (corner_depth == 0 || !any(is_corner))
            push!(squares, abcd)
        else
            for i in 1:4
                if !is_corner[i] 
                    push!(squares, sub_squares[i])
                else
                    dummy_is_edge = (false, false, false, false) # no further edge refinement needed
                    new_is_corner_mut = [false, false, false, false]
                    new_is_corner_mut[i] = true
                    _square_surface!(squares, sub_squares[i], dummy_is_edge, Tuple(new_is_corner_mut), 0, corner_depth - 1)
                end
            end
        end
    else
        for i in 1:4
            new_is_corner_mut = [false, false, false, false]
            new_is_corner_mut[i] = is_corner[i]
            new_is_corner = Tuple(new_is_corner_mut)

            j = mod1(i - 1, 4)
            new_is_edge_mut = [false, false, false, false]
            new_is_edge_mut[i] = is_edge[i]
            new_is_edge_mut[j] = is_edge[j]
            new_is_edge = Tuple(new_is_edge_mut)

            _square_surface!(squares, sub_squares[i], new_is_edge, new_is_corner, edge_depth - 1, corner_depth - 1)
        end
    end

    return squares
end

# first try to mesh a single box
