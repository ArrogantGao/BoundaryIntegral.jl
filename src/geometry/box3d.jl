# mesh a rectangle surface panel with tensor product Gauss-Legendre quadrature points
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

    cc = (a .+ b .+ c .+ d) ./ 4

    points = Vector{NTuple{3, T}}()
    for i in 1:length(ns)
        for j in 1:length(ns)
            # p = a .+ e_ab .* (ns[i] .+ 1) ./ 2 .* Lab .+ e_bc .* (ns[j] .+ 1) ./ 2 .* Lbc
            p = cc .+ (b .- a) .* (ns[i] / 2) .+ (d .- a) .* (ns[j] / 2)
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

    @assert n_adapt_edge <= n_adapt_corner "n_adapt_edge must be less than or equal to n_adapt_corner"

    squares = Vector{NTuple{4, NTuple{3, T}}}()

    # first handle the edges
    _square_surfaces!(squares, (a, b, c, d), is_edge, is_corner, n_adapt_edge, n_adapt_corner)

    panels = Vector{Panel{T, 3}}()
    for square in squares
        push!(panels, square_surface_uniform_panel(square..., ns, ws, normal))
    end

    return panels
end

function square_surface_adaptive_panels(a::NTuple{3, T}, b::NTuple{3, T}, c::NTuple{3, T}, d::NTuple{3, T}, n_quad_max::Int, n_quad_min::Int, normal::NTuple{3, T}, is_edge::NTuple{4, Bool}, is_corner::NTuple{4, Bool}, n_adapt_edge::Int, n_adapt_corner::Int) where T

    @assert n_adapt_edge <= n_adapt_corner "n_adapt_edge must be less than or equal to n_adapt_corner"


    if n_adapt_edge == 0
        ns, ws = gausslegendre(n_quad_min)
        return [square_surface_uniform_panel(a, b, c, d, ns, ws, normal)]
    end

    squares = Vector{NTuple{4, NTuple{3, T}}}()

    # first handle the edges
    _square_surfaces!(squares, (a, b, c, d), is_edge, is_corner, n_adapt_edge, n_adapt_corner)

    L_ab = norm(b .- a)
    panels = Vector{Panel{T, 3}}()

    L_max = zero(T)
    L_min = L_ab
    for square in squares
        L_max = max(L_max, norm(square[2] .- square[1]))
        L_min = min(L_min, norm(square[2] .- square[1]))
    end

    for square in squares
        L = norm(square[2] .- square[1])
        r = (L_max == L_min) ? 1 : (L - L_min) / (L_max - L_min)
        n_quad = ceil(Int, n_quad_min + (n_quad_max - n_quad_min) * r)
        ns, ws = gausslegendre(n_quad)
        push!(panels, square_surface_uniform_panel(square..., ns, ws, normal))
    end

    return panels
end

function _square_surfaces!(squares::Vector{NTuple{4, NTuple{3, T}}}, abcd::NTuple{4, NTuple{3, T}}, is_edge::NTuple{4, Bool}, is_corner::NTuple{4, Bool}, edge_depth::Int, corner_depth::Int) where T

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
                    _square_surfaces!(squares, sub_squares[i], dummy_is_edge, Tuple(new_is_corner_mut), 0, corner_depth - 1)
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

            _square_surfaces!(squares, sub_squares[i], new_is_edge, new_is_corner, edge_depth - 1, corner_depth - 1)
        end
    end

    return squares
end

# first try to mesh a single box
# each surfaces are first divided into n_boxexs * n_boxes sub_surfaces
function single_box3d(n_boxes::Int, n_quad_max::Int, n_quad_min::Int, n_adapt_edge::Int, n_adapt_corner::Int, ::Type{T} = Float64) where T
    t1 = one(T)
    t0 = zero(T)

    vertices = [
        ( t1,  t1,  t1),  # 1: A
        (-t1,  t1,  t1),  # 2: B
        (-t1, -t1,  t1),  # 3: C
        ( t1, -t1,  t1),  # 4: D
        ( t1,  t1, -t1),  # 5: E
        (-t1,  t1, -t1),  # 6: F
        (-t1, -t1, -t1),  # 7: G
        ( t1, -t1, -t1),  # 8: H
    ]

    faces = [
        (1, 2, 3, 4),  # z = +1  (A B C D)
        (5, 8, 7, 6),  # z = -1  (E H G F)
        (8, 5, 1, 4),  # x = +1  (H E A D)
        (7, 3, 2, 6),  # x = -1  (G C B F)
        (6, 2, 1, 5),  # y = +1  (F B A E)
        (7, 8, 4, 3),  # y = -1  (G H D C)
    ]

    normals = [
        ( t0,  t0,  t1),  # (1,2,3,4)
        ( t0,  t0, -t1),  # (5,8,7,6)
        ( t1,  t0,  t0),  # (8,5,1,4)
        (-t1,  t0,  t0),  # (7,3,2,6)
        ( t0,  t1,  t0),  # (6,2,1,5)
        ( t0, -t1,  t0),  # (7,8,4,3)
    ]

    panels = Vector{Panel{T, 3}}()
    for i in 1:6
        face = faces[i]
        a, b, c, d = [vertices[j] for j in face]
        r_ab = b .- a
        r_ad = d .- a
        d_ab = r_ab ./ n_boxes
        d_ad = r_ad ./ n_boxes
        normal = normals[i]
        for face_idx in 1:n_boxes
            for face_idy in 1:n_boxes
                af = a .+ d_ab .* (face_idx - 1) .+ d_ad .* (face_idy - 1)
                bf = af .+ d_ab
                cf = af .+ d_ab .+ d_ad
                df = af .+ d_ad

                is_edge_mut = [false, false, false, false]
                is_corner_mut = [false, false, false, false]

                is_edge_mut[4] = (face_idx == 1)
                is_edge_mut[2] = (face_idx == n_boxes)
                is_edge_mut[1] = (face_idy == 1)
                is_edge_mut[3] = (face_idy == n_boxes)

                is_corner_mut[1] = is_edge_mut[1] && is_edge_mut[4]
                is_corner_mut[2] = is_edge_mut[1] && is_edge_mut[2]
                is_corner_mut[3] = is_edge_mut[2] && is_edge_mut[3]
                is_corner_mut[4] = is_edge_mut[3] && is_edge_mut[4]

                # new_panels = square_surface_adaptive_panels(af, bf, cf, df, ns, ws, normal, Tuple(is_edge_mut), Tuple(is_corner_mut), n_adapt_edge, n_adapt_corner)
                new_panels = square_surface_adaptive_panels(af, bf, cf, df, n_quad_max, n_quad_min, normal, Tuple(is_edge_mut), Tuple(is_corner_mut), n_adapt_edge, n_adapt_corner)
                append!(panels, new_panels)
            end
        end
    end

    return Interface(length(panels), panels)
end

# generate a cubic box with left front bottom corner at lfd and right upper back corner at rbu
function cubic_box3d(lfd::NTuple{3, T}, rbu::NTuple{3, T}) where T
    xL, yF, zD = lfd  # left, front, down
    xR, yB, zU = rbu  # right, back, up

    @assert xR > xL && yB > yF && zU > zD "rbu must be the opposite corner of lfd (xR>xL, yB>yF, zU>zD)"

    # 1: (right, front, up)
    # 2: (left,  front, up)
    # 3: (left,  back,  up)
    # 4: (right, back,  up)
    # 5: (right, front, down)
    # 6: (left,  front, down)
    # 7: (left,  back,  down)
    # 8: (right, back,  down)
    vertices = NTuple{3,T}[
        (xR, yF, zU),  # 1
        (xL, yF, zU),  # 2
        (xL, yB, zU),  # 3
        (xR, yB, zU),  # 4
        (xR, yF, zD),  # 5
        (xL, yF, zD),  # 6
        (xL, yB, zD),  # 7
        (xR, yB, zD),  # 8
    ]

    faces = NTuple{4,Int}[
        (1, 2, 3, 4),  # z = zU  (top),    n = (0, 0, +1)
        (5, 8, 7, 6),  # z = zD  (bottom), n = (0, 0, -1)
        (8, 5, 1, 4),  # x = xR  (right),  n = (+1, 0, 0)
        (7, 3, 2, 6),  # x = xL  (left),   n = (-1, 0, 0)
        (6, 2, 1, 5),  # y = yF  (front),  n = (0, +1, 0)
        (7, 8, 4, 3),  # y = yB  (back),   n = (0, -1, 0)
    ]

    z0 = zero(T)
    o  = one(T)

    normals = NTuple{3,T}[
        ( z0,  z0,  o ),  # top    (+z)
        ( z0,  z0, -o ),  # bottom (-z)
        (  o,  z0,  z0),  # right  (+x)
        ( -o,  z0,  z0),  # left   (-x)
        ( z0,  o,  z0),   # front  (+y)
        ( z0, -o,  z0),   # back   (-y)
    ]

    return vertices, faces, normals
end