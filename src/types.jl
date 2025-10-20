
struct Surface{T, N}
    n::Int
    points::Vector{NTuple{N, T}}
    norms::Vector{NTuple{N, T}}
    weights::Vector{T}
end

Base.show(io::IO, s::Surface{T, N}) where {T, N} = print(io, "Surface in $N-dimensional space, with $(length(s.points)) points in $T")