# compute the radius of the Bernstein ellipse from the pole
function bernstein_rho_from_pole(z::Complex{T}) where T <: Real
    s  = sqrt(z*z - 1)          # complex sqrt (principal branch)
    w1 = z + s
    w2 = z - s
    return max(abs(w1), abs(w2))
end

# estimate the minimum rho for a 2d integral
# swap x and y from -1 to 1 to get the minimum one
# input x, y, z are the relative position of the pole in 3d space
function bernstein_rho_2d(x::T, y::T, z::T) where T <: Real
    rho_x_min = Inf
    rho_y_min = Inf
    # this is a very rough swap, but seems to work
    for t in - 1.0:0.1:1.0
        rho_x = bernstein_rho_from_pole(Complex{T}(x, sqrt((y - t)^2 + z^2)))
        rho_y = bernstein_rho_from_pole(Complex{T}(y, sqrt((x - t)^2 + z^2)))
        rho_x_min = min(rho_x_min, rho_x)
        rho_y_min = min(rho_y_min, rho_y)
    end
    return min(rho_x_min, rho_y_min)
end