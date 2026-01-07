using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using FastGaussQuadrature
using Test

@testset "barycentric_1d" begin
    f = x -> x^2 * sin(x) + exp(-x^2)

    # on the original grid
    x, w = gausslegendre(32)
    xt, _ = gausslegendre(128)
    E = BI.interp_matrix_1d_gl(x, w, xt)
    y = E * f.(x)
    @test norm(y - f.(xt)) < 1e-10
end

@testset "barycentric_2d_tensor" begin
    f = (x, y) -> x^2 * y^3 * sin(x) * sin(1.1 * y) + exp(-x^2 - y^2)
    x, w = gausslegendre(32)
    y, w = gausslegendre(32)

    xt, _ = gausslegendre(64)
    yt, _ = gausslegendre(64)

    E = BI.interp_matrix_2d_gl_tensor(x, w, y, w, xt, yt)

    val = zeros(length(x) * length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            val[i + (j-1) * length(x)] = f(x[i], y[j])
        end
    end

    y = E * val

    val_t = zeros(length(xt) * length(yt))
    for i in 1:length(xt)
        for j in 1:length(yt)
            val_t[i + (j-1) * length(xt)] = f(xt[i], yt[j])
        end
    end

    @test norm(y - val_t) < 1e-10
end

@testset "barycentric_2d_scattered" begin
    f = (x, y) -> x^2 * y^3 * sin(x) * sin(1.1 * y) + exp(-x^2 - y^2)
    x, w = gausslegendre(32)
    y, w = gausslegendre(32)

    xt = rand(20)
    yt = rand(20)

    E_scattered = BI.interp_matrix_2d_gl_scattered(x, w, y, w, xt, yt)

    val = zeros(length(x) * length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            val[i + (j-1) * length(x)] = f(x[i], y[j])
        end
    end

    y_scattered = E_scattered * val

    @test norm(y_scattered - f.(xt, yt)) < 1e-10
end