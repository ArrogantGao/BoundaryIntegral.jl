using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Random
using Test

@testset "laplace2d FMM2D comparisons" begin
    Random.seed!(0)
    interface = BI.single_dielectric_box2d(1.0, 1.0, 8, 0.2, 0.05, 5.0, 1.0, Float64)

    points = NTuple{2, Float64}[]
    normals = NTuple{2, Float64}[]
    weights = Float64[]
    for point in BI.eachpoint(interface)
        push!(points, point.panel_point.point)
        push!(normals, point.panel_point.normal)
        push!(weights, point.panel_point.weight)
    end

    n = length(points)
    charges = randn(n)
    tol = 1e-12

    # S: on-surface potential, skipping self interactions
    targets = zeros(Float64, 2, n)
    for i in 1:n
        targets[1, i] = points[i][1]
        targets[2, i] = points[i][2]
    end
    S_direct = BI.laplace2d_S(interface)
    S_fmm = BI.laplace2d_S_fmm2d(interface, tol)
    @test norm(S_direct * charges - S_fmm * charges) < 1e-10

    # D and DT: compare direct matrices with FMM2D linear maps
    D_direct = BI.laplace2d_D(interface)
    DT_direct = BI.laplace2d_DT(interface)
    D_fmm = BI.laplace2d_D_fmm2d(interface, tol)
    DT_fmm = BI.laplace2d_DT_fmm2d(interface, tol)
    @test norm(D_direct * charges - D_fmm * charges) < 1e-10
    @test norm(DT_direct * charges - DT_fmm * charges) < 1e-10

    # pottrg: off-surface targets
    m = n + 3
    targets_off = zeros(Float64, 2, m)
    for i in 1:m
        idx = ((i - 1) % n) + 1
        targets_off[1, i] = points[idx][1] + 0.1 * normals[idx][1]
        targets_off[2, i] = points[idx][2] + 0.1 * normals[idx][2]
    end
    pot_direct = BI.laplace2d_pottrg(interface, targets_off) * charges
    pot_fmm = BI.laplace2d_pottrg_fmm2d(interface, targets_off, tol) * charges
    @test norm(pot_direct - pot_fmm) < 1e-10
end

@testset "laplace2d D * I" begin
    # uni-box
    interface = BI.single_dielectric_box2d(1.0, 1.0, 8, 0.2, 0.05, 5.0, 1.0, Float64)
    D_direct = BI.laplace2d_D(interface)
    D_fmm = BI.laplace2d_D_fmm2d(interface, 1e-12)
    weights = BI.all_weights(interface)
    @test dot(abs.(D_direct * ones(size(D_direct, 1)) .- 0.5), weights) < 1e-3
    @test dot(abs.(D_fmm * ones(size(D_fmm, 1)) .- 0.5), weights) < 1e-3

    # rectangle box
    interface = BI.single_dielectric_box2d(2.0, 3.0, 8, 0.2, 0.05, 5.0, 1.0, Float64)
    D_direct = BI.laplace2d_D(interface)
    D_fmm = BI.laplace2d_D_fmm2d(interface, 1e-12)
    weights = BI.all_weights(interface)
    @test dot(abs.(D_direct * ones(size(D_direct, 1)) .- 0.5), weights) < 1e-3
    @test dot(abs.(D_fmm * ones(size(D_fmm, 1)) .- 0.5), weights) < 1e-3
end