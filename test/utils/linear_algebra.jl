using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Test

@testset "linear_algebra" begin
    A = rand(100, 100)
    b = rand(100)
    x = BI.solve_lu(A, b)
    @test norm(A * x - b) < 1e-10
end

@testset "solve_gmres" begin
    A = rand(100, 100)
    b = rand(100)
    x = BI.solve_gmres(A, b)
    @test norm(A * x - b) < 1e-10
end