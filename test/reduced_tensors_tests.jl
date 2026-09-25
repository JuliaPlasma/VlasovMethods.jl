using VlasovMethods
using GeometricBrackets: Arakawa, PoissonTensor
using Test

@testset "ReducedTensor" begin
    nx, nv = 5, 4
    N = nx * nv
    tensor = PoissonTensor(Float64, nx, nv, Arakawa(nx, nv, 1 / nx, 2 / nv))
    Pi = rand(N, 3)
    Pj = rand(N, 2)
    rt = ReducedTensor(tensor, Pi, Pj)

    @test size(rt) == (3, 2, N)

    # the double projection over every index pair, not only the stencil around k
    dense = [sum(tensor[m, n, k] * Pi[m, i] * Pj[n, j] for m in 1:N, n in 1:N)
             for i in 1:3, j in 1:2, k in 1:N]

    @test maximum(abs, dense) > 0
    @test [rt[i, j, k] for i in 1:3, j in 1:2, k in 1:N] ≈ dense
end
