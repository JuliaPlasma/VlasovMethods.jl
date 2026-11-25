using Test
using Random: Xoshiro
using SimpleSplines: SplineQuadrature, PeriodicBSplineBasis, BSplineBasis, UniformMesh,
                     Free,
                     Periodic, Dirichlet, (..)
using VlasovMethods
using VlasovMethods: evaluate

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

@testset "construction ($T)" for T in (Float64, Float32)
    xb, vb = spline_bases(T)
    A = randn(Xoshiro(1), T, 16, 100)
    dist = TensorSplineDistribution(xb, vb, A)
    @test dist isa VlasovMethods.DistributionFunction{T, 1, 1}
    @test dist.xbasis === xb && dist.vbasis === vb
    @test dist.coefficients === A
    @test dist.xquadrature isa SplineQuadrature{T} &&
          dist.vquadrature isa SplineQuadrature{T}
    @test dist.xquadrature.nq == SplineQuadrature(xb).nq
    @test dist.vquadrature.nq == SplineQuadrature(vb).nq

    # the x-basis is periodic
    clamped = BSplineBasis(UniformMesh(13, zero(T) .. T(4π)), 3, Free())
    @test_throws r"^ArgumentError: the x-basis of a TensorSplineDistribution is periodic" TensorSplineDistribution(
        clamped, vb, randn(Xoshiro(1), T, 16, 100))
    @test_throws r"^ArgumentError: the x-basis of a TensorSplineDistribution is periodic" TensorSplineDistribution(
        clamped, vb, (x, v) -> one(x))

    # the v-basis reproduces v²: Dirichlet, periodic and linear clamped bases do not
    for vbad in (BSplineBasis(UniformMesh(97, T(-25) .. T(25)), 3, Dirichlet()),
        BSplineBasis(UniformMesh(100, T(-25) .. T(25)), 3, Periodic()),
        BSplineBasis(UniformMesh(99, T(-25) .. T(25)), 1, Free()))
        err = try
            TensorSplineDistribution(xb, vbad, zeros(T, 16, nbasis(vbad)))
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("v²", err.msg)
        @test startswith(err.msg, "the v-basis of a TensorSplineDistribution reproduces v²")
        @test_throws r"^ArgumentError: the v-basis of a TensorSplineDistribution reproduces v²" TensorSplineDistribution(
            xb, vbad, (x, v) -> one(x))
    end
    # a quadratic clamped basis reproduces v² and is accepted
    vq2 = BSplineBasis(UniformMesh(98, T(-25) .. T(25)), 2, Free())
    @test TensorSplineDistribution(xb, vq2, zeros(T, 16, 100)) isa TensorSplineDistribution

    # one element type
    S = T === Float64 ? Float32 : Float64
    one_type = r"^ArgumentError: the bases and the coefficients of a TensorSplineDistribution have one element type"
    @test_throws one_type TensorSplineDistribution(xb, vb, zeros(S, 16, 100))
    xbS, _ = spline_bases(S)
    @test_throws one_type TensorSplineDistribution(xbS, vb, zeros(T, 16, 100))
    @test_throws one_type TensorSplineDistribution(xbS, vb, (x, v) -> one(x))

    # the size of the coefficients
    @test_throws r"^DimensionMismatch: the coefficients of a TensorSplineDistribution have the size \(16, 100\) of the bases; got \(16, 99\)$" TensorSplineDistribution(
        xb, vb, zeros(T, 16, 99))
    @test_throws r"^DimensionMismatch: the coefficients of a TensorSplineDistribution have the size \(16, 100\) of the bases; got \(15, 100\)$" TensorSplineDistribution(
        xb, vb, zeros(T, 15, 100))
end

@testset "projection and evaluation ($T)" for T in (Float64, Float32)
    xb, vb = spline_bases(T)

    # the span holds 1, v and v², so the projection reproduces f
    f(x, v) = 1 + v / 3 + v^2 / 5
    dist = TensorSplineDistribution(xb, vb, f)
    @test eltype(dist.coefficients) === T
    @test size(dist.coefficients) == (16, 100)
    rng = Xoshiro(2)
    for _ in 1:20
        x, v = T(4π) * rand(rng, T), T(50) * rand(rng, T) - T(25)
        y = evaluate(dist, x, v)
        @test y isa T
        @test abs(y - f(x, v)) ≤ 64eps(T) * abs(f(x, v))
    end

    # a function of Float64 literals gives coefficients of the bases' element type
    @test eltype(TensorSplineDistribution(xb, vb, (x, v) -> 0.5 * exp(-v^2)).coefficients) ===
          T

    # the projection is Mx⁻¹ F Mv⁻¹ with F[i, j] = ∫∫ f Bᵢ Cⱼ
    t = spline_tables(xb, vb)
    g = landau_state(T)
    P = TensorSplineDistribution(xb, vb, g).coefficients
    @test norm(P - spline_project(t, g)) ≤ 64eps(T) * norm(P)

    # evaluate is Σᵢⱼ A[i, j] Bᵢ(x) Cⱼ(v), x before v, on a state that is not symmetric
    A = randn(Xoshiro(3), T, 16, 100)
    D = TensorSplineDistribution(xb, vb, A)
    for (x, v) in ((T(1.3), T(-7.1)), (T(9.2), T(3.4)), (T(0.4), T(24.5)))
        ref = sum(A[i, j] * evaluate(xb, i, x) * evaluate(vb, j, v)
        for i in 1:16, j in 1:100)
        y = evaluate(D, x, v)
        @test y isa T
        @test abs(y - ref) ≤ 64eps(T) * sum(abs, A) / 100
        # x and v swapped give another value
        @test abs(evaluate(D, v + T(25), x) - ref) > 1000eps(T) * sum(abs, A) / 100
    end
end
