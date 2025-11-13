using BSplineKit
using LinearAlgebra
using QuadratureRules
using Random
using VlasovMethods
using Test

using VlasovMethods: cartesian_index, linear_index, evaluate_indices
using VlasovMethods: mass_matrix, order, remap_unit_interval, unique_knots


@testset "Spline Utilities" begin

    @test remap_unit_interval(0.0, 0.0, 2.0) == 0.0
    @test remap_unit_interval(0.5, 0.0, 2.0) == 1.0
    @test remap_unit_interval(1.0, 0.0, 2.0) == 2.0

end


@testset "Mass Matrix" begin

    # Knot vector with 5 knots
    nknots = 5
    sknots = collect(0.0:(nknots-1)) ./ 2

    # B-splines of order 2
    sorder = 2

    # B-spline basis with Dirichlet BCs
    sbasis = BSplineBasis(BSplineOrder(sorder), copy(sknots))
    nbasis = length(sbasis)

    # Dirichlet BCs with (exact) Gauss-Legendre quadrature
    smass = mass_matrix(sbasis, GaussLegendreQuadrature(2))
    rmass = galerkin_matrix(sbasis)

    @test all(isapprox.(smass, rmass; atol=2eps()))

    # Dirichlet BCs with trapezoidal quadrature
    smass = mass_matrix(sbasis, TrapezoidalQuadrature())
    rmass = Diagonal([0.5, repeat([1.0], nbasis - 2)..., 0.5] ./ 2)

    @test all(smass .== rmass)

    # B-spline basis with Periodic BCs
    sbasis = PeriodicBSplineBasis(BSplineOrder(sorder), copy(sknots))
    nbasis = length(sbasis)

    # Periodic BCs with (exact) Gauss-Legendre quadrature
    smass = mass_matrix(sbasis, GaussLegendreQuadrature(2))
    rmass = galerkin_matrix(sbasis)

    @test all(isapprox.(smass, rmass; atol=2eps()))

    # Periodic BCs with trapezoidal quadrature
    smass = mass_matrix(sbasis, TrapezoidalQuadrature())
    rmass = Matrix(1.0I, nbasis, nbasis) ./ 2

    @test all(smass .== rmass)


    # B-splines of order 3 to 6
    for sorder in 3:6

        # B-spline basis with Dirichlet BCs
        sbasis = BSplineBasis(BSplineOrder(sorder), copy(sknots))

        # Dirichlet BCs with (exact) Gauss-Legendre quadrature
        smass = mass_matrix(sbasis, GaussLegendreQuadrature(sorder))
        rmass = galerkin_matrix(sbasis)

        @test all(isapprox.(smass, rmass; atol=2eps()))


        # B-spline basis with Periodic BCs
        sbasis = PeriodicBSplineBasis(BSplineOrder(sorder), copy(sknots))

        # Periodic BCs with (exact) Gauss-Legendre quadrature
        smass = mass_matrix(sbasis, GaussLegendreQuadrature(sorder))
        rmass = galerkin_matrix(sbasis)

        @test all(isapprox.(smass, rmass; atol=2eps()))

    end

end


@testset "1-dimensional B-Spline" begin

    d = 1
    o = 2
    k = 0:0.1:2
    q = GaussLegendreQuadrature(2)

    l = -2:0.2:+2

    @test_nowarn SplineND(d, o, k, q)
    @test_nowarn SplineND(d, o, k, :Natural, q)
    @test_nowarn SplineND(d, o, k, :Dirichlet, q)
    @test_nowarn SplineND(d, o, k, :Periodic, q)


    ### B-spline with natural boundary conditions ###

    s = SplineND(d, o, k, :Natural, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k) + (o - 2),)
    @test length(s) == size(s)[1]
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin]) == s([k[begin]]) == s.coefficients[begin]
    @test s(k[end]) == s([k[end]]) == s.coefficients[end]

    @test evaluate_indices(s, [0.0]) == [2, 1]
    @test evaluate_indices(s, [0.1]) == [3, 2]
    @test evaluate_indices(s, [1.0]) == [12, 11]
    @test evaluate_indices(s, [1.89]) == [20, 19]
    @test evaluate_indices(s, [1.90]) == [21, 20]
    @test evaluate_indices(s, [1.99]) == [21, 20]
    @test evaluate_indices(s, [2.0]) == [21, 20]

    L2projection!(x -> 0.0, s)
    @test all(isapprox.(s.coefficients, 0; atol=8eps()))
    @test all(isapprox.(s.(k), 0; atol=8eps()))

    L2projection!(x -> 1.0, s)
    @test all(isapprox.(s.coefficients, 1; atol=8eps()))
    @test all(isapprox.(s.(k), 1; atol=8eps()))

    L2projection!(x -> x[1], s)
    @test all(isapprox.(s.(k), k; atol=8eps()))

    L2projection!(x -> sin(π * x[1]), s)
    @test all(isapprox.(s.(k), sin.(π * k); atol=1E-2))

    s = L2projection!(x -> sin(π * x[1]), SplineND(d, 10, k, :Natural, q))
    @test all(isapprox.(s.(k), sin.(π * k); atol=1E-7))


    ### B-spline with Dirichlet boundary conditions ###

    s = SplineND(d, o, k, :Dirichlet, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k) - 2,)
    @test length(s) == size(s)[1]
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin]) == s([k[begin]]) == 0
    @test s(k[end]) == s([k[end]]) == 0

    @test evaluate_indices(s, [0.0]) == [1, 0]
    @test evaluate_indices(s, [0.1]) == [2, 1]
    @test evaluate_indices(s, [1.0]) == [11, 10]
    @test evaluate_indices(s, [1.89]) == [19, 18]
    @test evaluate_indices(s, [1.90]) == [20, 19]
    @test evaluate_indices(s, [1.99]) == [20, 19]
    @test evaluate_indices(s, [2.0]) == [20, 19]

    L2projection!(x -> 0.0, s)
    @test all(isapprox.(s.coefficients, 0; atol=8eps()))
    @test all(isapprox.(s.(k), 0; atol=8eps()))

    s = L2projection!(x -> abs(x[1]) - 2, SplineND(d, o, l, :Dirichlet, q))
    @test all(isapprox.(s.(l), abs.(l) .- 2; atol=8eps()))

    s = L2projection!(x -> sin(π * x[1]), SplineND(d, o, k, :Dirichlet, q))
    @test all(isapprox.(s.(k), sin.(π * k); atol=1E-2))

    s = L2projection!(x -> sin(π * x[1]), SplineND(d, 10, k, :Dirichlet, q))
    @test all(isapprox.(s.(k), sin.(π * k); atol=1E-12))


    ### B-spline with periodic boundary conditions ###

    s = SplineND(d, o, k, :Periodic, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k) - 1,)
    @test length(s) == size(s)[1]
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin]) == s(k[end])

    @test evaluate_indices(s, [0.0]) == [2, 1]
    @test evaluate_indices(s, [0.1]) == [3, 2]
    @test evaluate_indices(s, [1.0]) == [12, 11]
    @test evaluate_indices(s, [1.89]) == [20, 19]
    @test evaluate_indices(s, [1.90]) == [1, 20]
    @test evaluate_indices(s, [1.99]) == [1, 20]
    @test evaluate_indices(s, [2.0]) == [2, 1]

    L2projection!(x -> 0.0, s)
    @test all(isapprox.(s.coefficients, 0; atol=8eps()))
    @test all(isapprox.(s.(k), 0; atol=8eps()))

    L2projection!(x -> 1.0, s)
    @test all(isapprox.(s.coefficients, 1; atol=8eps()))
    @test all(isapprox.(s.(k), 1; atol=8eps()))

    L2projection!(x -> sin(π * x[1]), s)
    @test all(isapprox.(s.(k), sin.(π * k); atol=1E-2))

    s = L2projection!(x -> sin(π * x[1]), SplineND(d, 10, k, :Periodic, q))
    @test all(isapprox.(s.(k), sin.(π * k); atol=1E-12))


    ### Check size, knots, and indices
    function check_indices(b)
        @test cartesian_index(b, 1) == (1,)
        @test cartesian_index(b, size(b)[1]) == size(b)
        @test linear_index(b, 1) == 1
        @test linear_index(b, size(b)...) == size(b)[1]
        @test cartesian_index(b, linear_index(b, size(b)...)) == size(b)
    end

    for o in 2:8
        b = SplineND(d, o, k, :Natural, q)
        @test size(b) == (length(k) + (o - 2),)
        @test unique_knots(b) == k
        check_indices(b)

        b = SplineND(d, o, k, :Dirichlet, q)
        @test size(b) == (length(k) + (o - 4),)
        @test unique_knots(b) == k
        check_indices(b)

        b = SplineND(d, o, k, :Periodic, q)
        @test size(b) == (length(k) - 1,)
        @test unique_knots(b) == k
        check_indices(b)
    end

end


@testset "2-dimensional B-Spline" begin

    d = 2
    o = 3
    k = -2:0.1:+1
    q = GaussLegendreQuadrature(3)

    x = vec(k * ones(axes(k))')
    y = vec(ones(axes(k)) * k')
    z = [[x[i], y[i]] for i in eachindex(x, y)]

    l = -2:0.1:+2
    m = vec(l * ones(axes(l))')
    n = vec(ones(axes(l)) * l')
    t = [[m[i], n[i]] for i in eachindex(m, n)]
    r = [m .* n for (m,n) in t]

    @test_nowarn SplineND(d, o, k, q)
    @test_nowarn SplineND(d, o, k, :Natural, q)
    @test_nowarn SplineND(d, o, k, :Dirichlet, q)
    @test_nowarn SplineND(d, o, k, :Periodic, q)

    ### B-spline with natural boundary conditions ###

    s = SplineND(d, o, k, :Natural, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k) + (o - 2), length(k) + (o - 2))
    @test length(s) == size(s)[1] * size(s)[2]
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin], k[begin]) == s.coefficients[begin, begin]
    @test s(k[begin], k[end]) == s.coefficients[begin, end]
    @test s(k[end], k[begin]) == s.coefficients[end, begin]
    @test s(k[end], k[end]) == s.coefficients[end, end]


    L2projection!(x -> 0.0, s)
    @test all(isapprox.(s.coefficients, 0; atol=8eps()))
    @test all(isapprox.(s.(z), 0; atol=8eps()))

    L2projection!(x -> 1.0, s)
    @test all(isapprox.(s.coefficients, 1; atol=32eps()))
    @test all(isapprox.(s.(z), 1; atol=16eps()))

    s = L2projection!(x -> x[1] * x[2], SplineND(d, 2, l, :Natural, q))
    @test all(isapprox.(s.(t), r; atol=32eps()))

    s = L2projection!(x -> sin(π * x[1]) * sin(π * x[2]), SplineND(d, o, l, :Natural, q))
    @test all(isapprox.(s.(t), sin.(π .* m) .* sin.(π .* n); atol=1E-2))



    ### B-spline with Dirichlet boundary conditions ###

    s = SplineND(d, o, k, :Dirichlet, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k) - 1, length(k) - 1)
    @test length(s) == size(s)[1] * size(s)[2]
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin], k[begin]) == 0
    @test s(k[begin], k[end]) == 0
    @test s(k[end], k[begin]) == 0
    @test s(k[end], k[end]) == 0

    L2projection!(x -> 0.0, s)
    @test all(isapprox.(s.coefficients, 0; atol=8eps()))
    @test all(isapprox.(s.(z), 0; atol=8eps()))

    s = L2projection!(x -> sin(π * x[1]) * sin(π * x[2]), SplineND(d, o, l, :Natural, q))
    @test all(isapprox.(s.(t), sin.(π .* m) .* sin.(π .* n); atol=1E-2))


    ### B-spline with periodic boundary conditions ###

    s = SplineND(d, o, k, :Periodic, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k) - 1, length(k) - 1)
    @test length(s) == size(s)[1] * size(s)[2]
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin], k[begin]) ≈ s(k[begin], k[end]) atol = 2eps()
    @test s(k[begin], k[begin]) ≈ s(k[end], k[begin]) atol = 2eps()
    @test s(k[begin], k[begin]) ≈ s(k[end], k[end]) atol = 2eps()

    L2projection!(x -> 0.0, s)
    @test all(isapprox.(s.coefficients, 0; atol=8eps()))
    @test all(isapprox.(s.(z), 0; atol=8eps()))

    L2projection!(x -> 1.0, s)
    @test all(isapprox.(s.coefficients, 1; atol=32eps()))
    @test all(isapprox.(s.(z), 1; atol=16eps()))

    s = L2projection!(x -> abs(x[1] * x[2]), SplineND(d, 2, l, :Periodic, q))
    @test all(isapprox.(s.(t), abs.(r); atol=32eps()))

    s = L2projection!(x -> sin(π * x[1]) * sin(π * x[2]), SplineND(d, o, l, :Periodic, q))
    L2projection!(x -> sin(π * x[1]) * sin(π * x[2]), s)
    @test all(isapprox.(s.(t), sin.(π .* m) .* sin.(π .* n); atol=1E-4))


    ### Check size and knots

    function check_indices(b)
        @test cartesian_index(b, 1) == (1, 1)
        @test cartesian_index(b, size(b)[1] * size(b)[2]) == size(b)
        @test linear_index(b, 1, 1) == 1
        @test linear_index(b, size(b)...) == size(b)[1] * size(b)[2]
        @test cartesian_index(b, linear_index(b, size(b)...)) == size(b)
    end

    for o in 2:8
        b = SplineND(d, o, k, :Natural, q)
        @test size(b) == (length(k) + (o - 2), length(k) + (o - 2))
        @test unique_knots(b) == k
        check_indices(b)

        b = SplineND(d, o, k, :Dirichlet, q)
        @test size(b) == (length(k) + (o - 4), length(k) + (o - 4))
        @test unique_knots(b) == k
        check_indices(b)

        b = SplineND(d, o, k, :Periodic, q)
        @test size(b) == (length(k) - 1, length(k) - 1)
        @test unique_knots(b) == k
        check_indices(b)
    end

end
