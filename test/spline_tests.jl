using BSplineKit
using LinearAlgebra
using QuadratureRules
using Random
using VlasovMethods
using Test

using VlasovMethods: cartesian_indices, mass_matrix, order, remap_unit_interval, unique_knots


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

    @test all(isapprox.(smass, rmass; atol = 2eps()))

    # Dirichlet BCs with trapezoidal quadrature
    smass = mass_matrix(sbasis, TrapezoidalQuadrature())
    rmass = Diagonal([0.5, repeat([1.0], nbasis-2)..., 0.5] ./ 2)

    @test all(smass .== rmass)

    # B-spline basis with Periodic BCs
    sbasis = PeriodicBSplineBasis(BSplineOrder(sorder), copy(sknots))
    nbasis = length(sbasis)

    # Periodic BCs with (exact) Gauss-Legendre quadrature
    smass = mass_matrix(sbasis, GaussLegendreQuadrature(2))
    rmass = galerkin_matrix(sbasis)

    @test all(isapprox.(smass, rmass; atol = 2eps()))

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

        @test all(isapprox.(smass, rmass; atol = 2eps()))


        # B-spline basis with Periodic BCs
        sbasis = PeriodicBSplineBasis(BSplineOrder(sorder), copy(sknots))
        
        # Periodic BCs with (exact) Gauss-Legendre quadrature
        smass = mass_matrix(sbasis, GaussLegendreQuadrature(sorder))
        rmass = galerkin_matrix(sbasis)

        @test all(isapprox.(smass, rmass; atol = 2eps()))

    end

end


@testset "1-dimensional B-Spline" begin

    d = 1
    o = 2
    k = 0:0.1:2
    q = GaussLegendreQuadrature(2)

    @test_nowarn SplineND(d, o, k, q)
    @test_nowarn SplineND(d, o, k, :Natural, q)
    @test_nowarn SplineND(d, o, k, :Dirichlet, q)
    @test_nowarn SplineND(d, o, k, :Periodic, q)


    ### B-spline with natural boundary conditions ###

    s = SplineND(d, o, k, :Natural, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k)+(o-2),)
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin]) == s([k[begin]]) == s.coefficients[begin]
    @test s(k[end]) == s([k[end]]) == s.coefficients[end]


    ### B-spline with Dirichlet boundary conditions ###

    s = SplineND(d, o, k, :Dirichlet, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k)-2,)
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin]) == s([k[begin]]) == 0
    @test s(k[end]) == s([k[end]]) == 0


    ### B-spline with periodic boundary conditions ###

    s = SplineND(d, o, k, :Periodic, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k)-1,)
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin]) == s(k[end])

    
    ### Check size, knots, and indices

    for o in 2:8
        b = SplineND(d, o, k, :Natural, q)
        @test size(b) == (length(k) + (o-2),)
        @test unique_knots(b) == k
        @test cartesian_indices(b, 1) == (1,)
        @test cartesian_indices(b, size(b)[1]) == size(b)
            
        b = SplineND(d, o, k, :Dirichlet, q)
        @test size(b) == (length(k) + (o-4),)
        @test unique_knots(b) == k
        @test cartesian_indices(b, 1) == (1,)
        @test cartesian_indices(b, size(b)[1]) == size(b)
        
        b = SplineND(d, o, k, :Periodic, q)
        @test size(b) == (length(k)-1,)
        @test unique_knots(b) == k
        @test cartesian_indices(b, 1) == (1,)
        @test cartesian_indices(b, size(b)[1]) == size(b)
    end

end


@testset "2-dimensional B-Spline" begin

    d = 2
    o = 3
    k = -2:0.1:+1
    q = GaussLegendreQuadrature(3)

    @test_nowarn SplineND(d, o, k, q)
    @test_nowarn SplineND(d, o, k, :Natural, q)
    @test_nowarn SplineND(d, o, k, :Dirichlet, q)
    @test_nowarn SplineND(d, o, k, :Periodic, q)

    ### B-spline with natural boundary conditions ###

    s = SplineND(d, o, k, :Natural, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k)+(o-2), length(k)+(o-2))
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin], k[begin]) == s.coefficients[begin, begin]
    @test s(k[begin], k[end]) == s.coefficients[begin, end]
    @test s(k[end], k[begin]) == s.coefficients[end, begin]
    @test s(k[end], k[end]) == s.coefficients[end, end]


    ### B-spline with Dirichlet boundary conditions ###

    s = SplineND(d, o, k, :Dirichlet, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k)-1, length(k)-1)
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin], k[begin]) == 0
    @test s(k[begin], k[end]) == 0
    @test s(k[end], k[begin]) == 0
    @test s(k[end], k[end]) == 0


    ### B-spline with periodic boundary conditions ###

    s = SplineND(d, o, k, :Periodic, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k)-1, length(k)-1)
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin], k[begin]) ≈ s(k[begin], k[end]) atol = 2eps()
    @test s(k[begin], k[begin]) ≈ s(k[end], k[begin]) atol = 2eps()
    @test s(k[begin], k[begin]) ≈ s(k[end], k[end])   atol = 2eps()


    ### Check size and knots

    for o in 2:8
        b = SplineND(d, o, k, :Natural, q)
        @test size(b) == (length(k) + (o-2), length(k) + (o-2))
        @test unique_knots(b) == k
        @test cartesian_indices(b, 1) == (1,1)
        @test cartesian_indices(b, size(b)[1] * size(b)[2]) == size(b)
        
        b = SplineND(d, o, k, :Dirichlet, q)
        @test size(b) == (length(k) + (o-4), length(k) + (o-4))
        @test unique_knots(b) == k
        @test cartesian_indices(b, 1) == (1,1)
        @test cartesian_indices(b, size(b)[1] * size(b)[2]) == size(b)
        
        b = SplineND(d, o, k, :Periodic, q)
        @test size(b) == (length(k)-1, length(k)-1)
        @test unique_knots(b) == k
        @test cartesian_indices(b, 1) == (1,1)
        @test cartesian_indices(b, size(b)[1] * size(b)[2]) == size(b)
    end

end
