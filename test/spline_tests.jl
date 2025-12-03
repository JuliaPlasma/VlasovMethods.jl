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

    @test evaluate_indices(s, [0.0]) == [2,1]
    @test evaluate_indices(s, [0.1]) == [3,2]
    @test evaluate_indices(s, [1.0]) == [12,11]
    @test evaluate_indices(s, [1.89]) == [20,19]
    @test evaluate_indices(s, [1.90]) == [21,20]
    @test evaluate_indices(s, [1.99]) == [21,20]
    @test evaluate_indices(s, [2.0]) == [21,20]

    s = SplineND(d, 3, k, :Natural, q)

    @test evaluate_indices(s, [0.0]) == [3,2,1]
    @test evaluate_indices(s, [0.1]) == [4,3,2]
    @test evaluate_indices(s, [1.0]) == [13,12,11]
    @test evaluate_indices(s, [1.89]) == [21,20,19]
    @test evaluate_indices(s, [1.90]) == [22,21,20]
    @test evaluate_indices(s, [1.99]) == [22,21,20]
    @test evaluate_indices(s, [2.0]) == [22,21,20]

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

    @test evaluate_indices(s, [0.0]) == [1,0]
    @test evaluate_indices(s, [0.1]) == [2,1]
    @test evaluate_indices(s, [1.0]) == [11,10]
    @test evaluate_indices(s, [1.89]) == [19,18]
    @test evaluate_indices(s, [1.90]) == [20,19]
    @test evaluate_indices(s, [1.99]) == [20,19]
    @test evaluate_indices(s, [2.0]) == [20,19]

    s = SplineND(d, 3, k, :Dirichlet, q)

    @test evaluate_indices(s, [0.0]) == [2,1,0]
    @test evaluate_indices(s, [0.1]) == [3,2,1]
    @test evaluate_indices(s, [1.0]) == [12,11,10]
    @test evaluate_indices(s, [1.89]) == [20,19,18]
    @test evaluate_indices(s, [1.90]) == [21,20,19]
    @test evaluate_indices(s, [1.99]) == [21,20,19]
    @test evaluate_indices(s, [2.0]) == [21,20,19]

    ### B-spline with periodic boundary conditions ###

    s = SplineND(d, o, k, :Periodic, q)

    @test eltype(s) == Float64
    @test ndims(s) == d
    @test order(s) == o
    @test size(s) == (length(k)-1,)
    @test unique_knots(s) == k

    rand!(s.coefficients)

    @test s(k[begin]) == s(k[end])

    @test evaluate_indices(s, [0.0]) == [2,1]
    @test evaluate_indices(s, [0.1]) == [3,2]
    @test evaluate_indices(s, [1.0]) == [12,11]
    @test evaluate_indices(s, [1.89]) == [20,19]
    @test evaluate_indices(s, [1.90]) == [1,20]
    @test evaluate_indices(s, [1.99]) == [1,20]
    @test evaluate_indices(s, [2.0]) == [2,1]

    
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
        @test size(b) == (length(k) + (o-2),)
        @test unique_knots(b) == k
        check_indices(b)
            
        b = SplineND(d, o, k, :Dirichlet, q)
        @test size(b) == (length(k) + (o-4),)
        @test unique_knots(b) == k
        check_indices(b)
        
        b = SplineND(d, o, k, :Periodic, q)
        @test size(b) == (length(k)-1,)
        @test unique_knots(b) == k
        check_indices(b)
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


    @test evaluate_indices(s, [-2.,-2.]) == [67,66,65,35,34,33,3,2,1]
    @test evaluate_indices(s, [-2.,-1.]) == [387,386,385,355,354,353,323,322,321]
    @test evaluate_indices(s, [-2., 0.]) == [707,706,705,675,674,673,643,642,641]
    @test evaluate_indices(s, [-2.,+1.]) == [995,994,993,963,962,961,931,930,929]

    @test evaluate_indices(s, [-1.,-2.]) == [77,76,75,45,44,43,13,12,11]
    @test evaluate_indices(s, [-1.,-1.]) == [397,396,395,365,364,363,333,332,331]
    @test evaluate_indices(s, [-1., 0.]) == [717,716,715,685,684,683,653,652,651]
    @test evaluate_indices(s, [-1.,+1.]) == [1005,1004,1003,973,972,971,941,940,939]

    @test evaluate_indices(s, [ 0.,-2.]) == [87,86,85,55,54,53,23,22,21]
    @test evaluate_indices(s, [ 0.,-1.]) == [407,406,405,375,374,373,343,342,341]
    @test evaluate_indices(s, [ 0., 0.]) == [727,726,725,695,694,693,663,662,661]
    @test evaluate_indices(s, [ 0.,+1.]) == [1015,1014,1013,983,982,981,951,950,949]

    @test evaluate_indices(s, [+1.,-2.]) == [96,95,94,64,63,62,32,31,30]
    @test evaluate_indices(s, [+1.,-1.]) == [416,415,414,384,383,382,352,351,350]
    @test evaluate_indices(s, [+1., 0.]) == [736,735,734,704,703,702,672,671,670]
    @test evaluate_indices(s, [+1.,+1.]) == [1024,1023,1022,992,991,990,960,959,958]

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

    # @test evaluate_indices(s, [-2.,-2.]) == []
    # @test evaluate_indices(s, [-2.,-1.]) == []
    # @test evaluate_indices(s, [-2., 0.]) == []
    # @test evaluate_indices(s, [-2.,+1.]) == []

    # @test evaluate_indices(s, [-1.,-2.]) == []
    @test evaluate_indices(s, [-1.,-1.]) == [342,341,340,312,311,310,282,281,280]
    @test evaluate_indices(s, [-1., 0.]) == [642,641,640,612,611,610,582,581,580]
    @test evaluate_indices(s, [-1.,+1.]) == [912,911,910,882,881,880,852,851,850]

    # @test evaluate_indices(s, [ 0.,-2.]) == []
    @test evaluate_indices(s, [ 0.,-1.]) == [352,351,350,322,321,320,292,291,290]
    @test evaluate_indices(s, [ 0., 0.]) == [652,651,650,622,621,620,592,591,590]
    @test evaluate_indices(s, [ 0.,+1.]) == [922,921,920,892,891,890,862,861,860]

    # @test evaluate_indices(s, [+1.,-2.]) == []
    # @test evaluate_indices(s, [+1.,-1.]) == []
    # @test evaluate_indices(s, [+1., 0.]) == []
    # @test evaluate_indices(s, [+1.,+1.]) == []


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

    function check_indices(b)
        @test cartesian_index(b, 1) == (1,1)
        @test cartesian_index(b, size(b)[1] * size(b)[2]) == size(b)
        @test linear_index(b, 1, 1) == 1
        @test linear_index(b, size(b)...) == size(b)[1] * size(b)[2]
        @test cartesian_index(b, linear_index(b, size(b)...)) == size(b)
    end

    for o in 2:8
        b = SplineND(d, o, k, :Natural, q)
        @test size(b) == (length(k) + (o-2), length(k) + (o-2))
        @test unique_knots(b) == k
        check_indices(b)

        b = SplineND(d, o, k, :Dirichlet, q)
        @test size(b) == (length(k) + (o-4), length(k) + (o-4))
        @test unique_knots(b) == k
        check_indices(b)

        b = SplineND(d, o, k, :Periodic, q)
        @test size(b) == (length(k)-1, length(k)-1)
        @test unique_knots(b) == k
        check_indices(b)
    end

end
