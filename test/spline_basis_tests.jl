using BSplineKit
using LinearAlgebra
using QuadratureRules
using VlasovMethods
using Test

using VlasovMethods: evaluate_basis, evaluate_basis_derivative


@testset "1D B-spline basis evaluation with natural boundary conditions" begin

    s = SplineND(1, 2, 0:0.5:1, :Natural, GaussLegendreQuadrature(1))

    @test evaluate_basis(s, [0.0], 1) == 1.0
    @test evaluate_basis(s, [0.0], 2) == 0.0
    @test evaluate_basis(s, [0.0], 3) == 0.0

    @test evaluate_basis(s, [0.5], 1) == 0.0
    @test evaluate_basis(s, [0.5], 2) == 1.0
    @test evaluate_basis(s, [0.5], 3) == 0.0

    @test evaluate_basis(s, [1.0], 1) == 0.0
    @test evaluate_basis(s, [1.0], 2) == 0.0
    @test evaluate_basis(s, [1.0], 3) == 1.0


    s = SplineND(1, 3, 0:0.5:1, :Natural, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.0], 1) == 1.0
    @test evaluate_basis(s, [0.0], 2) == 0.0
    @test evaluate_basis(s, [0.0], 3) == 0.0
    @test evaluate_basis(s, [0.0], 4) == 0.0

    @test evaluate_basis(s, [0.5], 1) == 0.0
    @test evaluate_basis(s, [0.5], 2) == 0.5
    @test evaluate_basis(s, [0.5], 3) == 0.5
    @test evaluate_basis(s, [0.5], 4) == 0.0

    @test evaluate_basis(s, [1.0], 1) == 0.0
    @test evaluate_basis(s, [1.0], 2) == 0.0
    @test evaluate_basis(s, [1.0], 3) == 0.0
    @test evaluate_basis(s, [1.0], 4) == 1.0

end


@testset "1D B-spline basis evaluation with Dirichlet boundary conditions" begin

    s = SplineND(1, 2, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.0 ], 1) == 0.0
    @test evaluate_basis(s, [0.0 ], 2) == 0.0
    @test evaluate_basis(s, [0.0 ], 3) == 0.0

    @test evaluate_basis(s, [0.25], 1) == 1.0
    @test evaluate_basis(s, [0.25], 2) == 0.0
    @test evaluate_basis(s, [0.25], 3) == 0.0

    @test evaluate_basis(s, [0.5 ], 1) == 0.0
    @test evaluate_basis(s, [0.5 ], 2) == 1.0
    @test evaluate_basis(s, [0.5 ], 3) == 0.0

    @test evaluate_basis(s, [0.75], 1) == 0.0
    @test evaluate_basis(s, [0.75], 2) == 0.0
    @test evaluate_basis(s, [0.75], 3) == 1.0

    @test evaluate_basis(s, [1.0 ], 1) == 0.0
    @test evaluate_basis(s, [1.0 ], 2) == 0.0
    @test evaluate_basis(s, [1.0 ], 3) == 0.0


    s = SplineND(1, 3, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.0 ], 1) == 0.0
    @test evaluate_basis(s, [0.0 ], 2) == 0.0
    @test evaluate_basis(s, [0.0 ], 3) == 0.0
    @test evaluate_basis(s, [0.0 ], 4) == 0.0

    @test evaluate_basis(s, [0.25], 1) == 0.5
    @test evaluate_basis(s, [0.25], 2) == 0.5
    @test evaluate_basis(s, [0.25], 3) == 0.0
    @test evaluate_basis(s, [0.25], 4) == 0.0

    @test evaluate_basis(s, [0.5 ], 1) == 0.0
    @test evaluate_basis(s, [0.5 ], 2) == 0.5
    @test evaluate_basis(s, [0.5 ], 3) == 0.5
    @test evaluate_basis(s, [0.5 ], 4) == 0.0

    @test evaluate_basis(s, [0.75], 1) == 0.0
    @test evaluate_basis(s, [0.75], 2) == 0.0
    @test evaluate_basis(s, [0.75], 3) == 0.5
    @test evaluate_basis(s, [0.75], 4) == 0.5

    @test evaluate_basis(s, [1.0 ], 1) == 0.0
    @test evaluate_basis(s, [1.0 ], 2) == 0.0
    @test evaluate_basis(s, [1.0 ], 3) == 0.0
    @test evaluate_basis(s, [1.0 ], 4) == 0.0

end


@testset "1D B-spline basis evaluation with periodic boundary conditions" begin

    s = SplineND(1, 2, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.0 ], 1) == 1.0
    @test evaluate_basis(s, [0.0 ], 2) == 0.0
    @test evaluate_basis(s, [0.0 ], 3) == 0.0
    @test evaluate_basis(s, [0.0 ], 4) == 0.0

    @test evaluate_basis(s, [0.25], 1) == 0.0
    @test evaluate_basis(s, [0.25], 2) == 1.0
    @test evaluate_basis(s, [0.25], 3) == 0.0
    @test evaluate_basis(s, [0.25], 4) == 0.0

    @test evaluate_basis(s, [0.5 ], 1) == 0.0
    @test evaluate_basis(s, [0.5 ], 2) == 0.0
    @test evaluate_basis(s, [0.5 ], 3) == 1.0
    @test evaluate_basis(s, [0.5 ], 4) == 0.0

    @test evaluate_basis(s, [0.75], 1) == 0.0
    @test evaluate_basis(s, [0.75], 2) == 0.0
    @test evaluate_basis(s, [0.75], 3) == 0.0
    @test evaluate_basis(s, [0.75], 4) == 1.0

    @test evaluate_basis(s, [1.0 ], 1) == 0.0
    @test evaluate_basis(s, [1.0 ], 2) == 0.0
    @test evaluate_basis(s, [1.0 ], 3) == 0.0
    @test evaluate_basis(s, [1.0 ], 4) == 0.0


    s = SplineND(1, 3, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.0 ], 1) == 0.5
    @test evaluate_basis(s, [0.0 ], 2) == 0.0
    @test evaluate_basis(s, [0.0 ], 3) == 0.0
    @test evaluate_basis(s, [0.0 ], 4) == 0.0

    @test evaluate_basis(s, [0.25], 1) == 0.5
    @test evaluate_basis(s, [0.25], 2) == 0.5
    @test evaluate_basis(s, [0.25], 3) == 0.0
    @test evaluate_basis(s, [0.25], 4) == 0.0

    @test evaluate_basis(s, [0.5 ], 1) == 0.0
    @test evaluate_basis(s, [0.5 ], 2) == 0.5
    @test evaluate_basis(s, [0.5 ], 3) == 0.5
    @test evaluate_basis(s, [0.5 ], 4) == 0.0

    @test evaluate_basis(s, [0.75], 1) == 0.0
    @test evaluate_basis(s, [0.75], 2) == 0.0
    @test evaluate_basis(s, [0.75], 3) == 0.5
    @test evaluate_basis(s, [0.75], 4) == 0.5

    @test evaluate_basis(s, [1.0 ], 1) == 0.0
    @test evaluate_basis(s, [1.0 ], 2) == 0.0
    @test evaluate_basis(s, [1.0 ], 3) == 0.0
    @test evaluate_basis(s, [1.0 ], 4) == 0.5

end


@testset "2D B-spline basis evaluation with natural boundary conditions" begin

    s = SplineND(2, 2, 0:0.5:1, :Natural, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.0, 0.0], (1,1)) == 1.0
    @test evaluate_basis(s, [0.0, 0.0], (1,2)) == 0.0
    @test evaluate_basis(s, [0.0, 0.0], (1,3)) == 0.0

    @test evaluate_basis(s, [0.0, 0.5], (1,1)) == 0.0
    @test evaluate_basis(s, [0.0, 0.5], (1,2)) == 1.0
    @test evaluate_basis(s, [0.0, 0.5], (1,3)) == 0.0

    @test evaluate_basis(s, [0.0, 1.0], (1,1)) == 0.0
    @test evaluate_basis(s, [0.0, 1.0], (1,2)) == 0.0
    @test evaluate_basis(s, [0.0, 1.0], (1,3)) == 1.0

    @test evaluate_basis(s, [0.5, 0.0], (2,1)) == 1.0
    @test evaluate_basis(s, [0.5, 0.0], (2,2)) == 0.0
    @test evaluate_basis(s, [0.5, 0.0], (2,3)) == 0.0

    @test evaluate_basis(s, [0.5, 0.5], (2,1)) == 0.0
    @test evaluate_basis(s, [0.5, 0.5], (2,2)) == 1.0
    @test evaluate_basis(s, [0.5, 0.5], (2,3)) == 0.0

    @test evaluate_basis(s, [0.5, 1.0], (2,1)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (2,2)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (2,3)) == 1.0

    @test evaluate_basis(s, [1.0, 0.0], (3,1)) == 1.0
    @test evaluate_basis(s, [1.0, 0.0], (3,2)) == 0.0
    @test evaluate_basis(s, [1.0, 0.0], (3,3)) == 0.0

    @test evaluate_basis(s, [1.0, 0.5], (3,1)) == 0.0
    @test evaluate_basis(s, [1.0, 0.5], (3,2)) == 1.0
    @test evaluate_basis(s, [1.0, 0.5], (3,3)) == 0.0

    @test evaluate_basis(s, [1.0, 1.0], (3,1)) == 0.0
    @test evaluate_basis(s, [1.0, 1.0], (3,2)) == 0.0
    @test evaluate_basis(s, [1.0, 1.0], (3,3)) == 1.0

    @test evaluate_basis(s, [0.0, 0.5], (2,1)) == 0.0
    @test evaluate_basis(s, [0.0, 1.0], (3,1)) == 0.0

    @test evaluate_basis(s, [0.5, 0.0], (1,2)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (3,2)) == 0.0

    @test evaluate_basis(s, [1.0, 0.0], (1,3)) == 0.0
    @test evaluate_basis(s, [1.0, 0.5], (2,3)) == 0.0


    s = SplineND(2, 3, 0:0.5:1, :Natural, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.0, 0.0], (1,1)) == 1.0
    @test evaluate_basis(s, [0.0, 0.0], (1,2)) == 0.0
    @test evaluate_basis(s, [0.0, 0.0], (1,3)) == 0.0
    @test evaluate_basis(s, [0.0, 0.0], (1,4)) == 0.0

    @test evaluate_basis(s, [0.0, 0.5], (1,1)) == 0.0
    @test evaluate_basis(s, [0.0, 0.5], (1,2)) == 0.5
    @test evaluate_basis(s, [0.0, 0.5], (1,3)) == 0.5
    @test evaluate_basis(s, [0.0, 0.5], (1,4)) == 0.0

    @test evaluate_basis(s, [0.0, 1.0], (1,1)) == 0.0
    @test evaluate_basis(s, [0.0, 1.0], (1,2)) == 0.0
    @test evaluate_basis(s, [0.0, 1.0], (1,3)) == 0.0
    @test evaluate_basis(s, [0.0, 1.0], (1,4)) == 1.0
    
    @test evaluate_basis(s, [0.5, 0.0], (2,1)) == 0.5
    @test evaluate_basis(s, [0.5, 0.0], (2,2)) == 0.0
    @test evaluate_basis(s, [0.5, 0.0], (2,3)) == 0.0
    @test evaluate_basis(s, [0.5, 0.0], (2,4)) == 0.0

    @test evaluate_basis(s, [0.5, 0.5], (2,1)) == 0.0
    @test evaluate_basis(s, [0.5, 0.5], (2,2)) == 0.25
    @test evaluate_basis(s, [0.5, 0.5], (2,3)) == 0.25
    @test evaluate_basis(s, [0.5, 0.5], (2,4)) == 0.0

    @test evaluate_basis(s, [0.5, 1.0], (2,1)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (2,2)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (2,3)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (2,4)) == 0.5

    @test evaluate_basis(s, [0.5, 0.0], (3,1)) == 0.5
    @test evaluate_basis(s, [0.5, 0.0], (3,2)) == 0.0
    @test evaluate_basis(s, [0.5, 0.0], (3,3)) == 0.0
    @test evaluate_basis(s, [0.5, 0.0], (3,4)) == 0.0

    @test evaluate_basis(s, [0.5, 0.5], (3,1)) == 0.0
    @test evaluate_basis(s, [0.5, 0.5], (3,2)) == 0.25
    @test evaluate_basis(s, [0.5, 0.5], (3,3)) == 0.25
    @test evaluate_basis(s, [0.5, 0.5], (3,4)) == 0.0

    @test evaluate_basis(s, [0.5, 1.0], (3,1)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (3,2)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (3,3)) == 0.0
    @test evaluate_basis(s, [0.5, 1.0], (3,4)) == 0.5

    @test evaluate_basis(s, [1.0, 0.0], (4,1)) == 1.0
    @test evaluate_basis(s, [1.0, 0.0], (4,2)) == 0.0
    @test evaluate_basis(s, [1.0, 0.0], (4,3)) == 0.0
    @test evaluate_basis(s, [1.0, 0.0], (4,4)) == 0.0

    @test evaluate_basis(s, [1.0, 0.5], (4,1)) == 0.0
    @test evaluate_basis(s, [1.0, 0.5], (4,2)) == 0.5
    @test evaluate_basis(s, [1.0, 0.5], (4,3)) == 0.5
    @test evaluate_basis(s, [1.0, 0.5], (4,4)) == 0.0

    @test evaluate_basis(s, [1.0, 1.0], (4,1)) == 0.0
    @test evaluate_basis(s, [1.0, 1.0], (4,2)) == 0.0
    @test evaluate_basis(s, [1.0, 1.0], (4,3)) == 0.0
    @test evaluate_basis(s, [1.0, 1.0], (4,4)) == 1.0

end


@testset "2D B-spline basis evaluation with Dirichlet boundary conditions" begin

    s = SplineND(2, 2, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.00, 0.0 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 0.0 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 0.0 ], (1,3)) == 0.0

    @test evaluate_basis(s, [0.00, 0.25], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 0.25], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 0.25], (1,3)) == 0.0

    @test evaluate_basis(s, [0.00, 0.5 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 0.5 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 0.5 ], (1,3)) == 0.0

    @test evaluate_basis(s, [0.00, 0.75], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 0.75], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 0.75], (1,3)) == 0.0

    @test evaluate_basis(s, [0.00, 1.0 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 1.0 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 1.0 ], (1,3)) == 0.0

    @test evaluate_basis(s, [0.25, 0.0 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.0 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.0 ], (1,3)) == 0.0

    @test evaluate_basis(s, [0.25, 0.25], (1,1)) == 1.0
    @test evaluate_basis(s, [0.25, 0.25], (1,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.25], (1,3)) == 0.0

    @test evaluate_basis(s, [0.25, 0.5 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.5 ], (1,2)) == 1.0
    @test evaluate_basis(s, [0.25, 0.5 ], (1,3)) == 0.0

    @test evaluate_basis(s, [0.25, 0.75], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (1,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (1,3)) == 1.0

    @test evaluate_basis(s, [0.25, 1.0 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 1.0 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.25, 1.0 ], (1,3)) == 0.0

    @test evaluate_basis(s, [0.50, 0.0 ], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.0 ], (2,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.0 ], (2,3)) == 0.0

    @test evaluate_basis(s, [0.50, 0.25], (2,1)) == 1.0
    @test evaluate_basis(s, [0.50, 0.25], (2,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.25], (2,3)) == 0.0

    @test evaluate_basis(s, [0.50, 0.5 ], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.5 ], (2,2)) == 1.0
    @test evaluate_basis(s, [0.50, 0.5 ], (2,3)) == 0.0

    @test evaluate_basis(s, [0.50, 0.75], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (2,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (2,3)) == 1.0

    @test evaluate_basis(s, [0.50, 1.0 ], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 1.0 ], (2,2)) == 0.0
    @test evaluate_basis(s, [0.50, 1.0 ], (2,3)) == 0.0

    @test evaluate_basis(s, [0.75, 0.0 ], (3,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.0 ], (3,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.0 ], (3,3)) == 0.0

    @test evaluate_basis(s, [0.75, 0.25], (3,1)) == 1.0
    @test evaluate_basis(s, [0.75, 0.25], (3,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.25], (3,3)) == 0.0

    @test evaluate_basis(s, [0.75, 0.5 ], (3,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.5 ], (3,2)) == 1.0
    @test evaluate_basis(s, [0.75, 0.5 ], (3,3)) == 0.0

    @test evaluate_basis(s, [0.75, 0.75], (3,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (3,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (3,3)) == 1.0

    @test evaluate_basis(s, [0.75, 1.0 ], (3,1)) == 0.0
    @test evaluate_basis(s, [0.75, 1.0 ], (3,2)) == 0.0
    @test evaluate_basis(s, [0.75, 1.0 ], (3,3)) == 0.0

    @test evaluate_basis(s, [1.00, 0.0 ], (3,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.0 ], (3,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.0 ], (3,3)) == 0.0

    @test evaluate_basis(s, [1.00, 0.25], (3,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.25], (3,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.25], (3,3)) == 0.0

    @test evaluate_basis(s, [1.00, 0.5 ], (3,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.5 ], (3,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.5 ], (3,3)) == 0.0

    @test evaluate_basis(s, [1.00, 0.75], (3,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.75], (3,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.75], (3,3)) == 0.0

    @test evaluate_basis(s, [1.00, 1.0 ], (3,1)) == 0.0
    @test evaluate_basis(s, [1.00, 1.0 ], (3,2)) == 0.0
    @test evaluate_basis(s, [1.00, 1.0 ], (3,3)) == 0.0


    s = SplineND(2, 3, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.25, 0.25], (1,1)) == 0.25
    @test evaluate_basis(s, [0.25, 0.25], (1,2)) == 0.25
    @test evaluate_basis(s, [0.25, 0.25], (1,3)) == 0.0
    @test evaluate_basis(s, [0.25, 0.25], (1,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.5 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.5 ], (1,2)) == 0.25
    @test evaluate_basis(s, [0.25, 0.5 ], (1,3)) == 0.25
    @test evaluate_basis(s, [0.25, 0.5 ], (1,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.75], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (1,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (1,3)) == 0.25
    @test evaluate_basis(s, [0.25, 0.75], (1,4)) == 0.25

    @test evaluate_basis(s, [0.50, 0.25], (2,1)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (2,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (2,3)) == 0.0
    @test evaluate_basis(s, [0.50, 0.25], (2,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.5 ], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.5 ], (2,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (2,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (2,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.75], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (2,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (2,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.75], (2,4)) == 0.25

    @test evaluate_basis(s, [0.50, 0.25], (3,1)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (3,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (3,3)) == 0.0
    @test evaluate_basis(s, [0.50, 0.25], (3,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.5 ], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.5 ], (3,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (3,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (3,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.75], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (3,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (3,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.75], (3,4)) == 0.25

    @test evaluate_basis(s, [0.75, 0.25], (4,1)) == 0.25
    @test evaluate_basis(s, [0.75, 0.25], (4,2)) == 0.25
    @test evaluate_basis(s, [0.75, 0.25], (4,3)) == 0.0
    @test evaluate_basis(s, [0.75, 0.25], (4,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.5 ], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.5 ], (4,2)) == 0.25
    @test evaluate_basis(s, [0.75, 0.5 ], (4,3)) == 0.25
    @test evaluate_basis(s, [0.75, 0.5 ], (4,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.75], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (4,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (4,3)) == 0.25
    @test evaluate_basis(s, [0.75, 0.75], (4,4)) == 0.25

end


@testset "2D B-spline basis evaluation with periodic boundary conditions" begin

    s = SplineND(2, 2, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.00, 0.0 ], (1,1)) == 1.0
    @test evaluate_basis(s, [0.00, 0.0 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 0.0 ], (1,3)) == 0.0
    @test evaluate_basis(s, [0.00, 0.0 ], (1,4)) == 0.0

    @test evaluate_basis(s, [0.00, 0.25], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 0.25], (1,2)) == 1.0
    @test evaluate_basis(s, [0.00, 0.25], (1,3)) == 0.0
    @test evaluate_basis(s, [0.00, 0.25], (1,4)) == 0.0

    @test evaluate_basis(s, [0.00, 0.5 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 0.5 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 0.5 ], (1,3)) == 1.0
    @test evaluate_basis(s, [0.00, 0.5 ], (1,4)) == 0.0

    @test evaluate_basis(s, [0.00, 0.75], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 0.75], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 0.75], (1,3)) == 0.0
    @test evaluate_basis(s, [0.00, 0.75], (1,4)) == 1.0

    @test evaluate_basis(s, [0.00, 1.0 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.00, 1.0 ], (1,2)) == 0.0
    @test evaluate_basis(s, [0.00, 1.0 ], (1,3)) == 0.0
    @test evaluate_basis(s, [0.00, 1.0 ], (1,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.0 ], (2,1)) == 1.0
    @test evaluate_basis(s, [0.25, 0.0 ], (2,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.0 ], (2,3)) == 0.0
    @test evaluate_basis(s, [0.25, 0.0 ], (2,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.25], (2,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.25], (2,2)) == 1.0
    @test evaluate_basis(s, [0.25, 0.25], (2,3)) == 0.0
    @test evaluate_basis(s, [0.25, 0.25], (2,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.5 ], (2,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.5 ], (2,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.5 ], (2,3)) == 1.0
    @test evaluate_basis(s, [0.25, 0.5 ], (2,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.75], (2,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (2,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (2,3)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (2,4)) == 1.0

    @test evaluate_basis(s, [0.25, 1.0 ], (2,1)) == 0.0
    @test evaluate_basis(s, [0.25, 1.0 ], (2,2)) == 0.0
    @test evaluate_basis(s, [0.25, 1.0 ], (2,3)) == 0.0
    @test evaluate_basis(s, [0.25, 1.0 ], (2,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.0 ], (3,1)) == 1.0
    @test evaluate_basis(s, [0.50, 0.0 ], (3,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.0 ], (3,3)) == 0.0
    @test evaluate_basis(s, [0.50, 0.0 ], (3,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.25], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.25], (3,2)) == 1.0
    @test evaluate_basis(s, [0.50, 0.25], (3,3)) == 0.0
    @test evaluate_basis(s, [0.50, 0.25], (3,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.5 ], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.5 ], (3,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.5 ], (3,3)) == 1.0
    @test evaluate_basis(s, [0.50, 0.5 ], (3,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.75], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (3,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (3,3)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (3,4)) == 1.0

    @test evaluate_basis(s, [0.50, 1.0 ], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 1.0 ], (3,2)) == 0.0
    @test evaluate_basis(s, [0.50, 1.0 ], (3,3)) == 0.0
    @test evaluate_basis(s, [0.50, 1.0 ], (3,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.0 ], (4,1)) == 1.0
    @test evaluate_basis(s, [0.75, 0.0 ], (4,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.0 ], (4,3)) == 0.0
    @test evaluate_basis(s, [0.75, 0.0 ], (4,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.25], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.25], (4,2)) == 1.0
    @test evaluate_basis(s, [0.75, 0.25], (4,3)) == 0.0
    @test evaluate_basis(s, [0.75, 0.25], (4,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.5 ], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.5 ], (4,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.5 ], (4,3)) == 1.0
    @test evaluate_basis(s, [0.75, 0.5 ], (4,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.75], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (4,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (4,3)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (4,4)) == 1.0

    @test evaluate_basis(s, [0.75, 1.0 ], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 1.0 ], (4,2)) == 0.0
    @test evaluate_basis(s, [0.75, 1.0 ], (4,3)) == 0.0
    @test evaluate_basis(s, [0.75, 1.0 ], (4,4)) == 0.0

    @test evaluate_basis(s, [1.00, 0.0 ], (1,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.0 ], (1,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.0 ], (1,3)) == 0.0
    @test evaluate_basis(s, [1.00, 0.0 ], (1,4)) == 0.0

    @test evaluate_basis(s, [1.00, 0.25], (2,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.25], (2,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.25], (2,3)) == 0.0
    @test evaluate_basis(s, [1.00, 0.25], (2,4)) == 0.0

    @test evaluate_basis(s, [1.00, 0.5 ], (3,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.5 ], (3,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.5 ], (3,3)) == 0.0
    @test evaluate_basis(s, [1.00, 0.5 ], (3,4)) == 0.0

    @test evaluate_basis(s, [1.00, 0.75], (4,1)) == 0.0
    @test evaluate_basis(s, [1.00, 0.75], (4,2)) == 0.0
    @test evaluate_basis(s, [1.00, 0.75], (4,3)) == 0.0
    @test evaluate_basis(s, [1.00, 0.75], (4,4)) == 0.0

    @test evaluate_basis(s, [1.00, 1.0 ], (1,1)) == 0.0
    @test evaluate_basis(s, [1.00, 1.0 ], (2,2)) == 0.0
    @test evaluate_basis(s, [1.00, 1.0 ], (3,3)) == 0.0
    @test evaluate_basis(s, [1.00, 1.0 ], (4,4)) == 0.0


    s = SplineND(2, 3, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis(s, [0.25, 0.25], (1,1)) == 0.25
    @test evaluate_basis(s, [0.25, 0.25], (1,2)) == 0.25
    @test evaluate_basis(s, [0.25, 0.25], (1,3)) == 0.0
    @test evaluate_basis(s, [0.25, 0.25], (1,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.5 ], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.5 ], (1,2)) == 0.25
    @test evaluate_basis(s, [0.25, 0.5 ], (1,3)) == 0.25
    @test evaluate_basis(s, [0.25, 0.5 ], (1,4)) == 0.0

    @test evaluate_basis(s, [0.25, 0.75], (1,1)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (1,2)) == 0.0
    @test evaluate_basis(s, [0.25, 0.75], (1,3)) == 0.25
    @test evaluate_basis(s, [0.25, 0.75], (1,4)) == 0.25

    @test evaluate_basis(s, [0.50, 0.25], (2,1)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (2,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (2,3)) == 0.0
    @test evaluate_basis(s, [0.50, 0.25], (2,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.5 ], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.5 ], (2,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (2,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (2,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.75], (2,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (2,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (2,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.75], (2,4)) == 0.25

    @test evaluate_basis(s, [0.50, 0.25], (3,1)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (3,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.25], (3,3)) == 0.0
    @test evaluate_basis(s, [0.50, 0.25], (3,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.5 ], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.5 ], (3,2)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (3,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.5 ], (3,4)) == 0.0

    @test evaluate_basis(s, [0.50, 0.75], (3,1)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (3,2)) == 0.0
    @test evaluate_basis(s, [0.50, 0.75], (3,3)) == 0.25
    @test evaluate_basis(s, [0.50, 0.75], (3,4)) == 0.25

    @test evaluate_basis(s, [0.75, 0.25], (4,1)) == 0.25
    @test evaluate_basis(s, [0.75, 0.25], (4,2)) == 0.25
    @test evaluate_basis(s, [0.75, 0.25], (4,3)) == 0.0
    @test evaluate_basis(s, [0.75, 0.25], (4,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.5 ], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.5 ], (4,2)) == 0.25
    @test evaluate_basis(s, [0.75, 0.5 ], (4,3)) == 0.25
    @test evaluate_basis(s, [0.75, 0.5 ], (4,4)) == 0.0

    @test evaluate_basis(s, [0.75, 0.75], (4,1)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (4,2)) == 0.0
    @test evaluate_basis(s, [0.75, 0.75], (4,3)) == 0.25
    @test evaluate_basis(s, [0.75, 0.75], (4,4)) == 0.25

end



@testset "1D B-spline basis evaluation with natural boundary conditions" begin

    s = SplineND(1, 2, 0:0.5:1, :Natural, GaussLegendreQuadrature(1))

    @test evaluate_basis_derivative(s, [0.0], 1) == [-2.0]
    @test evaluate_basis_derivative(s, [0.0], 2) == [+2.0]
    @test evaluate_basis_derivative(s, [0.0], 3) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.5], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.5], 2) == [-2.0]
    @test evaluate_basis_derivative(s, [0.5], 3) == [+2.0]

    @test evaluate_basis_derivative(s, [1.0], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0], 2) == [-2.0]
    @test evaluate_basis_derivative(s, [1.0], 3) == [+2.0]


    s = SplineND(1, 3, 0:0.5:1, :Natural, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0], 1) == [-4.0]
    @test evaluate_basis_derivative(s, [0.0], 2) == [+4.0]
    @test evaluate_basis_derivative(s, [0.0], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.0], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.5], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.5], 2) == [-2.0]
    @test evaluate_basis_derivative(s, [0.5], 3) == [+2.0]
    @test evaluate_basis_derivative(s, [0.5], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [1.0], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0], 3) == [-4.0]
    @test evaluate_basis_derivative(s, [1.0], 4) == [+4.0]

end


@testset "1D B-spline basis evaluation with Dirichlet boundary conditions" begin

    s = SplineND(1, 2, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0 ], 1) == [+4.0]
    @test evaluate_basis_derivative(s, [0.0 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.0 ], 3) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.25], 1) == [-4.0]
    @test evaluate_basis_derivative(s, [0.25], 2) == [+4.0]
    @test evaluate_basis_derivative(s, [0.25], 3) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.5 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.5 ], 2) == [-4.0]
    @test evaluate_basis_derivative(s, [0.5 ], 3) == [+4.0]

    @test evaluate_basis_derivative(s, [0.75], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 3) == [-4.0]

    @test evaluate_basis_derivative(s, [1.0 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 3) == [-4.0]


    s = SplineND(1, 3, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0 ], 1) == [+8.0]
    @test evaluate_basis_derivative(s, [0.0 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.0 ], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.0 ], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.25], 1) == [-4.0]
    @test evaluate_basis_derivative(s, [0.25], 2) == [+4.0]
    @test evaluate_basis_derivative(s, [0.25], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.25], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.5 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.5 ], 2) == [-4.0]
    @test evaluate_basis_derivative(s, [0.5 ], 3) == [+4.0]
    @test evaluate_basis_derivative(s, [0.5 ], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.75], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 3) == [-4.0]
    @test evaluate_basis_derivative(s, [0.75], 4) == [+4.0]

    @test evaluate_basis_derivative(s, [1.0 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 4) == [-8.0]

end


@testset "1D B-spline basis derivative evaluation with periodic boundary conditions" begin

    s = SplineND(1, 2, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0 ], 1) == [-4.0]
    @test evaluate_basis_derivative(s, [0.0 ], 2) == [+4.0]
    @test evaluate_basis_derivative(s, [0.0 ], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.0 ], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.25], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.25], 2) == [-4.0]
    @test evaluate_basis_derivative(s, [0.25], 3) == [+4.0]
    @test evaluate_basis_derivative(s, [0.25], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.5 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.5 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.5 ], 3) == [-4.0]
    @test evaluate_basis_derivative(s, [0.5 ], 4) == [+4.0]

    @test evaluate_basis_derivative(s, [0.75], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 4) == [-4.0]

    @test evaluate_basis_derivative(s, [1.0 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 4) == [-4.0]


    s = SplineND(1, 3, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0 ], 1) == [+4.0]
    @test evaluate_basis_derivative(s, [0.0 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.0 ], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.0 ], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.25], 1) == [-4.0]
    @test evaluate_basis_derivative(s, [0.25], 2) == [+4.0]
    @test evaluate_basis_derivative(s, [0.25], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.25], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.5 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.5 ], 2) == [-4.0]
    @test evaluate_basis_derivative(s, [0.5 ], 3) == [+4.0]
    @test evaluate_basis_derivative(s, [0.5 ], 4) == [ 0.0]

    @test evaluate_basis_derivative(s, [0.75], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [0.75], 3) == [-4.0]
    @test evaluate_basis_derivative(s, [0.75], 4) == [+4.0]

    @test evaluate_basis_derivative(s, [1.0 ], 1) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 2) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 3) == [ 0.0]
    @test evaluate_basis_derivative(s, [1.0 ], 4) == [-4.0]

end


@testset "2D B-spline basis derivative evaluation with natural boundary conditions" begin

    s = SplineND(2, 2, 0:0.5:1, :Natural, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0, 0.0], (1,1)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0], (1,2)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0], (2,1)) == [+2.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.5], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5], (1,2)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5], (1,3)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 1.0], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0], (1,2)) == [ 0.0, -2.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0], (1,3)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,1)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,2)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,1)) == [+2.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,2)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,3)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,2)) == [ 0.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,3)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.0], (3,1)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0], (3,2)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0], (2,1)) == [-2.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0], (1,3)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.5], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5], (3,2)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5], (3,3)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5], (1,3)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 1.0], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0], (3,2)) == [ 0.0, -2.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0], (3,3)) == [+2.0, +2.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0], (1,3)) == [ 0.0,  0.0]


    s = SplineND(2, 3, 0:0.5:1, :Natural, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0, 0.0], (1,1)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0], (1,2)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.5], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5], (1,2)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5], (1,3)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 1.0], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0], (1,3)) == [ 0.0, -4.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0], (1,4)) == [-4.0, +4.0]
    
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,1)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,2)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,2)) == [-1.0, -1.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,3)) == [-1.0, +1.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,3)) == [ 0.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,4)) == [-2.0, +2.0]

    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,1)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,2)) == [ 0.0, +2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,2)) == [+1.0, -1.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,3)) == [+1.0, +1.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,3)) == [ 0.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,4)) == [+2.0, +2.0]

    @test evaluate_basis_derivative(s, [1.0, 0.0], (4,1)) == [+4.0, -4.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0], (4,2)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.5], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5], (4,2)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5], (4,3)) == [+2.0, +2.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 1.0], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0], (4,3)) == [ 0.0, -4.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0], (4,4)) == [+4.0, +4.0]

end


@testset "2D B-spline basis derivative evaluation with Dirichlet boundary conditions" begin

    s = SplineND(2, 2, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.25], (1,1)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.25], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.25], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.25], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.25], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (1,2)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.75], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.75], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.75], (1,3)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.75], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.75], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (3,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.0 ], (2,1)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0 ], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0 ], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0 ], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.25], (2,1)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.25], (2,2)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.25], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.25], (3,1)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.25], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.5 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5 ], (2,2)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5 ], (2,3)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5 ], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.75], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.75], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.75], (2,3)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.75], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.75], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 1.0 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0 ], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0 ], (2,3)) == [ 0.0, -4.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0 ], (1,2)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (1,3)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.25], (3,1)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.25], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.25], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.25], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.25], (1,3)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (3,2)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (1,3)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.75], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.75], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.75], (3,3)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.75], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.75], (1,3)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (1,3)) == [ 0.0,  0.0]


    s = SplineND(2, 3, 0:0.25:1, :Dirichlet, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.0 ], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.25], (1,1)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.25], (1,2)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.25], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.25], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (1,2)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (1,3)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.5 ], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 0.75], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.75], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.75], (1,3)) == [+4.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 0.75], (1,4)) == [+4.0,  0.0]

    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.0, 1.0 ], (1,4)) == [ 0.0,  0.0]
    
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,1)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,2)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,3)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (2,4)) == [ 0.0, -4.0]

    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,1)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.0], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,2)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,3)) == [+2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.5, 0.5], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.5, 1.0], (3,4)) == [ 0.0, -4.0]

    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.0 ], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.25], (4,1)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.25], (4,2)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.25], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.25], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (4,2)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (4,3)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.5 ], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 0.75], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.75], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.75], (4,3)) == [-4.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 0.75], (4,4)) == [-4.0,  0.0]

    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.0, 1.0 ], (4,4)) == [ 0.0,  0.0]

end


@testset "2D B-spline basis derivative evaluation with periodic boundary conditions" begin

    s = SplineND(2, 2, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.00, 0.0 ], (1,1)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.00, 0.0 ], (1,2)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.00, 0.0 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 0.0 ], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.00, 0.25], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 0.25], (1,2)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.00, 0.25], (1,3)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.00, 0.25], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.00, 0.5 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 0.5 ], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 0.5 ], (1,3)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.00, 0.5 ], (1,4)) == [ 0.0, +4.0]

    @test evaluate_basis_derivative(s, [0.00, 0.75], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 0.75], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 0.75], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 0.75], (1,4)) == [-4.0, -4.0]

    @test evaluate_basis_derivative(s, [0.00, 1.0 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 1.0 ], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 1.0 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.00, 1.0 ], (1,4)) == [ 0.0, -4.0]

    @test evaluate_basis_derivative(s, [0.25, 0.0 ], (2,1)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.25, 0.0 ], (2,2)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.25, 0.0 ], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.0 ], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.25, 0.25], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.25], (2,2)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.25, 0.25], (2,3)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.25, 0.25], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (2,3)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (2,4)) == [ 0.0, +4.0]

    @test evaluate_basis_derivative(s, [0.25, 0.75], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.75], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.75], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.75], (2,4)) == [-4.0, -4.0]

    @test evaluate_basis_derivative(s, [0.25, 1.0 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 1.0 ], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 1.0 ], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 1.0 ], (2,4)) == [ 0.0, -4.0]

    @test evaluate_basis_derivative(s, [0.50, 0.0 ], (3,1)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.50, 0.0 ], (3,2)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.50, 0.0 ], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.0 ], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,2)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,3)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,3)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,4)) == [ 0.0, +4.0]

    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,4)) == [-4.0, -4.0]

    @test evaluate_basis_derivative(s, [0.50, 1.0 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 1.0 ], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 1.0 ], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 1.0 ], (3,4)) == [ 0.0, -4.0]

    @test evaluate_basis_derivative(s, [0.75, 0.0 ], (4,1)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.75, 0.0 ], (4,2)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.75, 0.0 ], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.0 ], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,2)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,3)) == [ 0.0, +4.0]
    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,3)) == [-4.0, -4.0]
    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,4)) == [ 0.0, +4.0]

    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,4)) == [-4.0, -4.0]

    @test evaluate_basis_derivative(s, [0.75, 1.0 ], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 1.0 ], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 1.0 ], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 1.0 ], (4,4)) == [ 0.0, -4.0]

    @test evaluate_basis_derivative(s, [1.00, 0.0 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.0 ], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.0 ], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.0 ], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.00, 0.25], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.25], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.25], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.25], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.00, 0.5 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.5 ], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.5 ], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.5 ], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [1.00, 0.75], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.75], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.75], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 0.75], (4,4)) == [-4.0,  0.0]

    @test evaluate_basis_derivative(s, [1.00, 1.0 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 1.0 ], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 1.0 ], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [1.00, 1.0 ], (4,4)) == [ 0.0,  0.0]


    s = SplineND(2, 3, 0:0.25:1, :Periodic, GaussLegendreQuadrature(2))

    @test evaluate_basis_derivative(s, [0.25, 0.25], (1,1)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.25, 0.25], (1,2)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.25, 0.25], (1,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.25], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (1,2)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (1,3)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.25, 0.5 ], (1,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.25, 0.75], (1,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.75], (1,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.25, 0.75], (1,3)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.25, 0.75], (1,4)) == [-2.0, +2.0]

    @test evaluate_basis_derivative(s, [0.50, 0.25], (2,1)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (2,2)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (2,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (2,2)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (2,3)) == [-2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (2,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.50, 0.75], (2,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (2,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (2,3)) == [-2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (2,4)) == [-2.0, +2.0]

    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,1)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,2)) == [+2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.25], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,2)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,3)) == [+2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.5 ], (3,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,3)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.50, 0.75], (3,4)) == [+2.0, +2.0]

    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,1)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,2)) == [+2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,3)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.25], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,2)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,3)) == [+2.0, +2.0]
    @test evaluate_basis_derivative(s, [0.75, 0.5 ], (4,4)) == [ 0.0,  0.0]

    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,1)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,2)) == [ 0.0,  0.0]
    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,3)) == [+2.0, -2.0]
    @test evaluate_basis_derivative(s, [0.75, 0.75], (4,4)) == [+2.0, +2.0]

end
