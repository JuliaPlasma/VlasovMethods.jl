struct SplineDistribution{XD, VD, ST, BT, CT, MT} <: DistributionFunction{XD,VD}
    spline::ST
    basis::BT
    coefficients::CT
    mass_matrix::MT

    function SplineDistribution(xdim, vdim, basis::BT, coefficients::CT) where {BT, CT}
        spline = Spline(basis, coefficients)
        mass_matrix = galerkin_matrix(basis)
        new{xdim, vdim, typeof(spline), typeof(basis), typeof(coefficients),typeof(mass_matrix) }(spline, basis, coefficients, mass_matrix)
    end
end

Base.size(dist::SplineDistribution) = size(dist.coefficients)


function SplineDistribution(xdim, vdim, nknots::KT, order::OT, domain::Tuple, bc::Symbol=:Dirichlet) where {KT, OT}
    knots = collect(LinRange(domain..., nknots))
    b = BSplineBasis(BSplineOrder(order), knots)

    if bc == :Dirichlet
        basis = RecombinedBSplineBasis(Derivative(0), b)
    else #TODO: add other boundary condition options here
        basis = b
    end

    coefficients = zeros(length(basis))

    return SplineDistribution(xdim, vdim, basis, coefficients)
end