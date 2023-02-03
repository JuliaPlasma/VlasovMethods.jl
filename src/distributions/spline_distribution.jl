
struct SplineDistribution{XD, VD, ST, BT, CT} <: DistributionFunction{XD,VD}
    spline::ST
    basis::BT
    coefficients::CT

    function SplineDistribution(xdim, vdim, basis::BT, coefficients::CT) where {BT, CT}
        spline = Spline(basis, coefficients)
        new{xdim, vdim, typeof(spline), BT, CT}(spline, basis, coefficients)
    end
end

Base.size(dist::SplineDistribution) = size(dist.coefficients)
