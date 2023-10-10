using VlasovMethods
using BSplineKit

nknots = 20
domain = [0., 10.]
ord = 4
knots = collect(LinRange(domain..., nknots))
b = BSplineBasis(BSplineOrder(ord), knots)
basis = RecombinedBSplineBasis(Derivative(0), b)
coefs = zeros(nknots, nknots)

S = TwoDSpline(basis, coefs)

S([1.,1.])
