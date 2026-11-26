@doc raw"""
    TensorSplineDistribution(xbasis, vbasis, coefficients::AbstractMatrix)
    TensorSplineDistribution(xbasis, vbasis, f)

The distribution function of a 1D1V phase space as a tensor product spline,

```math
f(x, v) = \sum_{i, j} A_{ij} \, B_i(x) \, C_j(v) ,
```

with the periodic B-spline basis ``B_i`` of `xbasis` in ``x``, the B-spline basis ``C_j`` of
`vbasis` in ``v``, and the `nbasis(xbasis) × nbasis(vbasis)` coefficient matrix ``A``. It is a
`DistributionFunction{T, 1, 1}`, with `T` the element type of the two bases. The layout of ``A``
is the one of the values of a [`GridDistribution`](@ref): `i` runs over ``x``, `j` over ``v``,
and the column-major order runs ``x`` fastest.

The second form holds the ``L^2`` projection of the function `f(x, v)`,
``A = \mathbb{M}_x^{-1} F \mathbb{M}_v^{-1}`` with ``F_{ij} = \int\!\!\int f B_i C_j \, dx \, dv``
by the tensor quadrature, in the element type `T`.

The fields are `xbasis` and `vbasis`, their default `SplineQuadrature`s `xquadrature` and
`vquadrature`, which hold the basis tables and the mass operators, and `coefficients`, the matrix
``A``, which the first form holds without a copy. `evaluate(dist, x, v)` returns ``f(x, v)``.

The constructors throw:

  - `ArgumentError` unless `xbasis` is a `PeriodicBSplineBasis`;
  - `ArgumentError` unless `vbasis` reproduces ``v^2``, `polynomial_reproduction(vbasis) ≥ 2`,
    so that its span holds the kinetic energy ``v^2 / 2`` of the Vlasov–Poisson model: a
    Dirichlet, a periodic and a linear clamped basis do not;
  - `ArgumentError` unless the two bases and the coefficients have one element type;
  - `DimensionMismatch` unless `size(coefficients) == (nbasis(xbasis), nbasis(vbasis))`.

A clamped (`Free()`) cubic basis reproduces ``1``, ``v`` and ``v^2``, and on a velocity domain
symmetric about 0 it is symmetric under the reflection ``v \mapsto -v``.
"""
struct TensorSplineDistribution{DT, XB <: PeriodicBSplineBasis{DT},
    VB <: AbstractBSplineBasis{DT}, XQ <: SplineQuadrature{DT}, VQ <: SplineQuadrature{DT},
    AT <: AbstractMatrix{DT}} <: DistributionFunction{DT, 1, 1}
    xbasis::XB
    vbasis::VB
    xquadrature::XQ
    vquadrature::VQ
    coefficients::AT

    function TensorSplineDistribution(xbasis::AbstractBSplineBasis,
            vbasis::AbstractBSplineBasis, coefficients::AbstractMatrix)
        _check_tensor_bases(xbasis, vbasis, eltype(coefficients))
        expected = (nbasis(xbasis), nbasis(vbasis))
        size(coefficients) == expected || throw(DimensionMismatch(
            "the coefficients of a TensorSplineDistribution have the size $expected of the " *
            "bases; got $(size(coefficients))"))
        xq, vq = SplineQuadrature(xbasis), SplineQuadrature(vbasis)
        new{eltype(xbasis), typeof(xbasis), typeof(vbasis), typeof(xq), typeof(vq),
            typeof(coefficients)}(xbasis, vbasis, xq, vq, coefficients)
    end
end

function TensorSplineDistribution(xbasis::AbstractBSplineBasis, vbasis::AbstractBSplineBasis, f)
    T = eltype(xbasis)
    dist = TensorSplineDistribution(xbasis, vbasis, zeros(T, nbasis(xbasis), nbasis(vbasis)))
    l2_projection!(dist.coefficients, TensorProductQuadrature(xbasis ⊗ vbasis), x -> f(x...))
    return dist
end

function _check_tensor_bases(xbasis, vbasis, T)
    xbasis isa PeriodicBSplineBasis || throw(ArgumentError(
        "the x-basis of a TensorSplineDistribution is periodic, a PeriodicBSplineBasis; " *
        "got a $(nameof(typeof(xbasis)))"))
    polynomial_reproduction(vbasis) ≥ 2 || throw(ArgumentError(
        "the v-basis of a TensorSplineDistribution reproduces v², so that the kinetic energy " *
        "v²/2 is in its span: a clamped (Free()) basis of degree 2 or more; got a basis that " *
        "reproduces polynomials of degree $(polynomial_reproduction(vbasis))"))
    eltype(xbasis) === eltype(vbasis) === T || throw(ArgumentError(
        "the bases and the coefficients of a TensorSplineDistribution have one element type; " *
        "got $(eltype(xbasis)), $(eltype(vbasis)) and $T"))
    return nothing
end

# Mx⁻¹ R Mv⁻¹, by the Kronecker product of the mass operators of the two quadratures
function _tensor_mass_solve(xq::SplineQuadrature, vq::SplineQuadrature, R::AbstractMatrix)
    KroneckerMass(mass_operator(xq), mass_operator(vq)) \ R
end

Base.eltype(::TensorSplineDistribution{DT}) where {DT} = DT
Base.size(dist::TensorSplineDistribution) = size(dist.coefficients)

function Base.show(io::IO, dist::TensorSplineDistribution{DT}) where {DT}
    print(io, "TensorSplineDistribution{", DT, "}(", size(dist.coefficients), ")")
end

function evaluate(dist::TensorSplineDistribution, x, v)
    # one number type for the point, whose element type sets the buffers of the evaluation
    evaluate(dist.xbasis ⊗ dist.vbasis, dist.coefficients, promote(only(x), only(v)))
end
