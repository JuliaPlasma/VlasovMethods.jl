@doc raw"""
    SplineDistribution{DT, XD, VD}

The distribution function represented as a spline in velocity,
``f_s(v) = \sum_I \hat{f}_I \, \varphi_I(v)``, together with everything needed to project onto
it and to integrate against it.

`VD` is the number of velocity dimensions. For `VD == 1` the basis is a one-dimensional
`SimpleSplines` B-spline basis and the coefficients are a vector; for `VD > 1` it is a
`TensorProductBasis` and the coefficients are a `VD`-dimensional array. Nothing else about
the two cases differs — the evaluation, the projection and the mass solve are the same calls.

# Fields

  - `spline` — the `SimpleSplines.Spline`, which **shares** `coefficients`. A projection
    writes that array, and every reference to the spline sees the new function without a
    rebuild.
  - `basis` — the spline basis.
  - `quadrature` — the assembly table: the basis tabulated at the Gauß-Legendre points, the
    weights, and the mass operator.
  - `coefficients` — the degrees of freedom ``\hat{f}``.

# The mass matrix is never formed as a Kronecker product

For `VD > 1` the mass matrix is ``\mathbb{M}^{(VD)} \otimes \dots \otimes \mathbb{M}^{(1)}``,
and a solve is `VD` one-dimensional solves applied along each axis. That is an identity, not
an approximation. The earlier implementation built `kron(M, M)` and took a dense `cholesky`
of it — 1681×1681 for the 41-element cubic basis of a two-dimensional velocity space, and
68921² in three dimensions, which does not fit.

# The boundary condition is not a free choice

The ``L^2`` projection reproduces ``\int \pi(v) f \, dv`` for a polynomial ``\pi`` exactly when
``\pi`` lies in the span of the basis, so the mass, momentum and energy of ``f_s`` agree with
those of the particles only if ``1``, ``v`` and ``v^2`` are in the span.
`polynomial_reproduction(basis)` reports the largest degree available: `p` for the clamped
basis (`Free()`), `0` for a periodic one — ``v`` is not periodic — and `-1` for a
Dirichlet-recombined one, where not even the constants survive.

Measured on a cubic basis over ``[-8,8]`` by `scripts/verify_conservation.jl`: `Free()`
reproduces all three to ``10^{-14}``, `Periodic()` gets ``v`` wrong by ``1.2`` and ``v^2`` by
``0.12``, and `Dirichlet()` gets even the constant wrong by ``4 \times 10^{-2}``.

!!! note "What this does and does not break"
    It does **not** by itself break the particle-level conservation of the conservative
    Lenard-Bernstein operator: its coefficients ``A_1, A_2`` are solved from the two
    conservation constraints, which are conditions on particle sums, so
    ``\sum_\alpha w_\alpha \dot{v}_\alpha`` and
    ``\sum_\alpha w_\alpha v_\alpha \dot{v}_\alpha`` vanish at round-off on any basis. What it
    breaks is the agreement between the particle and spline representations, and therefore
    every diagnostic computed from ``f_s`` — the entropy above all, which is what the
    H-theorem is a statement about. See [`check_conservation_basis`](@ref).
"""
struct SplineDistribution{DT, XD, VD, ST, BT, QT, CT} <: DistributionFunction{DT, XD, VD}
    spline::ST
    basis::BT
    quadrature::QT
    coefficients::CT

    function SplineDistribution{XD, VD}(basis::BT, quadrature::QT,
            coefficients::CT) where {XD, VD, BT, QT, CT <: AbstractArray}
        ndims(coefficients) == VD || throw(DimensionMismatch(
            "a $(VD)-dimensional velocity space needs a $(VD)-dimensional coefficient " *
            "array, got $(ndims(coefficients))"))
        spline = Spline(basis, coefficients)
        DT = eltype(coefficients)
        new{DT, XD, VD, typeof(spline), BT, QT, CT}(spline, basis, quadrature, coefficients)
    end
end

function SplineDistribution(xdim::Integer, vdim::Integer, basis, quadrature,
        coefficients::AbstractArray)
    SplineDistribution{Int(xdim), Int(vdim)}(basis, quadrature, coefficients)
end

Base.eltype(::SplineDistribution{DT}) where {DT} = DT
Base.length(dist::SplineDistribution) = length(dist.coefficients)
Base.size(dist::SplineDistribution) = size(dist.coefficients)
Base.ndims(::SplineDistribution{DT, XD, VD}) where {DT, XD, VD} = VD

basis(dist::SplineDistribution) = dist.basis
quadrature(dist::SplineDistribution) = dist.quadrature
coefficients(dist::SplineDistribution) = dist.coefficients

"""
    mass_operator(dist::SplineDistribution)

The mass operator of the basis — a factorised or circulant one for `VD == 1`, a
`KroneckerMass` for `VD > 1`. Solve with `mass_solve!(out, mass_operator(dist), rhs)`.
"""
SimpleSplines.mass_operator(dist::SplineDistribution) = mass_operator(dist.quadrature)

"""
    mass_matrix(dist::SplineDistribution)

The assembled mass matrix. For `VD > 1` this **forms** the Kronecker product, which the
representation exists to avoid; it is here for tests and diagnostics, not for solving.
"""
SimpleSplines.mass_matrix(dist::SplineDistribution) = mass_matrix(dist.quadrature)

function Base.similar(AT, s::SplineDistribution{DT, XD, VD}) where {DT, XD, VD}
    SplineDistribution{XD, VD}(s.basis, s.quadrature, zeros(AT, size(s.coefficients)))
end

# The concrete type of `similar(AT, s)`, computed from type parameters alone so that the cache
# lookup in `CacheDict` stays inferable. `Spline`'s own element type is the promotion of the
# basis's with the coefficients', which is what `promote_type` reproduces here.
function similar_type(AT, ::SplineDistribution{
        DT, XD, VD, ST, BT, QT, CT}) where {
        DT, XD, VD, ST, BT, QT, CT}
    NCT = Array{AT, VD}
    NST = Spline{promote_type(eltype(BT), AT), BT, NCT}
    SplineDistribution{AT, XD, VD, NST, BT, QT, NCT}
end

function Base.show(io::IO, dist::SplineDistribution{DT, XD, VD}) where {DT, XD, VD}
    print(io, "SplineDistribution{", DT, ", xdim=", XD, ", vdim=", VD, "}(",
        size(dist.coefficients), ")")
end

@doc raw"""
    check_conservation_basis(dist::SplineDistribution; moments = 2)

Throw unless the basis reproduces every polynomial of degree `≤ moments`, i.e. unless the
moments ``\int v^k f_s \, dv`` for `k ≤ moments` agree with the particles'.

Not called automatically: a diagnostic run may legitimately want a Dirichlet or periodic basis,
and refusing to build one would be wrong. Call it where a result depends on the two
representations agreeing — before reporting an entropy, or before comparing a spline moment
against a particle one.
"""
function check_conservation_basis(dist::SplineDistribution; moments::Integer = 2)
    m = polynomial_reproduction(dist.basis)
    m ≥ moments && return dist
    throw(ArgumentError(
        "this basis reproduces polynomials only up to degree $(m), so ∫ vᵏ f_s dv does not " *
        "equal Σ_α w_α v_α^k for k ≤ $(moments) and any diagnostic read off f_s describes a " *
        "distribution with the wrong moments. The boundary condition is the cause: `Free()` " *
        "gives degree p, `Periodic()` gives 0 (v is not periodic) and `Dirichlet()` gives -1 " *
        "(not even the constants survive)."))
end

# The boundary conditions the previous implementation named by symbol. `:nothing` selected the
# unconstrained basis by falling through an `else` branch rather than by saying so, and a
# mistyped symbol did the same silently; the mapping is written out here so that an
# unrecognised name reaches `SimpleSplines.BoundaryCondition` and throws.
_boundary_condition(bc::SimpleSplines.BoundaryCondition) = bc
_boundary_condition(bc::Tuple) = bc
function _boundary_condition(bc::Symbol)
    bc === :Dirichlet && return Dirichlet()
    bc === :Periodic && return Periodic()
    bc === :Natural && return Free()
    bc === :nothing && return Free()
    return SimpleSplines.BoundaryCondition(bc)
end

@doc raw"""
    SplineDistribution(xdim, vdim, nknots, s_order, domain, length_big_cell, bc = Free())

Build a `SplineDistribution` on `nknots` equally spaced breakpoints across `domain`, of spline
order `s_order`, with boundary condition `bc` on every velocity axis.

`s_order` is the **order** ``k = p+1``, the same convention as `SimpleSplines.order` — so
`s_order = 4` is the cubic basis both collision-operator manuscripts specify, and
`s_order = 3` is quadratic.

`length_big_cell > 0` appends one oversized cell of that width at each end of `domain`,
outside it, giving a `GeneralMesh`. That is the device for keeping a particle that strays out
of the resolved region inside the support of the basis; with `0` the mesh is uniform on
`domain` exactly.

`bc` may be a `SimpleSplines.BoundaryCondition`, a lowercase symbol, or a two-tuple for the
two ends. The default is `Free()` — the unconstrained clamped basis, which is the only choice
that reproduces ``1``, ``v`` and ``v^2`` and hence the only one on which the conservative
schemes conserve. The previous default was `:Dirichlet`, which reproduces nothing at all.

!!! note "The lumped-mass option is gone"
    The old signature took a trailing `compute_mass_galerkin::Bool`. Passing `false` assembled
    the mass matrix with a trapezoidal rule, which for a B-spline basis evaluates at the knots
    where ``\varphi_i(t_j) = \delta_{ij}`` and therefore produced a **diagonal**, mass-lumped
    matrix rather than ``\mathbb{M}_{ij} = \int \varphi_i \varphi_j``. That is a different
    discretisation, not a cheaper assembly of the same one, and it is not the one either
    manuscript describes. The mass matrix is now always the exact Galerkin one. A call passing
    the old flag is a `MethodError` rather than being silently reinterpreted.
"""
function SplineDistribution(xdim::Integer, vdim::Integer, nknots::Integer,
        s_order::Integer, domain, length_big_cell = 0, bc = Free())
    a, b = first(domain), last(domain)
    s_order ≥ 1 || throw(ArgumentError(
        "the spline order must be at least 1, got s_order = $(s_order)"))
    nknots ≥ 2 || throw(ArgumentError(
        "at least two breakpoints are needed, got nknots = $(nknots)"))

    p = s_order - 1

    ts = collect(LinRange(a, b, nknots))
    mesh = if length_big_cell > 0
        GeneralMesh([a - length_big_cell; ts; b + length_big_cell])
    else
        UniformMesh(nknots - 1, a .. b)
    end

    condition = _boundary_condition(bc)
    axis = BSplineBasis(mesh, p, condition)

    if vdim == 1
        q = SplineQuadrature(axis)
        coefficients = zeros(Float64, nbasis(axis))
    else
        # One independent basis per velocity axis. They are built from the same mesh, degree
        # and condition here because that is what this constructor's arguments describe; a
        # basis differing per axis is made by calling `TensorProductBasis` directly and
        # passing it to the other constructor.
        B = TensorProductBasis(ntuple(_ -> BSplineBasis(mesh, p, condition), vdim))
        q = TensorProductQuadrature(B)
        axis = B
        coefficients = zeros(Float64, size(B)...)
    end

    return SplineDistribution{Int(xdim), Int(vdim)}(axis, q, coefficients)
end

## Evaluation

function evaluate(dist::SplineDistribution{DT, XD, 1}, ::Any, v) where {DT, XD}
    dist.spline(v isa Number ? v : v[1])
end

function evaluate(dist::SplineDistribution{DT, XD, VD}, ::Any, v) where {DT, XD, VD}
    dist.spline(v)
end
