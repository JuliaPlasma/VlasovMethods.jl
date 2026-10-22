"""
    ElectricField

The electric field of a one-dimensional Vlasov–Poisson problem on particles.

A field deposits particle positions with their weights, solves for the potential, and evaluates
the field ``E(x) = -\\partial_x \\phi(x)`` at particle positions. Its interface is two functions
and three call forms:

```julia
update!(f, x, w, t)     # deposit the positions `x` with the weights `w` at the time `t`
efield!(f, e, x)        # write `E` at each position of `x` into `e`
f(e, x, w, t)           # `update!` followed by `efield!`
f(e, x)                 # `efield!` alone
f(x, w, t)              # both, into an `e` that call allocates
```

`update!` and `efield!` are called once per time step of a march, so a concrete field holds the
buffers it needs and neither function allocates. `energy(f)` and `coefficients(f)` read the field
energy and the potential coefficients of `f`. A reduced field, such as `ReducedElectricField`,
holds a full-order field and lifts into it.

`update!` and `efield!` are not exported, so a caller reaches them through the three call forms
above, as `VlasovMethods.update!` and `VlasovMethods.efield!`, or with
`using VlasovMethods: update!, efield!`. Neither is the `PoissonSolvers.update!` that a
[`PoissonField`](@ref) calls on its own potential: that one solves the potential's linear system,
and it is a different function that happens to carry the same name.
"""
abstract type ElectricField end

(f::ElectricField)(e, x, w, t) = (update!(f, x, w, t); efield!(f, e, x))
(f::ElectricField)(e, x) = efield!(f, e, x)
(f::ElectricField)(x, w, t) = (e = zero(x); f(e, x, w, t))

# The deposit of one vector of positions with one vector of weights, which the lifted positions
# of a reduced field are. It is the loop of `src/projections/potential.jl`, reached through the
# row view a `transpose` of a vector is: that view wraps the caller's buffers, so the deposit
# copies nothing, and one loop serves both a state matrix and a lifted vector.
function _deposit!(potential::PoissonSolvers.Potential, vals::AbstractVector,
        x::AbstractVector, w::AbstractVector)
    return _deposit!(potential, vals, transpose(x), transpose(w))
end

"""
    PoissonField(potential)

The full-order electric field of the PoissonSolvers `Potential` `potential`, the field a particle
`VlasovPoisson` model integrates.

The deposit of [`update!`](@ref) is `_deposit!`, the field of [`efield!`](@ref) is
`_electric_field`, and the solve is the potential's own `PoissonSolvers.update!`.

`potential` is on a B-spline basis, as the potential of a particle model is. The constructor holds
the buffer of the `local_width(basis)` basis values that the deposit and the field evaluation
fill, and the quadrature of the basis that [`energy`](@ref) integrates ``(\\partial_x\\phi)^2``
with, so that neither allocates on a call.
"""
struct PoissonField{
    PT <: PoissonSolvers.Potential, VT <: AbstractVector, XT <: AbstractVector,
    WT <: AbstractVector} <: ElectricField
    potential::PT
    work::VT
    nodes::XT
    weights::WT
end

function PoissonField(potential::PoissonSolvers.Potential)
    quadrature = SplineQuadrature(basis(potential))
    return PoissonField(
        potential,
        _deposit_buffer(potential),
        quadrature_nodes(quadrature),
        quadrature_weights(quadrature))
end

"""
    update!(f::PoissonField, x, w, t)

Deposit the positions `x` with the weights `w` into the right-hand side of `f.potential` and
re-solve it in place. The time `t` is accepted for the `ElectricField` interface and ignored:
the field has no explicit time dependence.
"""
function update!(f::PoissonField, x::AbstractVector, w::AbstractVector, t)
    _deposit!(f.potential, f.work, x, w)
    PoissonSolvers.update!(f.potential)
    return nothing
end

"""
    efield!(f::PoissonField, e, x)

Write the field ``E(x_i) = -\\partial_x\\phi(x_i)`` at each position of `x` into `e`, and return
`e`.
"""
function efield!(f::PoissonField, e::AbstractVector, x::AbstractVector)
    for i in eachindex(e, x)
        e[i] = _electric_field(f.potential, f.work, x[i])
    end
    return e
end

"""
    energy(f::PoissonField)

The field energy ``\\frac{1}{2}\\int (\\partial_x\\phi)^2 \\,\\mathrm{d}x``, by the quadrature of
the potential's basis that the constructor tabulates. It is the diagnostic a snapshot march saves
at every state, so it is type-stable and allocates nothing.
"""
function energy(f::PoissonField)
    e² = zero(eltype(f.weights))
    for (x, ω) in zip(f.nodes, f.weights)
        ∂ϕ = _electric_field(f.potential, f.work, x)
        e² += ω * ∂ϕ^2
    end
    return e² / 2
end

coefficients(f::PoissonField) = PoissonSolvers.coefficients(f.potential)

"""
    ScaledField(field::ElectricField, χ::Real)

The full-order `field` as the rescaled state of a march sees it: `efield!` divides the field of
`field` by `χ²`, and [`energy`](@ref) divides its energy by the same `χ²`. That is not the energy
of the divided field, which is quadratic in the field and would carry `χ⁴`. [`update!`](@ref) and
`coefficients` forward to `field` unchanged, because the deposit and the solve carry no amplitude.

It is an `ElectricField`, so a reduced march wraps it as
`ReducedElectricField(ScaledField(field, χ), Ψ, X)`. The `χ` stored here must be the `params.χ` of
the march the field is handed to: the march takes its step scale and the `χ²` of the field from
`params.χ` alone, so a `ScaledField`, or a reduced field over one, at another `χ` is refused with
an `ArgumentError` that names both values rather than rescaled in silence (see
[`integrate_vp`](@ref)). The stored `χ` therefore serves a call outside a march, and the
construction of one at the scale the march will take.

`χ` must be finite and positive; the constructor throws an `ArgumentError` otherwise, which
refuses `0`, a negative and `NaN` alike. Each division converts `χ` to the element type of the
value it divides, so a `Float32` field with a `Float64` `χ` keeps `Float32`.
"""
struct ScaledField{FT <: ElectricField, T <: Real} <: ElectricField
    field::FT
    χ::T

    function ScaledField(field::FT, χ::T) where {FT <: ElectricField, T <: Real}
        isfinite(χ) && χ > 0 || throw(ArgumentError(
            "χ must be finite and positive, got χ = $χ"))
        return new{FT, T}(field, χ)
    end
end

function update!(f::ScaledField, x::AbstractVector, w::AbstractVector, t)
    return update!(f.field, x, w, t)
end

function efield!(f::ScaledField, e::AbstractVector, x::AbstractVector)
    efield!(f.field, e, x)
    χ = convert(eltype(e), f.χ)
    e ./= χ^2
    return e
end

function energy(f::ScaledField)
    E = energy(f.field)
    return E / convert(typeof(E), f.χ)^2
end

coefficients(f::ScaledField) = coefficients(f.field)

# The field of a march at the scale `χ` of its parameters, which are the one source of the scale: a
# `ScaledField` is rebuilt over the same inner field with `χ`, and any other field is its own. A
# stored `χ` that the rebuild would discard is refused, because the march's step and its field would
# then carry two scales. `src/particles/electric_field.jl` gives a reduced field over a
# `ScaledField` the same rule.
_with_scale(f::ElectricField, χ) = f
function _with_scale(f::ScaledField, χ)
    f.χ == χ || throw(ArgumentError(
        "the field is a ScaledField at χ = $(f.χ); a march takes its step and the scale of its " *
        "field from params.χ = $χ, so the two would differ"))
    return ScaledField(f.field, χ)
end
