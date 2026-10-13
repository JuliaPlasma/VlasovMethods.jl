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
"""
abstract type ElectricField end

(f::ElectricField)(e, x, w, t) = (update!(f, x, w, t); efield!(f, e, x))
(f::ElectricField)(e, x) = efield!(f, e, x)
(f::ElectricField)(x, w, t) = (e = zero(x); f(e, x, w, t))

# The deposit of one vector of positions with one vector of weights, which the lifted positions
# of a reduced field are. The `_deposit!` of `src/projections/potential.jl` reads the positions
# and the weights from the first row of a state matrix, which a reduced field does not have: it
# lifts its reduced coordinates into a vector of positions alone. The loop is that one, indexed
# by the single axis.
function _deposit!(potential::PoissonSolvers.Potential, vals::AbstractVector,
        x::AbstractVector, w::AbstractVector)
    b = basis(potential)
    ρ = PoissonSolvers.rhs(potential)
    ρ .= 0

    for i in eachindex(x, w)
        j₀ = evaluate_all!(vals, b, x[i])
        for (t, value) in pairs(vals)
            # the padding of a recombined basis, as in `projection!`
            iszero(value) && continue
            ρ[basis_index(b, j₀ + t - 1)] += w[i] * value
        end
    end

    return potential
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
re-solve it in place. The time `t` is accepted for the `ElectricField` interface and ignored: the
field of a `PoissonField` is time-independent.
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
the potential's basis that the constructor tabulated.
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
