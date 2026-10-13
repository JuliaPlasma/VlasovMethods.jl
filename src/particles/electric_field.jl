"""
    ReducedElectricField(field, Ψₓ, Ψₚ, Ψₑ, X)
    ReducedElectricField(field, Ψ, X)

The reduced electric field of the full-order `field`: the lift, solve and project path, with no
hyper-reduction. The reduced coordinates `X` are lifted to the particle positions by `Ψₓ`, the
field is deposited from them and solved, evaluated at the positions lifted by `Ψₚ`, and projected
back by `Ψₑ`:

```julia
update!(f, X, w, t)     # deposit  Ψₓ X  and solve
efield!(f, E, X)        # E = Ψₑ E(Ψₚ X)
```

A `DEIMElectricField` is the same path with `Ψₚ = ΠₑᵀΨₓ` and `Ψₑ = ΨₓᵀΨₑ(ΠₑᵀΨₑ)⁻¹`, which lifts
the field evaluation to the interpolation points of the DEIM matrix instead of to every particle.

The three-argument form is the full lift and the transpose projection,
`ReducedElectricField(field, Ψ, X) = ReducedElectricField(field, Ψ, Ψ, Ψ', X)`, which is exact at
full basis size: the reduced field then equals the projection of the full-order field onto the
basis. See [`ElectricField`](@ref) for the interface, and [`DEIMElectricField`](@ref) for the
hyper-reduced path.
"""
struct ReducedElectricField{
    FT <: ElectricField, ΨₓT <: AbstractMatrix, ΨₚT <: AbstractMatrix,
    ΨₑT <: AbstractMatrix, AT <: AbstractArray} <: ElectricField
    field::FT
    Ψₓ::ΨₓT
    Ψₚ::ΨₚT
    Ψₑ::ΨₑT
    x::AT
    y::AT
    e::AT

    function ReducedElectricField(field, Ψₓ, Ψₚ, Ψₑ, X)
        x = zero(Ψₓ * X)
        y = zero(Ψₚ * X)
        e = zero(Ψₚ * X)
        new{typeof(field), typeof(Ψₓ), typeof(Ψₚ), typeof(Ψₑ), typeof(x)}(
            field, Ψₓ, Ψₚ, Ψₑ, x, y, e)
    end
end

ReducedElectricField(field, Ψ, X) = ReducedElectricField(field, Ψ, Ψ, Ψ', X)

function VlasovMethods.update!(f::ReducedElectricField, X::AbstractArray, w::AbstractArray, t)
    mul!(f.x, f.Ψₓ, X)
    update!(f.field, f.x, w, t)
end

function VlasovMethods.efield!(f::ReducedElectricField, E::AbstractArray, X::AbstractArray)
    mul!(f.y, f.Ψₚ, X)
    efield!(f.field, f.e, f.y)
    mul!(E, f.Ψₑ, f.e)
end

VlasovMethods.energy(f::ReducedElectricField) = energy(f.field)

VlasovMethods.coefficients(f::ReducedElectricField) = coefficients(f.field)

"""
    DEIMElectricField(field, Ψₓ, Ψₑ, Πₑ, X)

The hyper-reduced reduced electric field of the full-order `field`: a `ReducedElectricField` with
`Ψₚ = ΠₑᵀΨₓ` and `Ψₑ = ΨₓᵀΨₑ(ΠₑᵀΨₑ)⁻¹`, where `Πₑ` is the DEIM interpolation matrix of the
field basis `Ψₑ`. The lift of the field evaluation is

```math
E_{\\text{DEIM}} = \\Psi_x^T \\Psi_e (\\Pi_e^T \\Psi_e)^{-1} E(\\Pi_e^T \\Psi_x x),
```

so the field is evaluated at the ``k_e`` interpolation points of `Πₑ` rather than at every
particle, and the reduced acceleration agrees with the non-hyper-reduced one to the DEIM error
bound of the paper `ReducedBasisMethods.get_DEIM_interpolation_matrix` cites. `X` is the reduced
coordinate the buffers of the field are sized from.
"""
function DEIMElectricField(field, Ψₓ, Ψₑ, Πₑ, X)
    ΠₑᵀΨₓ = Πₑ' * Ψₓ
    ΨₓᵀPₑ = Ψₓ' * Ψₑ * inv(Πₑ' * Ψₑ)
    ReducedElectricField(field, Ψₓ, ΠₑᵀΨₓ, ΨₓᵀPₑ, X)
end
