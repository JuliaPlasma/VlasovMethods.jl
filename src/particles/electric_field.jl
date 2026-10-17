"""
    ReducedElectricField(field, Ψₓ, Ψₚ, Ψₑ, X)
    ReducedElectricField(field, Ψ, X)

The reduced electric field of the full-order `field`: the lift, solve and project path, with no
hyper-reduction. The reduced coordinates `X` are lifted to the particle positions by `Ψₓ`, the
field is deposited from them and solved, evaluated at the positions lifted by `Ψₚ`, and projected
back by `Ψₑ`:

```julia
update!(f, X, w, t)     # deposit  Ψₓ X  and solve
efield!(f, e, X)        # e = Ψₑ E(Ψₚ X)
```

A `DEIMElectricField` of a field basis `Ψₑ` with DEIM row indices `p` is the same path, lifted by
`Ψₓ[p, :]` and projected by `ΨₓᵀΨₑ Ψₑ[p, :]⁻¹`, which evaluates the field at the `kₑ`
interpolation points `p` instead of at every particle.

For an orthonormal basis `Ψ`, the three-argument form is the full lift and the transpose projection,
`ReducedElectricField(field, Ψ, X) = ReducedElectricField(field, Ψ, Ψ, Ψ', X)`. At full basis size,
the reduced field equals the projection of the full-order field onto the basis. See
[`ElectricField`](@ref) for the interface, and [`DEIMElectricField`](@ref) for the
hyper-reduced path.
"""
struct ReducedElectricField{
    FT <: ElectricField, ΨₓT <: AbstractMatrix, ΨₚT <: AbstractMatrix,
    ΨₑT <: AbstractMatrix, AT <: AbstractVector} <: ElectricField
    field::FT
    Ψₓ::ΨₓT
    Ψₚ::ΨₚT
    Ψₑ::ΨₑT
    x::AT
    y::AT
    e::AT

    function ReducedElectricField(field, Ψₓ, Ψₚ, Ψₑ, X::AbstractVector)
        T = promote_type(eltype(Ψₓ), eltype(Ψₚ), eltype(X))
        x = similar(X, T, size(Ψₓ, 1))
        y = similar(X, T, size(Ψₚ, 1))
        e = similar(y)
        new{typeof(field), typeof(Ψₓ), typeof(Ψₚ), typeof(Ψₑ), typeof(x)}(
            field, Ψₓ, Ψₚ, Ψₑ, x, y, e)
    end
end

ReducedElectricField(field, Ψ, X) = ReducedElectricField(field, Ψ, Ψ, Ψ', X)

function update!(f::ReducedElectricField, X::AbstractVector, w::AbstractVector, t)
    mul!(f.x, f.Ψₓ, X)
    update!(f.field, f.x, w, t)
end

function efield!(f::ReducedElectricField, E::AbstractVector, X::AbstractVector)
    mul!(f.y, f.Ψₚ, X)
    efield!(f.field, f.e, f.y)
    mul!(E, f.Ψₑ, f.e)
end

energy(f::ReducedElectricField) = energy(f.field)

# The reduced field over the `ScaledField` of a march, at the scale `χ`: the same bases and a new
# `ScaledField` over the same inner field, so the buffers of the result are its own and the
# potential is the one the caller's field holds.
function _with_scale(f::ReducedElectricField{<:ScaledField}, χ)
    return ReducedElectricField(_with_scale(f.field, χ), f.Ψₓ, f.Ψₚ, f.Ψₑ, f.x)
end

coefficients(f::ReducedElectricField) = coefficients(f.field)

"""
    DEIMElectricField(field, Ψₓ, Ψₑ, p, X)
    DEIMElectricField(field, Ψₓ, Ψₑ, Πₑ, X)

The hyper-reduced electric field of the full-order `field`: a `ReducedElectricField` lifted by
`Ψₓ[p, :]` and projected by `ΨₓᵀΨₑ Ψₑ[p, :]⁻¹`, where `p` are the DEIM row indices of the field
basis `Ψₑ`, from `ReducedBasisMethods.get_DEIM_indices`. The second form takes the DEIM
interpolation matrix `Πₑ` of `get_DEIM_interpolation_matrix`, reads the row of the nonzero entry
of each column, and calls the first; it does not check that `Πₑ` is a selection matrix. The
buffers of the field follow the element type of `Ψₓ` and `X`, and `X` is the reduced coordinate
they are sized from. With `Πₑ = I[:, p]`, the lift of the field evaluation is

```math
E_{\\text{DEIM}} = \\Psi_x^T \\Psi_e (\\Pi_e^T \\Psi_e)^{-1} E(\\Pi_e^T \\Psi_x X),
```

so the field is evaluated at the ``k_e`` interpolation points rather than at every particle.
For orthonormal position and field bases, the reduced acceleration agrees with the
non-hyper-reduced one to the DEIM error bound of the paper
`ReducedBasisMethods.get_DEIM_interpolation_matrix` cites.
"""
function DEIMElectricField(field, Ψₓ, Ψₑ, p::AbstractVector{<:Integer}, X)
    ΨₓᵀPₑ = (Ψₓ' * Ψₑ) / Ψₑ[p, :]
    ReducedElectricField(field, Ψₓ, Ψₓ[p, :], ΨₓᵀPₑ, X)
end

function DEIMElectricField(field, Ψₓ, Ψₑ, Πₑ::AbstractMatrix, X)
    DEIMElectricField(field, Ψₓ, Ψₑ, [findfirst(!iszero, c) for c in eachcol(Πₑ)], X)
end
