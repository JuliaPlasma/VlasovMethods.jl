
struct ReducedElectricField{
    FT <: ElectricField, ΨₓT <: AbstractMatrix, ΨₚT <: AbstractMatrix,
    ΨₑT <: AbstractMatrix, AT <: AbstractArray} <: ElectricField
    field::FT
    Ψₓ::ΨₓT
    Ψₚ::ΨₑT
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

# $ \Psi_x^T \Psi_{DEIM} \nabla \phi ( \Pi_{DEIM} \Psi_x x) $
function DEIMElectricField(field, Ψₓ, Ψₑ, Πₑ, X)
    ΠₑᵀΨₓ = Πₑ' * Ψₓ
    ΨₓᵀPₑ = Ψₓ' * Ψₑ * inv(Πₑ' * Ψₑ)
    ReducedElectricField(field, Ψₓ, ΠₑᵀΨₓ, ΨₓᵀPₑ, X)
end
