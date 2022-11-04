
abstract type ElectricField end

function (f::ElectricField)(e::AbstractArray, x::AbstractArray, w::AbstractArray, t)
    update!(f, x, w, t)
    efield!(f, e, x)
end

function (f::ElectricField)(e::AbstractArray, x::AbstractArray)
    efield!(f, e, x)
end

function (f::ElectricField)(x::AbstractArray, w::AbstractArray, t)
    e = zero(x)
    f(e, x, w, t)
    return e
end



struct ScaledField{FT <: ElectricField, χT} <: ElectricField
    field::FT
    χ::χT
end

function efield!(f::ScaledField, e::AbstractArray, x::AbstractArray)
    efield!(f.field, e, x)
    e ./= f.χ^2
end

update!(f::ScaledField, x::AbstractArray, w::AbstractArray, t) = update!(f.field, x, w, t)

energy(f::ScaledField) = energy(f.field) / f.χ^2

coefficients(f::ScaledField) = coefficients(f.field)



struct PoissonField{PT <: PoissonSolver} <: ElectricField
    poisson::PT
end

efield!(f::PoissonField, e::AbstractArray, x::AbstractArray) = eval_field!(e, f.poisson, x)

update!(f::PoissonField, x::AbstractArray, w::AbstractArray, t) = solve!(f.poisson, x, w)

energy(f::PoissonField) = dot(f.poisson.ϕ, f.poisson.S, f.poisson.ϕ) / 2

coefficients(f::PoissonField) = f.poisson.ϕ

ScaledPoissonField(poisson::PoissonSolver, χ) = ScaledField(PoissonField(poisson), χ)



mutable struct ExternalField{PT <: PoissonSolver, DT, TT, AT <: AbstractArray{DT}} <: ElectricField
    poisson::PT
    coeffs::OffsetMatrix{DT,AT}
    Δt::TT
    ts::Int

    ExternalField(poisson::PoissonSolver, coeffs::OffsetMatrix{DT,AT}, Δt::TT) where {DT,TT,AT} = new{typeof(poisson), DT, TT, AT}(poisson, coeffs, Δt, 0)
end

ExternalField(poisson::PoissonSolver, coeffs::AbstractMatrix, Δt) = ExternalField(poisson, OffsetArray(coeffs, axes(coeffs,1), axes(coeffs,2) .- 1), Δt)

function update!(f::ExternalField, x::AbstractArray, w::AbstractArray, t)
    f.ts = round(t / f.Δt)
    f.poisson.ϕ .= f.coeffs[:, f.ts]
end

efield!(f::ExternalField, e::AbstractArray, x::AbstractArray) = eval_field!(e, f.poisson, x)

coefficients(f::ExternalField) = f.poisson.ϕ

energy(f::ExternalField) = dot(f.poisson.ϕ, f.poisson.S, f.poisson.ϕ) / 2

ScaledExternalField(poisson::PoissonSolver, coeffs, Δt, χ) = ScaledField(ExternalField(poisson, coeffs, Δt), χ)
