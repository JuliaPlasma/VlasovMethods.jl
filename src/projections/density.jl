
@doc raw"""
Particle-sampled moments of the spline distribution function.

Each of these is a sum **over the particles** of a moment weighted by the spline evaluated at
the particle,

```math
\sum_\alpha \pi(v_\alpha) \, f_s(v_\alpha) ,
```

not the integral ``\int \pi(v) f(v) \, dv``. The two are different quantities and the
distinction matters: it is the particle-sampled form that appears in the coefficient system of
the conservative Lenard-Bernstein operator, where the same sum appears on both sides and the
sampling error cancels.
"""

# convenience function for computing the first three particle-sampled moments of f over v
function compute_f_densities(distribution::SplineDistribution, vp::AbstractArray{VT}) where {VT}
    n = projection_density(distribution, vp)
    μ = projection_momentum(distribution, vp)
    ε = projection_energy(distribution, vp)

    return n, μ, ε
end

function compute_df_densities(distribution::SplineDistribution, vp::AbstractArray{VT}) where {VT}
    n = projection_density(distribution, vp; isDerivative = true)
    μ = projection_momentum(distribution, vp; isDerivative = true)

    return n, μ
end

# Σ_α f_s(v_α)
function projection_density(distribution::SplineDistribution, vp::AbstractArray{VT}; kwargs...) where {VT}
    @inline f(v) = one(eltype(v))

    return density = projection(f, distribution, vp; kwargs...)
end

# Σ_α v_α f_s(v_α)
function projection_momentum(distribution::SplineDistribution, vp::AbstractArray{VT}; kwargs...) where {VT}
    @inline f(v) = v

    return momentum = projection(f, distribution, vp; kwargs...)
end

# Σ_α v_α² f_s(v_α)
function projection_energy(distribution::SplineDistribution, vp::AbstractArray{VT}; kwargs...) where {VT}
    @inline f(v) = v .^ 2

    return energy = projection(f, distribution, vp; kwargs...)
end

function projection(moment::Function, distribution::SplineDistribution,
        vp::AbstractArray{VT}; isDerivative::Bool = false) where {VT}
    # `derivative` shares the coefficient array with the spline, so this needs no rebuild and
    # sees whatever the last projection wrote. The earlier version tested
    # `typeof(spline) <: Spline` against a name that was never imported, so the derivative
    # branch raised `UndefVarError` rather than taking the derivative.
    fs = isDerivative ? derivative(distribution.spline) : distribution.spline

    return sum(moment.(vp) .* fs.(vp))
end
