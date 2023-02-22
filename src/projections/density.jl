
# function projection!(out, xdensity, distribution) end
# function projection!(out, vdensity, distribution) end

# convenience function for computing first three moments of f over v
function compute_f_densities(distribution::SplineDistribution, vp::AbstractArray{VT}) where {VT}

    n = projection_density(distribution, vp)
    μ = projection_momentum(distribution, vp)
    ε = projection_energy(distribution, vp)

    return n, μ, ε
end

function compute_df_densities(distribution::SplineDistribution, vp::AbstractArray{VT}) where {VT}
    n = projection_density(distribution, vp; isDerivative=true)
    μ = projection_momentum(distribution, vp; isDerivative=true)

    return n, μ
end

# compute density, i.e. n = ∫f dv
function projection_density(distribution::SplineDistribution, vp::AbstractArray{VT}; kwargs...) where {VT}
    @inline f(v) = one(eltype(v))

    return density = projection(f, distribution, vp; kwargs...)
end

# compute mean momentum, i.e. nu = ∫vf dv
function projection_momentum(distribution::SplineDistribution, vp::AbstractArray{VT}; kwargs...) where {VT}
    @inline f(v) = v

    return momentum = projection(f, distribution, vp; kwargs...)
end

# compute mean energy density, i.e. nε = ∫v²f dv
function projection_energy(distribution::SplineDistribution, vp::AbstractArray{VT}; kwargs...) where {VT}
    @inline f(v) = v.^2

    return energy = projection(f, distribution, vp; kwargs...)
end

function projection(moment::Function, distribution::SplineDistribution, vp::AbstractArray{VT}; isDerivative::Bool=false) where {VT}
    if !(isDerivative)
        out = sum(moment.(vp) .* distribution.spline.(vp))
    elseif (isDerivative)
        df = Derivative(1) * distribution.spline
        out = sum(moment.(vp) .* df.(vp))
    end

    return out
end