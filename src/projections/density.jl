
# function projection!(out, xdensity, distribution) end
# function projection!(out, vdensity, distribution) end

# convenience function for computing first three moments of f over v
function compute_f_densities(distribution::SplineDistribution, particle_dist::ParticleDistribution)
    # n = zero(eltype(distribution))
    # μ = copy(n)
    # ε = copy(n)

    n = projection_density(distribution, particle_dist)
    μ = projection_momentum(distribution, particle_dist)
    ε = projection_energy(distribution, particle_dist)

    return n, μ, ε
end

function compute_df_densities(distribution::SplineDistribution, particle_dist::ParticleDistribution)
    n = projection_density(distribution, particle_dist; isDerivative=true)
    μ = projection_momentum(distribution, particle_dist; isDerivative=true)

    return n, μ
end

# compute density, i.e. n = ∫f dv
function projection_density(distribution::SplineDistribution, particle_dist::ParticleDistribution; kwargs...)
    @inline f(v) = one(eltype(v))

    return density = projection(f, distribution, particle_dist; kwargs...)
end

# compute mean momentum, i.e. nu = ∫vf dv
function projection_momentum(distribution::SplineDistribution, particle_dist::ParticleDistribution; kwargs...)
    @inline f(v) = v

    return momentum = projection(f, distribution, particle_dist; kwargs...)
end

# compute mean energy density, i.e. nε = ∫v²f dv
function projection_energy(distribution::SplineDistribution, particle_dist::ParticleDistribution; kwargs...)
    @inline f(v) = v.^2

    return energy = projection(f, distribution, particle_dist; kwargs...)
end

function projection(moment::Function, distribution::SplineDistribution, particle_dist::ParticleDistribution; isDerivative::Bool=false)
    if !(isDerivative)
        out = sum(moment.(particle_dist.particles.v) .* distribution.spline.(particle_dist.particles.v))
    elseif (isDerivative)
        df = Derivative(1) * distribution.spline
        out = sum(moment.(particle_dist.particles.v) .* df.(particle_dist.particles.v))
    end

    return out
end