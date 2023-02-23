
struct UniformDistribution{XT, VT} <: VlasovExample
    xdomain::XT
    vdomain::VT
    function UniformDistribution(xdomain::XT = (0.0, 1.0), vdomain::VT = (-2.0, 2.0)) where {XT <: Tuple, VT <: Tuple}
        new{XT, VT}(xdomain, vdomain)
    end
end


function initialize!(dist::ParticleDistribution, params::UniformDistribution, ::SamplingMethod = NoSampling())
    # number of particles
    npart = length(dist.particles)

    # random initial conditions
    x₀ = rand(npart)
    v₀ = rand(npart)

    # shift x₀ to domain
    x₀ .*= params.xdomain[end] - params.xdomain[begin]
    x₀ .+= params.xdomain[begin]

    # shift v₀ to domain
    v₀ .*= params.vdomain[end] - params.vdomain[begin]
    v₀ .+= params.vdomain[begin]

    # write to particle distribution
    dist.particles.x[1,:] .= x₀
    dist.particles.v[1,:] .= v₀
    dist.particles.w[1,:] .= 1 / npart

    return dist
end
