
struct ShiftedUniformDistribution{XT, VT, ST} <: VlasovExample
    xdomain::XT
    vdomain::VT
    shift::ST
    function ShiftedUniformDistribution(xdomain::XT = (0.0, 1.0), vdomain::VT = (-2.0, 2.0), shift::ST = 2.0) where {XT <: Tuple, VT <: Tuple, ST}
        new{XT, VT, ST}(xdomain, vdomain, shift)
    end
end

function initialize!(dist::ParticleDistribution, params::ShiftedUniformDistribution, ::SamplingMethod = NoSampling())
    # number of particles
    npart = length(dist.particles)

    # random initial conditions
    x₀ = rand(npart)    # sample uniform dist for x₀
    v₀ = rand(npart)   # sample normal dist for v₀

    # shift x₀ to domain
    x₀ .*= params.xdomain[end] - params.xdomain[begin]
    x₀ .+= params.xdomain[begin]

    # shift v₀ to domain
    v₀ .*= params.vdomain[end] - params.vdomain[begin]
    v₀ .+= params.vdomain[begin]

    # shift v₀ to the interval (domain[1] + shift, domain[2] + shift)
    v₀ .+= params.shift

    # write to particle distribution
    dist.particles.x[1,:] .= x₀
    dist.particles.v[1,:] .= v₀
    dist.particles.w[1,:] .= 1 / npart

    return dist

end
