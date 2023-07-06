struct ShiftedNormalV{DT,ST} <: VlasovExample
    domain::DT      # domain for x-dist
    shift::ST       # shift for v-dist
    function ShiftedNormalV(domain::DT = (-5.0, 5.0), shift::ST = 2.0) where {DT <: Tuple, ST}
        new{DT, ST}(domain, shift)
    end
end

function initialize!(dist::ParticleDistribution, params::ShiftedNormalV, ::SamplingMethod = NoSampling())
    # number of particles
    npart = length(dist.particles)

    # random initial conditions
    x₀ = rand(npart)    # sample uniform dist for x₀
    v₀ = randn(npart)   # sample normal dist for v₀

    # # shift x₀ to the interval [0,1]
    # xmax = ceil(maximum(abs.(x₀)))
    # x₀ .+= xmax
    # x₀ ./= 2*xmax

    # shift x₀ to domain
    x₀ .*= params.domain[end] - params.domain[begin]
    x₀ .+= params.domain[begin]

    # shift v₀ to the interval (domain[1] + shift, domain[2] + shift)
    v₀ .+= params.shift

    # write to particle distribution
    dist.particles.x[1,:] .= x₀
    dist.particles.v[1,:] .= v₀
    dist.particles.w[1,:] .= 1 / npart

    return dist

end
