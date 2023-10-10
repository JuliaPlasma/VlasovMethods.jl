struct DoubleMaxwellian{DT,ST} <: VlasovExample
    domain::DT      # domain for x-dist
    shift::ST       # shift away from 0 for both peaks of v-dist
    function DoubleMaxwellian(domain::DT = (-5.0, 5.0), shift::ST = 3.0) where {DT <: Tuple, ST}
        new{DT, ST}(domain, shift)
    end
end

function initialize!(dist::ParticleDistribution, params::DoubleMaxwellian, ::SamplingMethod = NoSampling())
    # number of particles
    npart = length(dist.particles)
    vdim = size(dist.particles.v)[1]

    # random initial conditions
    x₀ = rand(npart)    # sample uniform dist for x₀
    v₀ = randn(vdim, Int(floor(npart / 2)))   # sample normal dist for v₀
    v₁ = randn(vdim, Int(ceil(npart / 2)))   # sample normal dist for v₀

    # # shift x₀ to the interval [0,1]
    # xmax = ceil(maximum(abs.(x₀)))
    # x₀ .+= xmax
    # x₀ ./= 2*xmax

    # shift x₀ to domain
    x₀ .*= params.domain[end] - params.domain[begin]
    x₀ .+= params.domain[begin]

    # shift v₀ to the interval (domain[1] + shift, domain[2] + shift)
    v₀ .+= params.shift
    v₁ .-= params.shift
    v = hcat(v₀, v₁)
    # v = vcat(v₀, v₁)

    # write to particle distribution
    dist.particles.x[1,:] .= x₀
    dist.particles.v .= v
    dist.particles.w[1,:] .= 1 / npart

    return dist

end
