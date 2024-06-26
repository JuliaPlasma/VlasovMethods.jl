
struct NormalDistribution{DT} <: VlasovExample
    domain::DT
    function NormalDistribution(domain::DT = (0.0, 1.0)) where {DT <: Tuple}
        new{DT}(domain)
    end
end


function initialize!(dist::ParticleDistribution, params::NormalDistribution, ::SamplingMethod = NoSampling())
    # number of particles
    npart = length(dist.particles)
    vsize = size(dist.particles.v)

    # random initial conditions
    x₀ = randn(npart)
    # v₀ = randn(npart)
    Random.seed!(1234)
    
    v₀ = randn(vsize)


    # shift x₀ to the interval [0,1]
    xmax = ceil(maximum(abs.(x₀)))
    x₀ .+= xmax
    x₀ ./= 2*xmax

    # shift x₀ to domain
    x₀ .*= params.domain[end] - params.domain[begin]
    x₀ .+= params.domain[begin]

    # concatenate x₀ and v₀
    # z₀ = collect(hcat(x₀, v₀)')

    # write to particle distribution
    dist.particles.x[1,:] .= x₀
    dist.particles.v .= v₀
    # dist.particles.v[1,:] .= v₀
    dist.particles.w[1,:] .= 1 / npart

    return dist
end
