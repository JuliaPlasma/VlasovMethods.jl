struct SumMaxwellian{T} <: VlasovExample
    mean_one::T
    var_one::T
    mean_two::T
    var_two::T
    function SumMaxwellian(; mean_one = 0.0, var_one = 1.0, mean_two = 0.0, var_two = 2.0)
        new{typeof(mean_one)}(mean_one, var_one, mean_two, var_two)
    end
end

function initialize!(dist::ParticleDistribution, params::SumMaxwellian, ::SamplingMethod = NoSampling())
    # number of particles
    npart = length(dist.particles)
    vdim = size(dist.particles.v)[1]

    # random initial conditions
    v₀ = randn(vdim, Int(floor(npart / 2)))   # sample normal dist for v₀
    v₁ = randn(vdim, Int(ceil(npart / 2)))   # sample normal dist for v₀

    # adjust variance of the two distributions
    v₀ .*= sqrt(params.var_one)
    v₁ .*= sqrt(params.var_two)

    # adjust means of the two distributions
    v₀ .+= params.mean_one
    v₁ .+= params.mean_two

    v = hcat(v₀, v₁)
    # v = vcat(v₀, v₁)

    # write to particle distribution
    dist.particles.v .= v
    dist.particles.w[1, :] .= 1 / npart

    return dist
end
