
struct BumpOnTail{T} <: VlasovExample
    ε::T    # amplitude of perturbation
    κ::T    # wave number of perturbation
    α::T    # fraction of particles in bump
    σ::T    # standard deviation of bump
    v₀::T   # offset velocity of bump

    function BumpOnTail(;
            ε::T = 0.03, κ::T = 0.3, α::T = 0.1, σ::T = 0.5, v₀::T = 4.5) where {T}
        new{T}(ε, κ, α, σ, v₀)
    end
end

# x-part of distribution function
function fₓ(x::T, params::BumpOnTail) where {T}
    @unpack κ, ε = params
    @. (1 - ε * cos(κ*x))
end

# v-part of distribution function
function fᵥ(v::T, params::BumpOnTail) where {T}
    @unpack α, v₀, σ = params
    @. ((1-α) * exp(- v^2 / 2) + α / σ * exp(- (v-v₀)^2 / (2σ^2))) / √(2π)
end

# distribution function
function f(x::T, v::T, params::BumpOnTail) where {T}
    fₓ(x, params) .* fᵥ(v, params)
end

function initialize!(dist::ParticleDistribution, params::BumpOnTail,
        sampling::SamplingMethod = AcceptRejectSampling())
    draw!(dist, fₓ, params, sampling)
end

"""
    draw!(dist, fₓ, params::BumpOnTail, ::AcceptRejectSampling)

Fill the 1D1V distribution `dist` with the bump-on-tail particles of `params`, and return it. Each
position is a Sobol proposal, uniform on `[0, 2π/κ]`, accepted with the probability
`fₓ(x)/(1 + ε)`, and each particle has the weight `2π/κ/N`. Each velocity is drawn from the bulk
`N(0, 1)` by the inverse CDF, and moved into the tail, `σv + v₀`, with the probability `α`.
"""
function draw!(dist::ParticleDistribution{<:Any, 1, 1}, fₓ::Base.Callable,
        params::BumpOnTail, ::AcceptRejectSampling)
    @unpack ε, κ, α, σ, v₀ = params

    # number of particles
    N = length(dist.particles)

    # Sobol sampling
    s = SobolSeq(2)
    skip(s, 2N)

    n = 1
    while n < N+1
        ###
        y = Sobol.next!(s)
        x₀ = y[1]*2π/κ  # proposal, uniform

        # accept-reject in x
        if rand(1)[1] > ((fₓ(x₀, params))/(1 + ε))
            continue                            # reject
        end
        dist.particles.x[1, n] = x₀                               # accept
        dist.particles.w[1, n] = 2π / κ / N

        # inverse CDF sampling in v
        vₛ = 2y[2] - 1
        dist.particles.v[1, n] = √2 * erfinv(vₛ)  # draw from bulk
        if rand(1)[1] > 1-α
            dist.particles.v[1, n] = dist.particles.v[1, n] * σ + v₀  # draw from tail
        end
        n = n+1
    end

    return dist
end

"""
    draw!(dist, fₓ, params::BumpOnTail, ::ImportanceSampling)

Fill the 1D1V distribution `dist` with the bump-on-tail particles of `params`, and return it. Each
position is a Sobol proposal, uniform on `[0, 2π/κ]`, with the weight `fₓ(x)·2π/κ/N`. The
velocities are drawn as for `AcceptRejectSampling`.
"""
function draw!(dist::ParticleDistribution{<:Any, 1, 1}, fₓ::Base.Callable,
        params::BumpOnTail, ::ImportanceSampling)
    @unpack κ, ε, α, v₀, σ = params

    # number of particles
    N = length(dist.particles)

    # Sobol sampling
    s = SobolSeq(2)
    skip(s, 2N)

    n = 1
    while n < N+1
        ###
        y = Sobol.next!(s)
        x₀ = y[1] * 2π / κ  # proposal, uniform

        # importance sampling in x
        dist.particles.x[1, n] = x₀
        dist.particles.w[1, n] = fₓ(x₀, params) .* 2π ./ κ ./ N    # *L is equivalent to 1/(1/L) i.e. uniform proposal. f is NOT normalized!

        # inverse CDF sampling in v
        vₛ = 2y[2] - 1
        dist.particles.v[1, n] = √2 * erfinv(vₛ)  # draw from bulk
        if rand(1)[1] > 1-α
            dist.particles.v[1, n] = dist.particles.v[1, n]*σ + v₀  # draw from tail
        end
        n = n+1
    end

    return dist
end
