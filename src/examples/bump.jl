
struct Bump{T} <: VlasovExample
    a::T # steepness of bump
    b::T # width parameters

    function Bump(; a::T = 10.0, b::T = 0.5) where {T}
        new{T}(a, b)
    end
end

# v-part of distribution function
function fᵥ(v::T, params::Bump) where {T}
    @unpack a, b = params
    @. (1/(4b) * (tanh(a * (v + b)) - tanh(a * (v - b)))) # normalised to integrate to 1.
end

function initialize!(dist::ParticleDistribution, params::Bump,
        sampling::SamplingMethod = AcceptRejectSampling())
    draw!(dist, fᵥ, params, sampling)
end

"""
    draw!(dist, fᵥ, params::Bump, ::AcceptRejectSampling)

Fill the velocities of the 1D1V distribution `dist` with the bump of `params`, and return it.
Each velocity is a uniform proposal on `[-b, b]`, accepted with the probability
`2b fᵥ(v) / 1.5`, and each particle has the weight `1/N`. The positions keep their values.
"""
function draw!(dist::ParticleDistribution{<:Any, 1, 1}, fᵥ::Base.Callable, params::Bump,
        ::AcceptRejectSampling)
    @unpack a, b = params

    # number of particles
    N = length(dist.particles)

    # # Sobol sampling
    # s = SobolSeq(2); skip(s, 2N)

    n = 1
    while n < N+1
        y = (rand(1)[1] - 0.5) * 2 * b

        ###
        # y = Sobol.next!(s)
        # x₀ = y[1] - b  # proposal, uniform
        # x₀ = y  # proposal, uniform shifted to center around zero

        # accept-reject in v
        M = 1.5
        if rand(1)[1] > ((fᵥ(y, params))/(M / (2*b)))
            continue                            # reject
        end
        dist.particles.v[1, n] = y                             # accept
        dist.particles.w[1, n] = 1 / N

        n = n+1
    end

    return dist
end
