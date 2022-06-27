
using Parameters
using Random
using Sobol
using SpecialFunctions


"""
Returns particles drawn from g for the bump-on-tail case.

input: nb. pf particles Nₚ, x-marginal function gₓ, parameters μ
output: Particles struct
"""
function draw_g_accept_reject(N::Int, fₓ::Function, params::NamedTuple) where {T}
    @unpack κ, ε, a, v₀, σ = params
    
    x = zeros(N); v = zeros(N); w = zeros(N)

    # Sobol sampling
    s = SobolSeq(2); skip(s, 2N)

    n = 1
    while n < N+1
        ###
        y = next!(s)
        x₀ = y[1]*2π/κ  # proposal, uniform

        # accept-reject in x
        if rand(1)[1] > ( (fₓ(x₀, params))/(1 + ε) )
            continue                            # reject
        end
        x[n] = x₀                               # accept
        w[n] = 2π / κ / N

        # inverse CDF sampling in v
        vₛ = 2y[2] - 1
        v[n] = √2 * erfinv(vₛ)  # draw from bulk
        if rand(1)[1] > 1-a
            v[n] = v[n] * σ + v₀  # draw from tail
        end
        n = n+1
    end

    return ParticleList(x, v, w)
end

"""
Returns particles drawn from g for the bump-on-tail case.

input: nb. pf particles Nₚ, x-marginal function gₓ, parameters μ
output: Particles struct
"""
function draw_g_importance_sampling(N::Int, fₓ::Function, params::NamedTuple) where {T}
    @unpack κ, ε, a, v₀, σ = params

    x = zeros(N); v = zeros(N); w = zeros(N)

    # Sobol sampling
    s = SobolSeq(2); skip(s, 2N)

    n = 1
    while n < N+1
        ###
        y = next!(s)
        x₀ = y[1] * 2π / κ  # proposal, uniform

        # importance sampling in x
        x[n] = x₀
        w[n] = fₓ(x₀, params) .* 2π ./ κ ./ N    # *L is equivalent to 1/(1/L) i.e. uniform proposal. f is NOT normalized!

        # inverse CDF sampling in v
        vₛ = 2y[2] - 1
        v[n] = √2 * erfinv(vₛ)  # draw from bulk
        if rand(1)[1] > 1-a
            v[n] = v[n]*σ + v₀  # draw from tail
        end
        n = n+1
    end

    return ParticleList(x, v, w)
end

"""
Returns re-weighted particles according to new parameters.

input: Particles struct P₀, proposal distribution function g, target distribution function f
output: Particles struct
"""
function weight_f(P₀, g, f)
    Px = zero(P₀.x); Pv = zero(P₀.v); Pw = zero(P₀.w)
    Px .= P₀.x; Pv .= P₀.v; Pw .= P₀.w
    P = ParticleList(Px, Pv, Pw)
    P.w .*= f(P₀.x, P₀.v) ./ g(P₀.x, P₀.v)
    return P
end
