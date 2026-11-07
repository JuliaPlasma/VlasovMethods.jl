@doc raw"""
    DoubleMaxwellian2D(; weights = (0.5, 0.5), variances = ((1.0, 1.0), (1.0, 1.0)),
                       means = ((0.0, 0.0), (0.0, 0.0)), seed = 20261107)

The sum of two displaced centred Maxwellians of the weights `weights`, the variances `variances`
and the means `means`,

```math
f(v) = w_1 \, g_1(v - m_1) + w_2 \, g_2(v - m_2) ,
```

with `g_c` the Gaussian of the diagonal covariance `variances[c]` and the mean `means[c]`, each a
pair of two velocity components: the 2D reference case whose second moments relax to those of the
mixture. Its exact mean is `w_1 m_1 + w_2 m_2` and its exact variances are `Σ_c w_c (σ²_c + m²_c)`
less the square of that mean.

`initialize!` gives the first component the share `round(weights[1] * N)` of the `N` particles,
draws each from `Random.MersenneTwister(seed)` and gives the rest to the second, so that the
realised weights are those counts over `N` and the same seed gives the same sample within one
Julia version. The positions are not touched, and the weights are set to `1/N`. The constructor
refuses what a mixture of two components cannot be, each with an `ArgumentError`: other than two
weights, variances or means, a variance or mean pair of other than two velocity components, a
weight that is negative or that does not sum with its twin to one, a variance that is not positive
and finite, and a mean that is not finite.

# Example

```julia
dist = initialize!(ParticleDistribution(1, 2, 10_000),
                   DoubleMaxwellian2D(weights = (0.6, 0.4),
                                      variances = ((0.5, 0.3), (0.3, 0.5)),
                                      means = ((0.4, 0.0), (-0.4, 0.0))))
```
"""
struct DoubleMaxwellian2D{WT, VT, MT} <: VlasovExample
    weights::WT
    variances::VT
    means::MT
    seed::Int

    function DoubleMaxwellian2D(weights::WT, variances::VT, means::MT,
            seed::Integer = 20261107) where {WT <: Tuple, VT <: Tuple, MT <: Tuple}
        _check_components(weights, variances, means)
        new{WT, VT, MT}(weights, variances, means, seed)
    end
end

# The mixture is the sum of two components, and its weights are the shares of the particles the two
# take: the first component's count is the weight's share of them and the second takes the rest, so
# that weights of `(0.6, 0.6)` would realise `(0.6, 0.4)` and a negative one would take every
# particle. A variance and a mean are the numbers each component's density is defined by, which
# `_check_variance` and `_check_mean` refuse where they are not.
function _check_components(weights, variances, means)
    all(==(2), (length(weights), length(variances), length(means))) &&
    all(pair -> length(pair) == 2, (variances..., means...)) ||
        throw(ArgumentError("the double Maxwellian is the sum of two components, each of one " *
                            "weight, one pair of variances and one pair of means, and this one " *
                            "has $(length(weights)) weights, $(length(variances)) variance " *
                            "pairs of lengths $(length.(variances)) and $(length(means)) mean " *
                            "pairs of lengths $(length.(means))"))
    all(weight -> weight ≥ 0, weights) && isapprox(sum(weights), 1) ||
        throw(ArgumentError("the weights of the two components must be non-negative and sum to " *
                            "one, and those of this one sum to $(sum(weights))"))
    for c in 1:2
        for d in 1:2
            _check_variance(variances[c][d], "variances[$c][$d]")
            _check_mean(means[c][d], "means[$c][$d]")
        end
    end
    return nothing
end

function DoubleMaxwellian2D(; weights = (0.5, 0.5),
        variances = ((1.0, 1.0), (1.0, 1.0)), means = ((0.0, 0.0), (0.0, 0.0)),
        seed::Integer = 20261107)
    return DoubleMaxwellian2D(weights, variances, means, seed)
end

function initialize!(dist::ParticleDistribution, params::DoubleMaxwellian2D,
        ::SamplingMethod = NoSampling())
    v = dist.particles.v
    w = dist.particles.w
    _check_velocity_dimension(v)
    T = eltype(v)
    npart = size(v, 2)
    rng = MersenneTwister(params.seed)
    # the number of particles of the first component is its weight's share, rounded, and the second
    # component takes the rest, so that the realised weights are these two counts over `npart`
    nfirst = round(Int, params.weights[1] * npart)
    for α in 1:npart
        c = α ≤ nfirst ? 1 : 2
        for d in 1:2
            v[d, α] = T(params.means[c][d]) +
                      sqrt(T(params.variances[c][d])) * randn(rng)
        end
    end
    w[1, :] .= 1 / npart
    return dist
end
