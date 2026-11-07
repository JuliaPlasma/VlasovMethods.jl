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
realised weights are those counts over `N` and the same seed gives the same sample. The positions
are not touched, and the weights are set to `1/N`.

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
        new{WT, VT, MT}(weights, variances, means, seed)
    end
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
