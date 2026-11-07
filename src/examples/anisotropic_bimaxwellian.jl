@doc raw"""
    AnisotropicBiMaxwellian(σ1², σ2²; seed = 20261107)

The centred bi-Maxwellian of the two temperatures `σ1²` and `σ2²`,

```math
f(v) = \frac{1}{2 \pi \sigma_1 \sigma_2}
  \exp\left(-\frac{v_1^2}{2 \sigma_1^2} - \frac{v_2^2}{2 \sigma_2^2}\right) ,
```

the 2D reference case whose second moments isotropise: at the Maxwell molecules of `γ = 0` the
difference `T₁₁ − T₂₂` of its second moments decays at the rate `8`, which
`scripts/reference_forms.jl` derives from the moment equations.

`initialize!` fills the two velocity components of a `ParticleDistribution` with a sample of `f`
drawn from `Random.MersenneTwister(seed)`, one stream of `N` normal variates per component, so that
the same seed gives the same sample. The weights and the positions of the distribution are not
touched: the second moments of a sample are those of the weights it already holds. The sample's
means, variances and covariance approach the exact `0`, `0`, `σ1²`, `σ2²` and `0` at the standard
errors `σ/√N`, `σ²√(2/N)` and `σ₁σ₂/√N` of its `N` particles.
"""
struct AnisotropicBiMaxwellian{T} <: VlasovExample
    variance_one::T
    variance_two::T
    seed::Int

    function AnisotropicBiMaxwellian(variance_one::T, variance_two::T;
            seed::Integer = 20261107) where {T <: Real}
        new{T}(variance_one, variance_two, seed)
    end
end

# The three reference cases of `src/examples/` are two-dimensional. A distribution of one velocity
# component cannot hold them, and one of three would keep its third component at zero, which reads
# as a temperature of zero: both are refused before any draw.
@noinline function _throw_velocity_dimension(vdim)
    throw(DimensionMismatch("the 2D reference cases fill two velocity components, and this " *
                            "distribution has $vdim"))
end
function _check_velocity_dimension(v::AbstractMatrix)
    size(v, 1) == 2 || _throw_velocity_dimension(size(v, 1))
    return nothing
end

function initialize!(dist::ParticleDistribution, params::AnisotropicBiMaxwellian,
        ::SamplingMethod = NoSampling())
    v = dist.particles.v
    _check_velocity_dimension(v)
    T = eltype(v)
    npart = size(v, 2)
    rng = MersenneTwister(params.seed)
    v[1, :] .= sqrt(T(params.variance_one)) .* randn(rng, npart)
    v[2, :] .= sqrt(T(params.variance_two)) .* randn(rng, npart)
    return dist
end
