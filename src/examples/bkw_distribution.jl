@doc raw"""
    BKWDistribution(t₀; seed = 20261107)

The two-dimensional BKW solution of the homogeneous Landau equation at the time `t₀`, the third 2D
reference case: at the Maxwell molecules of `γ = 0` it is an exact solution, and its departure from
the Maxwellian decays at a rate that `scripts/reference_forms.jl` derives from the moment equations.

At `d = 2` the profile of Carrillo, Hu, Wang and Wu (2020), Appendix A, is

```math
f(v) = \frac{1}{2 \pi K} \exp\left(-\frac{v^2}{2K}\right)
  \left[ \frac{4K - 2}{2K} + \frac{1 - K}{2K^2} \, v^2 \right] ,
\qquad K(t) = 1 - \frac{1}{2} e^{-2t} ,
```

of mass one and temperature one, whose second moment is the identity at every `t`. The paper
states `K(t) = 1 - e^{-t/8}/2` for its kernel `A(z) = (|z|² I - z ⊗ z)/16`; this package's kernel
`U(z) = |z|² (I - ẑ ẑᵀ)` is `16 A`, and the sixteenfold kernel multiplies every rate by sixteen.

`initialize!` fills the two velocity components of a `ParticleDistribution` with a sample of the
profile drawn from `Random.MersenneTwister(seed)`. The profile is a Gaussian of covariance `K I` of
the fraction `(2K - 1)/K`, and that Gaussian weighted by `|v|²` for the rest, whose squared radius is
`2K` times a `Gamma(2, 1)` variate: a sum of two standard exponentials. The sample's means are `0`
and its variances `1` to the standard errors `1/√N` and `√(2/N)` of its `N` particles, and its
fourth moment is that of the profile rather than the Gaussian's `8`. The positions are not touched,
and the weights are set to `1/N`. A non-finite `t₀` throws an `ArgumentError`, and a negative one
too: the quadratic factor of the profile is negative near `v = 0` there, and the profile is no
distribution.
"""
struct BKWDistribution{T} <: VlasovExample
    time::T
    seed::Int

    function BKWDistribution(time::T = 0.2; seed::Integer = 20261107) where {T <: Real}
        isfinite(time) ||
            throw(ArgumentError("t₀ must be finite, got $time"))
        time ≥ 0 ||
            throw(ArgumentError("the BKW profile is a distribution only at the times " *
                                "t₀ ≥ 0, and t₀ = $time is negative"))
        new{T}(time, seed)
    end
end

"""
    bkw_parameter(t) -> Float64

The parameter `K(t) = 1 - exp(-2t)/2` of the BKW profile at the time `t`, at the kernel
`U(z) = |z|² (I - ẑ ẑᵀ)` of `γ = 0`: the paper's `K(t) = 1 - exp(-t/8)/2` at a rate sixteen times
its own, since `U = 16 A` for its kernel `A`. `K(0) = 1/2` is the profile that vanishes
quadratically at `v = 0`, and `K(t) → 1` is the Maxwellian.
"""
bkw_parameter(t) = 1 - exp(-2 * t) / 2

function initialize!(dist::ParticleDistribution, params::BKWDistribution,
        ::SamplingMethod = NoSampling())
    v = dist.particles.v
    w = dist.particles.w
    _check_velocity_dimension(v)
    T = eltype(v)
    npart = size(v, 2)
    rng = MersenneTwister(params.seed)
    K = T(bkw_parameter(params.time))
    # the profile is (1/(2πK)) exp(-v²/(2K)) (A + B v²) with A = (2K - 1)/K and
    # B = (1 - K)/(2K²); of its mass the fraction A is the Gaussian of covariance K I and the rest
    # is that Gaussian weighted by v²/(2K), whose v² is 2K times a Gamma(2, 1) variate
    A = (2K - 1) / K
    for α in 1:npart
        if rand(rng) < A
            v[1, α] = sqrt(K) * randn(rng)
            v[2, α] = sqrt(K) * randn(rng)
        else
            θ = 2 * T(π) * rand(rng)
            r = sqrt(2 * K * (-log(rand(rng)) - log(rand(rng))))
            v[1, α] = r * cos(θ)
            v[2, α] = r * sin(θ)
        end
    end
    w[1, :] .= 1 / npart
    return dist
end
