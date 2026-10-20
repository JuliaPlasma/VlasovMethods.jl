# The code the three bump-on-tail scripts share: the parameters of the model, the wave numbers of
# a sweep, the field of a sample, the particle cloud and the name of a run.
#
# It uses VlasovMethods and the standard library only, so that a test can `include` it.

using Printf: @sprintf
using Random
using VlasovMethods: BumpOnTail, ParticleDistribution, PoissonField, ScaledField,
                     coefficients,
                     initialize!

"""
    κ_ref

The reference wave number of the bump-on-tail model, `0.3`: the `χ` of a sample is its wave number
divided by this, as the training notebook computes it.
"""
const κ_ref = 0.3

"""
    bump_on_tail_parameters()

The parameters of the bump-on-tail model: the wave number `κ = κ_ref`, the amplitude `ε = 0.03`, the
share `α = 0.1` of the fast particles, the temperature `σ = 0.5` and the offset velocity `v₀ = 4.5`
of `src/examples/bumpontail.jl`.
"""
bump_on_tail_parameters() = BumpOnTail(ε = 0.03, κ = κ_ref, α = 0.1, σ = 0.5, v₀ = 4.5)

"""
    domain()

The domain of the bump-on-tail model: one period of the perturbation at the reference wave number.
"""
domain() = (0.0, 2π / κ_ref)

"""
    march_parameters(κ)

The parameters of a march of the sample at the wave number `κ`: its `χ` scale, which is `κ` divided
by the reference wave number. The `χ` is the one source of the scale of the march.
"""
march_parameters(κ) = (χ = κ / κ_ref,)

"""
    step_count(final_time, tstep)

The number of time steps of a march from zero to `final_time`: `round` and not the
floating-point `div` of the old script, which gave 249 for `T = 25` and `dt = 0.1`. The
integrator parameters need `nₛ - 1` to divide `nₜ`, and `nₛ = nₜ + 1` divides it.
"""
step_count(final_time, tstep) = round(Int, final_time / tstep)

"""
    vp_field(potential, params)

The electric field of a sample: a `ScaledField` over the `PoissonField` of the `potential`, at the
`χ` of the parameters. A march takes the scale from `params.χ` and refuses a field that carries
another one, so the field is built at the march's `χ`.
"""
vp_field(potential, params) = ScaledField(PoissonField(potential), params.χ)

"""
    potential_degrees(potential)

The number of degrees of freedom of the `potential`: the length of the coefficient vector that a
march saves into the snapshots.
"""
potential_degrees(potential) = length(coefficients(PoissonField(potential)))

"""
    initial_distribution(n; seed = 1234)

The particle cloud of the bump-on-tail model: `n` particles drawn by the accept–reject method of
`src/examples/bumpontail.jl` from the seeded generator.
"""
function initial_distribution(n; seed = 1234)
    Random.seed!(seed)
    dist = ParticleDistribution(1, 1, n)
    initialize!(dist, bump_on_tail_parameters())
    return dist
end

"""
    κ_sweep(κ_min, κ_max, n)

The `n` wave numbers of a training sweep, evenly spaced between `κ_min` and `κ_max`.
"""
κ_sweep(κ_min, κ_max, n) = collect(range(κ_min, κ_max; length = n))

"""
    κ_test(κ_min, κ_max, n; seed = 1234)

The `n` wave numbers of a test sweep, drawn uniformly from `[κ_min, κ_max]` and sorted.
"""
function κ_test(κ_min, κ_max, n; seed = 1234)
    Random.seed!(seed)
    return sort(rand(n) .* (κ_max - κ_min) .+ κ_min)
end

"""
    bump_on_tail_runid(npart, nsamples, final_time)

The name of a bump-on-tail run: the particle count, the wave-number range of the sweep, the number
of samples and the final time. The name carries the parameters it names, so that a run of other
parameters is not mistaken for this one.
"""
function bump_on_tail_runid(npart, nsamples, final_time)
    return "BoT_Np$(npart)_k_010_050_ns$(nsamples)_T" * @sprintf("%.4g", final_time)
end

# The defaults of the three scripts, so that they name the same run: the training script marches
# with them, and the other two read the files it writes.
const DEFAULT_NPART = 1000
const DEFAULT_NSAMPLES = 10
const DEFAULT_FINAL_TIME = 25.0
const DEFAULT_RUNID = bump_on_tail_runid(DEFAULT_NPART, DEFAULT_NSAMPLES, DEFAULT_FINAL_TIME)

"""
    repository_commit(directory = dirname(@__DIR__))

The commit of the VlasovMethods checkout the script runs in, or `"unknown"` when it runs outside
one. The `directory` is the root of the checkout, which is the parent of the `scripts` directory.
"""
function repository_commit(directory = dirname(@__DIR__))
    try
        return readchomp(`git -C $directory rev-parse HEAD`)
    catch
        return "unknown"
    end
end
