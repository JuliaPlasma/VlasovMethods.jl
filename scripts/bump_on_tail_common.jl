# What the drivers of the paper's two example sets share: the parameters of the bump-on-tail set,
# its two sweeps of the wave number, its initial condition, the field of one sample and the
# accelerations that field applies, the stride of the analytic set's saves, and the commit and the
# command that a written file records.
#
# It uses VlasovMethods and the standard library only, so that
# `test/integration/bump_on_tail_scripts.jl` includes it without the environment of `scripts/`.

using LinearAlgebra
using Random
using VlasovMethods

"The wave number `κ` of the reference problem, on whose domain every sample of a sweep runs."
const KAPPA_REF = 0.3

"The bump-on-tail initial condition of ReducedBasisMethods' scripts: the amplitude, the fast
share, the width and the offset velocity of the bump."
const BUMP = (ε = 0.03, α = 0.1, σ = 0.5, v₀ = 4.5)

"""
    BUMP_ON_TAIL_DEFAULTS

The defaults of the three bump-on-tail drivers, which are the parameters of ReducedBasisMethods'
scripts: `10⁴` particles, `16` cells of a cubic periodic spline, `dt = 0.1` up to `T = 25`, ten
samples of `κ ∈ [0.1, 0.5]`, the seed `1234`, and the tolerances `1e-4` of the particle basis and
`1e-2` of the field basis, the pair its projection script runs.
"""
const BUMP_ON_TAIL_DEFAULTS = (runid = "bump_on_tail", npart = 10_000, ncells = 16,
    degree = 3, tstep = 0.1, final_time = 25.0, kappa_min = 0.1, kappa_max = 0.5,
    nsamples = 10, seed = 1234, particle_tol = 1e-4, field_tol = 1e-2)

"The periodic domain `(0, 2π/κ)` of one wavelength of the perturbation of wave number `κ`."
domain(κ) = (0.0, 2π / κ)

"""
    step_count(final_time, tstep)

The number of steps of `tstep` up to `final_time`, `round(Int, final_time / tstep)`.
ReducedBasisMethods' training script takes `Int(div(25, 0.1))`, which is `249`, because `div` of
two floats rounds the quotient `249.99…` down.
"""
step_count(final_time, tstep) = round(Int, final_time / tstep)

"""
    training_kappas(kappa_min, kappa_max, nsamples)

The wave numbers of the training sweep: `nsamples` values evenly spaced over
`[kappa_min, kappa_max]`, both ends included.
"""
function training_kappas(kappa_min, kappa_max, nsamples)
    return collect(range(kappa_min, kappa_max; length = nsamples))
end

"""
    testing_kappas(kappa_min, kappa_max, nsamples, seed)

The wave numbers of the test sweep: `nsamples` values uniform in `[kappa_min, kappa_max)`, drawn
from `Xoshiro(seed)` and sorted, as ReducedBasisMethods' test script draws them after
`Random.seed!(seed)`.
"""
function testing_kappas(kappa_min, kappa_max, nsamples, seed)
    κ = rand(Xoshiro(seed), nsamples) .* (kappa_max - kappa_min) .+ kappa_min
    return sort!(κ)
end

"""
    march_parameters(κ)

The parameters of the march of the sample at the wave number `κ`: `(χ = κ / KAPPA_REF,)`, the
scale by which the sample's time and field differ from the reference problem's.
"""
march_parameters(κ) = (χ = κ / KAPPA_REF,)

"""
    initial_distribution(npart, seed)

The `npart` particles of the bump-on-tail initial condition at `KAPPA_REF`, drawn by the
accept–reject method of `BumpOnTail` after `Random.seed!(seed)`. The draw reads the global
generator, so the seed is set here and not passed on.
"""
function initial_distribution(npart, seed)
    Random.seed!(seed)
    return initialize!(ParticleDistribution(1, 1, npart),
        BumpOnTail(; ε = BUMP.ε, κ = KAPPA_REF, α = BUMP.α, σ = BUMP.σ, v₀ = BUMP.v₀))
end

"""
    vp_field(potential, params)

The field of one sample: the Poisson field of `potential`, scaled by `params.χ`. Each sample
builds its own, from the parameters that its march takes, so the march finds the `χ` it marches
at and the accelerations of [`accelerations!`](@ref) are the ones the march applied, `E/χ²`.
"""
vp_field(potential, params) = ScaledField(PoissonField(potential), params.χ)

"""
    accelerations!(SS, field, w, t, p)

Write the acceleration of `field` at each saved position of the sample `p` into `SS.A`, with the
particle weights `w` and the saved times `t`. The marches write `X`, `V`, `Φ`, `W`, `K` and `M`,
and not `A`, which the field basis is built from.
"""
function accelerations!(SS, field, w, t, p)
    for ts in axes(SS.X, 3)
        x = SS.X[1, :, ts, p]
        VlasovMethods.update!(field, x, w, t[ts])
        VlasovMethods.efield!(field, view(SS.A, 1, :, ts, p), x)
    end
    return SS
end

"""
    save_stride(save_step, dt, rdt)

The number of reduced steps of `rdt` between two saves, which are `save_step` steps of `dt`
apart: the ratio `save_step * dt / rdt`, which must be a whole number. Otherwise the saved times
of the reduced run fall between those of the full model, and an `ArgumentError` names the ratio.

The ratio is not computed exactly. `dt` and `rdt` are each the rounding of a decimal, and the
product and the quotient are each rounded once more: four roundings of at most `u = eps/2`
relative each, so the computed ratio is within `4u = 2eps` of the exact one, relative. A ratio
within `2 * eps(Float64) * stride` of the whole number `stride` is therefore taken as `stride`,
because the saved times then differ from `stride * rdt` only by the round-off of the times
themselves. `save_stride(3, 0.7, 0.084)` computes the ratio `24.999999999999993`, which is
`1.28 eps` from `25`, relative, and is `25`.

A ratio that is not finite, such as the ratio of `rdt = 0` or of a `NaN` input, and a ratio
whose magnitude an `Int` does not hold are refused by an `ArgumentError` of their own.
"""
function save_stride(save_step, dt, rdt)
    ratio = save_step * dt / rdt
    abs(ratio) < typemax(Int) || throw(ArgumentError(
        "a save of $(save_step) steps of dt = $(dt) is $(ratio) reduced steps of " *
        "rdt = $(rdt), which is no count of reduced steps"))
    stride = round(Int, ratio)
    stride ≥ 1 || throw(ArgumentError(
        "a save of $(save_step) steps of dt = $(dt) is $(ratio) reduced steps of " *
        "rdt = $(rdt), and a save is at least one reduced step"))
    abs(ratio - stride) ≤ 2 * eps(Float64) * stride || throw(ArgumentError(
        "a save of $(save_step) steps of dt = $(dt) is $(ratio) reduced steps of " *
        "rdt = $(rdt), which is not a whole number"))
    return stride
end

"""
    repository_commit()

The commit of the VlasovMethods checkout that holds this file, `git rev-parse HEAD`. A failure of
`git` is the caller's error: a written file is not to record a commit it does not know.
"""
repository_commit() = readchomp(`git -C $(dirname(@__DIR__)) rev-parse HEAD`)

"""
    command_line(script, args)

The command of a run: the path of the script under `scripts/` and its arguments, as a written
file records it.
"""
function command_line(script, args)
    join(["julia", "--project=scripts", "scripts/" * script, args...], " ")
end
