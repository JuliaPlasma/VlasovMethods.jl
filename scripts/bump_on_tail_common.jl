# The code the four bump-on-tail scripts share: the field every march takes, the two sweeps, the
# growth-rate fit of the paper's item 28, and the writer of P39's two input files.
#
# It loads no Makie package, so that `test/integration/bump_on_tail_scripts.jl` can `include` it
# with the test dependencies alone; a script that draws loads CairoMakie itself.
#
# The paper's model is the bump-on-tail of `κ = 0.3` on the periodic domain `2π/κ`, and its sweep
# parameter is `χ = κₛ/κ` of the reference's `bump_on_tail_1_training.jl`: one physical run, the
# same particles and the same potential, marched at ten time scales.

using HDF5
using LinearAlgebra
using PoissonSolvers
using Random
using ReducedBasisMethods
using ReducedComplexityModeling: ParameterSpace, Parameter
using VlasovMethods
using VlasovMethods: efield!, update!   # the field interface, which the package does not export

# The paper's parameters, from the reference's `bump_on_tail_1_training.jl`.
const BOT_KAPPA = 0.3       # the wave number of the physical run
const BOT_DT = 0.1          # timestep of the paper command
const BOT_T = 25.0          # final time of the paper command
const BOT_NP = 10_000       # particles of the paper command
const BOT_NH = 16           # periodic spline cells of the paper command
const BOT_ORDER = 3         # spline degree of the paper command
const BOT_SEED = 1234       # the seed of the particle draw and of the test sweep

# The ten wave numbers of the paper's sweeps, the `[0.1, 0.5]` of the reference's
# `χ = Parameter(:χ, 0.1/κ, 0.5/κ, 10)`, written as an explicit list of values rather than as a
# `Parameter` range, so that each sample's `χ` is computed from a number this file names (K56
# records why a one-sample range cannot carry one at all).
const TRAINING_KAPPA = [0.1 + 0.4 * j / 9 for j in 0:9]
const KAPPA_RANGE = (0.1, 0.5)

# The fit of item 28: the reference's `get_regression_αβ(t, W, 5)` with `m = m[1:4]` when four
# maxima are found, and every maximum when fewer are.
const GROWTH_WINDOW = 5
const GROWTH_MAXIMA = 4

"The number of steps of a march of the final time `T` at the timestep `dt`, `round(Int, T/dt)`."
step_count(T, dt) = round(Int, T / dt)

"""
    test_sweep(seed = BOT_SEED) -> (κ, χ)

The ten test wave numbers of the reference's `bump_on_tail_3_testing.jl`, drawn uniform in
`[0.1, 0.5]` after `Random.seed!(seed)` and sorted, and the time scales `χ = κ/0.3` they stand for.

The seed is a keyword because the draw consumes the global generator: a sweep of another seed is
another set of wave numbers.
"""
function test_sweep(seed = BOT_SEED)
    Random.seed!(seed)
    lo, hi = KAPPA_RANGE
    κ = sort(lo .+ (hi - lo) .* rand(10))
    return (; κ = κ, χ = κ ./ BOT_KAPPA)
end

"""
    vp_field(potential, params)

The field of a bump-on-tail march at the scale `χ` of `params`: the [`ScaledField`](@ref) of the
[`PoissonField`](@ref) of `potential`. Every march of these scripts takes its field from here, with
the same `params` it marches, so that the step's scale and the field's `χ²` are one number — the
march refuses a `ScaledField` at any other `χ` (D-P38-8).
"""
vp_field(potential, params) = ScaledField(PoissonField(potential), params.χ)

"The potential of the paper's run: `nh` periodic spline cells of degree `order` on the one period
`2π/κ`."
function bot_potential(nh, order; κ = BOT_KAPPA)
    Potential(PeriodicBasisSpline((0.0, 2π / κ), order, nh))
end

"""
    bump_on_tail_particles(np, seed = BOT_SEED) -> ParticleDistribution

The paper's initial condition: `np` particles of the bump-on-tail of `BOT_KAPPA`, drawn by the
accept–reject sampler of `src/examples/bumpontail.jl` after `Random.seed!(seed)`. It is drawn once
and marched at every `χ` of a sweep, as the reference does.
"""
function bump_on_tail_particles(np, seed = BOT_SEED)
    Random.seed!(seed)
    return initialize!(ParticleDistribution(1, 1, np), BumpOnTail(κ = BOT_KAPPA))
end

# A number as a file name: `0.1` is `0p1`, so that a run id carries no second dot.
_number(x::Real) = replace(string(x), "." => "p")

"""
    run_id(params) -> String

The name of a run, from the parameters that change its results: the particle and cell counts, the
timestep and the final time, and the wave-number range of the sweep. Two runs of different
parameters never share a file.
"""
function run_id(params)
    return "BoT_np$(params.np)_nh$(params.nh)_dt$(_number(params.dt))_T$(_number(params.T))" *
           "_k$(_number(params.κlo))-$(_number(params.κhi))"
end

# The number of maxima of `W` at the fit's window, which selects between `1:4` and fewer.
_maxima_count(W) = length(VlasovMethods._find_maxima(W, GROWTH_WINDOW))

"""
    fit_growth_rate(t, W)

The growth rate `γ` of the field energy history `W` at the times `t`, by the fit of the reference's
`VP-ROM-Plotting.ipynb` cell 36: half of the least-squares slope of `log W` through the first four
maxima of `W` found with a window of five samples, and through every maximum when fewer than four
are found. It is `growth_rate(t, W; window = 5, maxima = 1:m)` with that `m`; `growth_rate`'s own
defaults, `window = 2` and `maxima = 2:5`, are not the reference's fit.

The fit of a history with fewer than two maxima is `growth_rate`'s own `ArgumentError`, which names
the count found.
"""
function fit_growth_rate(t, W)
    m = min(GROWTH_MAXIMA, _maxima_count(W))
    return growth_rate(t, W; window = GROWTH_WINDOW, maxima = 1:m)
end

"""
    fill_accelerations!(A, X, potential, χ, params, w)

Write the accelerations of every saved state of the positions `X` into `A`, because
[`integrate_vp`](@ref) saves none: the positions of each saved state are deposited with the weights
`w` into [`vp_field`](@ref)`(potential, (χ = χ,))` again, and the field it solves is written into
the sample. The values are those of the scaled state, `E/χ²`, as the march's kicks are.

The potential is left holding the solution of the last saved state.
"""
function fill_accelerations!(A, X, potential, χ, params, w)
    field = vp_field(potential, (χ = χ,))
    a = zeros(size(X, 2))
    for p in axes(X, 4), ts in axes(X, 3)

        x = @view X[1, :, ts, p]
        update!(field, x, w, (ts - 1) * params.dt)
        efield!(field, a, x)
        A[1, :, ts, p] .= a
    end
    return A
end

"""
    training_sweep(params; κ = TRAINING_KAPPA) -> (snapshots, kappa, chi, potentials, particles)

The ten training samples of the paper's bump-on-tail set: the particles of `params.seed`, one
potential on the one period, and one full-order march per entry of `κ`, each at the time scale
`χ = κ/0.3` and saved at every step of `params.dt` over `[0, params.T]`.

The accelerations of the saved states are filled by [`fill_accelerations!`](@ref), which the march
does not write.
"""
function training_sweep(params; κ = TRAINING_KAPPA)
    nₜ = step_count(params.T, params.dt)
    nₛ = nₜ + 1
    χ = κ ./ BOT_KAPPA
    space = ParameterSpace(Parameter(:χ, χ))
    SS = Snapshots(Float64, 1, params.np, params.nh, nₛ, space)
    ip = IntegratorParameters(params.dt, nₜ, nₛ, params.nh, params.np, length(κ))

    dist = bump_on_tail_particles(params.np, params.seed)
    potentials = [bot_potential(params.nh, params.order) for _ in κ]

    for p in eachindex(κ)
        march_params = (χ = χ[p],)
        integrate_vp(dist.particles, march_params, vp_field(potentials[p], march_params),
            SS, ip, p)
    end

    for p in eachindex(κ)
        fill_accelerations!(SS.A, @view(SS.X[:, :, :, p:p]), potentials[p], χ[p], params,
            vec(dist.particles.w))
    end

    return (; snapshots = SS, kappa = κ, chi = χ, potentials, particles = dist)
end

"""
    test_sample(params, χ) -> (snapshots, field, particles, potential, ip)

One sample of the test sweep at the time scale `χ`: the particles of `params.seed`, their
potential, and the full-order march over `[0, params.T]` saved at every step. `field` is the
[`vp_field`](@ref) of the sample, which the reduced march of the same sample takes; the snapshot
holds the one sample, and `ip` is its integrator parameters.
"""
function test_sample(params, χ)
    nₜ = step_count(params.T, params.dt)
    nₛ = nₜ + 1
    space = ParameterSpace(Parameter(:χ, [χ]))
    SS = Snapshots(Float64, 1, params.np, params.nh, nₛ, space)
    ip = IntegratorParameters(params.dt, nₜ, nₛ, params.nh, params.np, 1)

    dist = bump_on_tail_particles(params.np, params.seed)
    potential = bot_potential(params.nh, params.order)
    field = vp_field(potential, (χ = χ,))
    integrate_vp(dist.particles, (χ = χ,), field, SS, ip, 1)

    return (; snapshots = SS, field, particles = dist, potential, ip)
end

# The three root attributes of a P39 file, written into the file `io` rather than into a group.
# `attributes` is HDF5's and not Makie's, which a drawing script loads as well.
function _attributes!(io, commit, command, seed)
    attrs = HDF5.attributes(io)
    attrs["commit"] = commit
    attrs["command"] = command
    attrs["seed"] = seed
    return io
end

"""
    write_p39_files(dir, params, particle, field, X, A; commit, command, seed) -> (ppath, vpath)

Write P39's two input files of the run `params` into `dir`:

- `<runid>_projections.h5`, the particle basis at `/particle` and the field basis at `/field`;
- `<runid>_vectorfield.h5`, the reduced positions `X` and accelerations `A`, each
  `(k, nₛ, nparam)` for the particle basis's `k` columns.

Both files carry the root attributes `commit`, `command` and `seed`, so that a reader can say which
checkout and which command produced them.

# Throws

- `ArgumentError` when `A` does not have the shape of `X`, and when `X` does not hold one row per
  column of the particle basis.
"""
function write_p39_files(dir, params, particle, field, X, A; commit, command, seed)
    size(X) == size(A) || throw(ArgumentError(
        "the reduced positions hold $(size(X)) and the reduced accelerations $(size(A))"))
    size(X, 1) == particle.k || throw(ArgumentError(
        "the reduced state holds $(size(X, 1)) rows; the particle basis holds " *
        "$(particle.k) columns"))

    id = run_id(params)
    ppath = joinpath(dir, "$(id)_projections.h5")
    vpath = joinpath(dir, "$(id)_vectorfield.h5")

    h5open(ppath, "w") do io
        h5save(io, particle; path = "particle")
        h5save(io, field; path = "field")
        _attributes!(io, commit, command, seed)
    end

    h5open(vpath, "w") do io
        io["X"] = X
        io["A"] = A
        _attributes!(io, commit, command, seed)
    end

    return ppath, vpath
end

"The commit of the checkout at `root`, as `git rev-parse HEAD` prints it."
function git_commit(root)
    return strip(read(`git -C $root rev-parse HEAD`, String))
end

"The name of the figure of an acceptance item and a quantity: `item24_trajectory.png`."
function item_figure(item::Integer, what::AbstractString)
    "item" * lpad(item, 2, '0') * "_" * what * ".png"
end
