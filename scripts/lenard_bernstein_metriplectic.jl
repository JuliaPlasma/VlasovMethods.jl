# The three tests of the metriplectic Lenard–Bernstein manuscript (MP-LB): the relaxation of a
# shifted normal, a double Maxwellian and a uniform distribution under
# `MetriplecticLenardBernstein`, with either metric bracket, at the parameters of the
# Lenard–Bernstein manuscript's figures (PD-LB), and the rate check of the temperature under
# frozen coefficients.
#
#   julia --project=scripts scripts/lenard_bernstein_metriplectic.jl [--quick] \
#       [--case shifted-normal|double-maxwellian|uniform] [--bracket dougherty|jeyakumar] \
#       [--npart n] [--nbasis 41] [--tstep h] [--final-time t] [--nu 1.0] [--tolerance 0.05] \
#       [--seed 1] [--runs-dir dir] [--results-dir dir]
#
# `--case` picks the initial condition and its defaults, those of the PD-LB driver of the same
# distribution: `shifted-normal` draws 1000 velocities from the normal distribution of mean 2 and
# variance 1 (`ShiftedNormalV`, which seeds the random number generator with 1234 itself) and runs
# to `t = 1`; `double-maxwellian` draws 1000 velocities, half from the normal distribution of mean
# 2 and half from that of mean −2 (`DoubleMaxwellian`), and runs to `t = 5`; both step with
# `h = 8e-4`. `uniform` draws 200 velocities from the uniform distribution on [−2, 2]
# (`UniformDistribution`) and runs to `t = 1` with `h = 1e-4`. The last two seed the global random
# number generator with `--seed` first; `shifted-normal` takes no `--seed`. Every case takes 41
# cubic basis functions on [−10, 10] and `ν = 1`. `--bracket` picks the metric bracket,
# `DoughertyBracket` (the default) or `JeyakumarBracket`. A case or a bracket that is not one of
# these throws an `ArgumentError` that names the flag and the value, before any run.
#
# Each step is the model's implicit midpoint step, `Picard_iterate_over_particles!`, solved to
# the absolute residual `SOLVER_SCALE · eps · ‖v⁰‖₂`. The bracket conserves the energy, so `‖v‖₂`
# stays at its initial value and this is a fixed multiple of the round-off of the residual's norm.
#
# A run writes, with the file stem `lenard_bernstein_metriplectic_<case>_<bracket>`, the record
# `<stem>.hdf5` into `--runs-dir`, with the times `t`, the velocities `z`, one column per time,
# and the mass, the momentum, the energy and the entropy `S_h`, and three figures into
# `--results-dir`: `<stem>_initial_final_distribution.pdf` (`plot_collision_distribution`),
# `<stem>_conservation.pdf`, the relative drift of the mass, the momentum and the energy
# (`plot_drift`), and `<stem>_entropy_log.pdf`, the relative change `|S_h(t) − S_h(0)| / |S_h(0)|`
# of the entropy on a logarithmic axis (`plot_drift`). Where a step throws — the `DomainError` of
# the projection or the `ErrorException` of a Picard solve that does not converge — the run stops:
# the record is written with the attribute `stopped_at`, the time of that step, the steps reached
# and, after a rate check, its datasets, no figure is written, and the script exits non-zero.
#
# The rate check runs with `--case shifted-normal`, before the run itself. Its bracket keeps the
# density `n₀` and the mean `u₀` of the initial particles and the temperature `T₀ = 2` in place of
# the moments of the current ones: the frozen coefficients of the Ornstein–Uhlenbeck limit. In the
# continuum limit, where `S_γ / w_γ` is `∂ᵥ log f` at the particle, the bracket's equation of motion
# is then `v̇ = −ν T₀ ∂ᵥ log f − ν (v − u₀)` for the Dougherty form and the same times `n₀ / T₀`
# for the Jeyakumar form, so the temperature relaxes to `T₀` as `T(t) − T₀ ∝ exp(−λ t)` with the
# closed form `λ = 2ν` for the Dougherty form and `λ = 2ν n₀ / T₀ = ν n₀` for the Jeyakumar form.
# The weights of the initial particles sum to `n₀ = 1`, so the two closed forms differ by two; at
# `T₀ = n₀` they would be equal. The check dilates the initial particles about their mean to the
# temperature `T₀`, contracts them by `1 / √2` to `T₀ / 2`, takes the steps of the run with the
# implicit midpoint rule, fits `λ` to `ln|T(t) − T₀|` by least squares, and prints the fitted
# value, the closed form and their relative difference, the control, the other form's rate, and
# its relative difference to the fit, and the count of the steps whose fixed-point solve stops at
# its iteration cap (`relax_frozen`). The script exits non-zero, after it has written its files,
# where the relative difference exceeds `--tolerance` or the control's does not. The table goes to
# `<stem>_rate.md` in `--results-dir`, the values to attributes of the record, and the times and
# the temperatures of the check to its datasets `rate_t` and `rate_temperature`. The other two
# cases run no rate check: on the uniform case the fit of one bracket misses its closed form by
# more than `5 %` (`docs/src/collisions.md`).
#
# `--quick` takes 200 particles and the final time 0.05 for `shifted-normal`, the final time 0.05
# for `double-maxwellian`, and the step 1e-3 and the final time 0.05 for `uniform`. The
# double Maxwellian keeps its 1000 particles: with 200, the Picard solve of the Jeyakumar form
# does not converge in the step to `t = 0.016`.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

using LinearAlgebra

const NAME = "lenard_bernstein_metriplectic"

const CASES = Dict(
    "shifted-normal" =>
        (example = ShiftedNormalV((-5.0, 5.0), 2.0), xdomain = (-5.0, 5.0),
            paper = (npart = 1000, nbasis = 41, tstep = 8e-4, final_time = 1.0, nu = 1.0),
            quick = (npart = 200, final_time = 0.05)),
    "double-maxwellian" =>
        (example = DoubleMaxwellian(; domain = (-5.0, 5.0), shift = 2.0),
            xdomain = (-5.0, 5.0),
            paper = (npart = 1000, nbasis = 41, tstep = 8e-4, final_time = 5.0, nu = 1.0,
                seed = 1),
            quick = (final_time = 0.05,)),
    "uniform" =>
        (example = UniformDistribution((0.0, 1.0), (-2.0, 2.0)), xdomain = (0.0, 1.0),
            paper = (npart = 200, nbasis = 41, tstep = 1e-4,
                final_time = 1.0, nu = 1.0, seed = 1),
            quick = (tstep = 1e-3, final_time = 0.05)))

const BRACKETS = Dict("dougherty" => DoughertyBracket(), "jeyakumar" => JeyakumarBracket())

# The case of the rate check, and the temperature `T₀` at which it freezes the coefficients: away
# from the density `n₀ = 1` of the unit weights, so that the two closed forms differ by two.
const RATE_CASE = "shifted-normal"
const RATE_TEMPERATURE = 2.0

# The multiple of `eps · ‖v⁰‖₂` that the residual of each Picard solve must reach.
const SOLVER_SCALE = 8

# The fixed-point solve of a step of the rate check: its damping, the residual relative to the
# iterate at which it stops, and the most iterations it takes.
const RATE_SOLVER_DAMPING = 0.5
const RATE_SOLVER_TOLERANCE = 1e-12
const RATE_SOLVER_ITERATIONS = 1000

"""
    flag_value(args, flag, default)

The argument after the first `flag` in `args`, or `default` where `flag` is not given or is the
last argument. `parse_arguments` checks the command line itself afterwards.
"""
function flag_value(args, flag, default)
    i = findfirst(==(flag), args)
    return i === nothing || i == length(args) ? default : args[i + 1]
end

"""
    march(case, bracket, options; record, data = (;)) -> NamedTuple

The run of `MetriplecticLenardBernstein` with `bracket` from `options.npart` particles of
`case.example`, seeded with `options.seed` where the case takes one, on the cubic spline of
`options.nbasis` functions, for `nsteps(options.final_time, options.tstep)` steps. Returns the
initial and the final particles, the times `t`, the velocities `z` and the `diagnostics`.

Where a step throws a `DomainError` or an `ErrorException`, the record `record` is written with
the attribute `stopped_at`, the time of that step, the times and the velocities of the steps
reached, and the datasets `data`, and the error is rethrown.
"""
function march(case, bracket, options; record, data = (;))
    dist = initial_distribution(case.example, options.npart; seed = get(options, :seed, nothing))
    initial = deepcopy(dist)
    sdist = spline_distribution(options.nbasis)
    model = MetriplecticLenardBernstein(dist, CollisionEntropy(sdist); ν = options.nu, bracket)

    n = nsteps(options.final_time, options.tstep)
    z = zeros(options.npart, n + 1)
    z[:, 1] = dist.particles.v[1, :]
    abstol = SOLVER_SCALE * eps() * norm(z[:, 1])
    history = zeros(options.npart, 2)
    v, vn, vn_minus_one = similar(z[:, 1]), z[:, 1], z[:, 1]
    reached = 0
    try
        for k in 1:n
            Picard_iterate_over_particles!(
                v, vn, vn_minus_one, history, k, k * options.tstep,
                options.tstep, abstol, 1e-50, model)
            z[:, k + 1] = v
            v, vn, vn_minus_one = vn_minus_one, v, vn
            reached = k
        end
    catch e
        e isa Union{DomainError, ErrorException} || rethrow()
        stopped_at = (reached + 1) * options.tstep
        @printf("  the operator stops at t = %g, after %d steps\n", stopped_at, reached)
        write_record(record, merge(options, (; stopped_at));
            t = collect(0:reached) .* options.tstep, z = z[:, 1:(reached + 1)], data...)
        rethrow()
    end
    dist.particles.v[1, :] .= z[:, end]
    return (initial, final = dist, t = collect(0:n) .* options.tstep, z,
        diagnostics = diagnostics(z, dist, sdist))
end

"""
    frozen_rhs!(v̇, v, model, n₀, u₀, ε₀)

The equation of motion of the bracket of `model` at the velocities `v`, with the density, the
mean and the energy moment of the bracket held at `n₀`, `u₀` and `ε₀` in place of the moments of
`v`, written into `v̇`.
"""
function frozen_rhs!(v̇, v, model, n₀, u₀, ε₀)
    cache = model.cache[Float64]
    VlasovMethods._projection!(cache.rhs, cache.values, v, model.dist, cache.sdist)
    VlasovMethods.compute_J!(cache.J, cache.g, cache.values, cache.sdist, model)
    VlasovMethods.compute_dS!(
        cache.dS, cache.values, cache.J, v, cache.sdist, model, model.dist)
    w = view(model.dist.particles.w, 1, :)
    VlasovMethods.bracket_rhs!(v̇, v, w, n₀, u₀, ε₀, cache.dS, sum(cache.dS),
        dot(v, cache.dS), model.bracket)
    return v̇ .*= model.ν
end

"""
    temperature(v, w)

The temperature `Σ w (v − u)² / Σ w` of the velocities `v` with the weights `w`, `u` their mean.
"""
function temperature(v, w)
    n = sum(w)
    u = dot(w, v) / n
    return sum(w .* (v .- u) .^ 2) / n
end

"""
    relax_frozen(initial, bracket, options) -> NamedTuple

The rate check's run: the particles of `initial` dilated about their mean to the temperature
`T₀ = RATE_TEMPERATURE`, then contracted about it by `1 / √2` to `T₀ / 2`, and marched by the
implicit midpoint rule on `frozen_rhs!` of `bracket` with the density and the mean of `initial`
and the temperature `T₀`, for the steps of `options`.

Each step solves `x = vⁿ + h v̇((x + vⁿ) / 2)` by the fixed-point iteration damped by
`RATE_SOLVER_DAMPING`, `x ← x + β r` with the residual `r = vⁿ + h v̇((x + vⁿ) / 2) − x`, until
`‖r‖₂ ≤ RATE_SOLVER_TOLERANCE ‖x‖₂`. Where a particle in a sparse tail sits where `f_s` is small,
the step is stiff and the iteration can stall above that tolerance: after
`RATE_SOLVER_ITERATIONS` iterations the step takes the last iterate, and the run counts the step
and keeps the largest residual of such a step.

Returns the times `t`, the temperature `T` at each time, the density `n₀` of `initial`, the
temperature `T₀`, the count `stalled` of the steps that stopped at the iteration cap and the
largest residual `residual` of those, zero where none did.
"""
function relax_frozen(initial, bracket, options)
    dist = deepcopy(initial)
    w = dist.particles.w[1, :]
    v₀ = dist.particles.v[1, :]
    n₀ = sum(w)
    u₀ = dot(w, v₀) / n₀
    T₀ = RATE_TEMPERATURE
    v₀ = u₀ .+ (v₀ .- u₀) .* sqrt(T₀ / temperature(v₀, w))
    model = MetriplecticLenardBernstein(dist,
        CollisionEntropy(spline_distribution(options.nbasis)); ν = options.nu, bracket)

    n = nsteps(options.final_time, options.tstep)
    h = options.tstep
    v = u₀ .+ (v₀ .- u₀) ./ sqrt(2)
    T = [temperature(v, w)]
    v̇ = similar(v)
    r = similar(v)
    stalled = 0
    residual = 0.0
    for _ in 1:n
        x = copy(v)
        converged = false
        for _ in 1:RATE_SOLVER_ITERATIONS
            frozen_rhs!(v̇, (x .+ v) ./ 2, model, n₀, u₀, T₀ + u₀^2)
            r .= v .+ h .* v̇ .- x
            converged = norm(r) ≤ RATE_SOLVER_TOLERANCE * norm(x)
            converged && break
            x .+= RATE_SOLVER_DAMPING .* r
        end
        if !converged
            stalled += 1
            residual = max(residual, norm(r))
        end
        v = x
        push!(T, temperature(v, w))
    end
    return (t = collect(0:n) .* h, T, n₀, T₀, stalled, residual)
end

"""
    fit_rate(t, T, T₀)

The rate `λ` of `T(t) − T₀ ∝ exp(−λ t)`: minus the least-squares slope of `ln|T(t) − T₀|`
against `t`. Throws an `ArgumentError` where the temperature at the first sample is not positive.
"""
function fit_rate(t, T, T₀)
    first(T) > 0 || throw(ArgumentError(
        "the temperature at the first sample of the rate check is $(first(T)), not positive"))
    y = log.(abs.(T .- T₀))
    t̄ = sum(t) / length(t)
    ȳ = sum(y) / length(y)
    return -sum((t .- t̄) .* (y .- ȳ)) / sum((t .- t̄) .^ 2)
end

"""
    closed_rate(bracket, ν, n₀, T₀)

The closed form of the rate of the temperature under the frozen coefficients: `2ν` for
`DoughertyBracket`, and `2ν n₀ / T₀` for `JeyakumarBracket`, whose equation of motion is the
Dougherty one times `n₀ / T₀`. At the rate check's `T₀ = 2` the second is `ν n₀`.
"""
closed_rate(::DoughertyBracket, ν, n₀, T₀) = 2ν
closed_rate(::JeyakumarBracket, ν, n₀, T₀) = 2ν * n₀ / T₀

other_bracket(::DoughertyBracket) = JeyakumarBracket()
other_bracket(::JeyakumarBracket) = DoughertyBracket()

"""
    rate_check(initial, bracket, options) -> NamedTuple

The rate check of `bracket` from the particles `initial`: the run of `relax_frozen`, the fitted
rate, the closed form of `bracket`, the control, the closed form of the other bracket, the
relative difference of the fit to each, and whether the check passes at `options.tolerance`.
"""
function rate_check(initial, bracket, options)
    run = relax_frozen(initial, bracket, options)
    fitted = fit_rate(run.t, run.T, run.T₀)
    closed = closed_rate(bracket, options.nu, run.n₀, run.T₀)
    control = closed_rate(other_bracket(bracket), options.nu, run.n₀, run.T₀)
    difference = abs(fitted - closed) / closed
    control_difference = abs(fitted - control) / control
    passes = difference ≤ options.tolerance && control_difference > options.tolerance
    return (; run..., fitted, closed, control, difference, control_difference, passes)
end

function write_rate_table(name, bracket_name, rate, options)
    path = results_path(name)
    open(path, "w") do io
        println(io, "# The rate check of the temperature under frozen coefficients")
        println(io)
        println(
            io, "`", COMMAND, "`, commit `", COMMIT, "`, `T₀ = ", rate.T₀, "`, tolerance ",
            options.tolerance, ". ", rate.stalled, " of ", length(rate.t) - 1,
            " steps stop at the iteration cap, the largest residual ", rate.residual, ".")
        println(io)
        println(io, "| bracket | fitted λ | closed form | relative difference | control |",
            " control's relative difference |")
        println(io, "|:--|--:|--:|--:|--:|--:|")
        @printf(io,
            "| %s | %.6f | %.6f | %.3e | %.6f | %.3e |\n", bracket_name, rate.fitted,
            rate.closed, rate.difference, rate.control, rate.control_difference)
    end
    @printf("  wrote %s\n", path)
    return path
end

function main(args)
    case_name = flag_value(args, "--case", "shifted-normal")
    case = choose(CASES, "--case", case_name)
    bracket_name = flag_value(args, "--bracket", "dougherty")
    bracket = choose(BRACKETS, "--bracket", bracket_name)
    quick, options = pd_lb_options(args;
        paper = (; case = "shifted-normal", bracket = "dougherty", case.paper...,
            tolerance = 0.05),
        quick = case.quick)

    @printf("%s: %s run, case %s, bracket %s, %d particles, M = %d, h = %g, t = %g\n",
        NAME,
        quick ? "quick" : "paper", case_name, bracket_name, options.npart, options.nbasis,
        options.tstep, options.final_time)
    stem = "$(NAME)_$(case_name)_$(bracket_name)"

    rate = nothing
    attributes = options
    if case_name == RATE_CASE
        initial = initial_distribution(case.example, options.npart;
            seed = get(options, :seed, nothing))
        rate = rate_check(initial, bracket, options)
        @printf("  rate check: fitted λ = %.6f, closed form %.6f, relative difference %.3e\n",
            rate.fitted, rate.closed, rate.difference)
        @printf("  control: the other form's rate %.6f, relative difference %.3e\n",
            rate.control, rate.control_difference)
        @printf("  %d of %d steps of the rate check stop at the iteration cap, %s %.3e\n",
            rate.stalled, length(rate.t) - 1, "the largest residual", rate.residual)
        write_rate_table(stem * "_rate.md", bracket_name, rate, options)
        attributes = merge(options,
            (rate_fitted = rate.fitted, rate_closed = rate.closed,
                rate_control = rate.control,
                rate_stalled = rate.stalled, rate_residual = rate.residual))
    end

    rate_data = rate === nothing ? (;) :
                (rate_t = rate.t, rate_temperature = rate.T)
    relaxed = march(case, bracket, attributes; record = stem * ".hdf5", data = rate_data)
    d = relaxed.diagnostics
    write_record(stem * ".hdf5", attributes; relaxed.t, relaxed.z, d..., rate_data...)
    save_distribution(stem * "_initial_final_distribution.pdf", relaxed, case.xdomain)
    save_figure(stem * "_conservation.pdf", plot_drift(relaxed.t, (
        d.mass, d.momentum, d.energy)))
    save_figure(stem * "_entropy_log.pdf",
        plot_drift(relaxed.t, (d.entropy,); labels = ("entropy",)))

    report_time()

    if rate !== nothing && !rate.passes
        @printf("%s: the rate check fails at the tolerance %g\n", NAME, options.tolerance)
        exit(1)
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
