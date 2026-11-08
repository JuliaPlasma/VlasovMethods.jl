# The three 2D reference cases of the Landau operator under either of the package's two methods:
# the anisotropic bi-Maxwellian, the double Maxwellian and the BKW solution, run by `PairwiseLandau`
# with its Gonzalez discrete gradient step or by `Landau`, the spline method of the Landau
# manuscript's appendix, with its implicit midpoint step; and the rate check of the pairwise method
# at the Maxwell molecules.
#
#   julia --project=scripts scripts/landau_cases.jl --method pairwise|spline \
#       [--case bimaxwellian|double-maxwellian|bkw] [--gamma -3] [--rate-check] [--quick] \
#       [--seed 20261107] [--npart 1000] [--nsteps 100] [--tstep h] [--ncell 8] \
#       [--tolerance 0.02] [--runs-dir dir] [--results-dir dir]
#
# `--method` has no default. `pairwise` takes `PairwiseLandau` and `discrete_gradient_step!`, with
# the step `h = 1e-5`; `spline` takes `Landau` and `Picard_iterate_Landau_nls!`, with `h = 1e-8`.
# Both solve each step to the absolute residual `ABSTOL`. `--case` picks the initial condition,
# drawn from `Random.MersenneTwister(--seed)` by the type of `src/examples/`:
#
# - `bimaxwellian`, the default: `AnisotropicBiMaxwellian(0.5, 0.3)`;
# - `double-maxwellian`: `DoubleMaxwellian2D` of two halves of the variances `(0.5, 0.5)` and
#   the means `(1, 0)` and `(−1, 0)`;
# - `bkw`: `BKWDistribution(0.2)`, the BKW profile at `t₀ = 0.2`.
#
# Every particle has the weight `1 / N`. `--gamma` is the exponent of the kernel `|u|^{γ+2}`,
# `−3` for the Coulomb case by default, and `ν = 1`. The spline is quadratic, on `--ncell` cells
# per axis of the square `[−L, L]²`, with `L` 1.5 times the largest velocity component of the
# initial particles. A method or a case that is not one of these throws an `ArgumentError` that
# names the flag and the value, before any run and before any directory is created, and a missing
# `--method` one that names the flag.
#
# A run takes `--nsteps` steps of `h` and writes, with the file stem
# `landau_<case>_<method>_gamma<γ>`, the record `<stem>.hdf5` into `--runs-dir`, with the times
# `t`, the velocities `v` (`2 × N` per time), the mass, the momentum (two components per time), the
# energy `Σ w |v|² / 2` and the entropy `S_h` of the projection, and two figures into
# `--results-dir`: `<stem>_conservation.pdf`, the relative drift of the mass and the energy
# (`plot_drift`), and `<stem>_entropy.pdf`, `S_h` against `t` (`plot_entropy`). The momentum of the
# samples is near zero, so its relative drift is not defined; the run prints its largest absolute
# drift instead. Where a step throws — the `DomainError` of the projection or the `ErrorException`
# of a step whose solve does not converge — the run stops: the record is written with the attribute
# `stopped_at`, the time of that step, and the steps reached, no figure is written, and the script
# exits non-zero.
#
# `--rate-check` first runs the rate check of the pairwise method on case (a). At `γ = 0` the
# difference `T₁₁ − T₂₂` of the second moments decays at the closed rate `8` at `ν = 1` and density
# one, whatever the temperatures, which `scripts/reference_forms.jl` derives; the check reads that
# rate from its function `rate_checks`. It takes the weighted lattice of `RATE_LATTICE²` particles on
# `[−RATE_EXTENT, RATE_EXTENT]²`, weighted by the bi-Maxwellian of the variances `0.5` and `0.3` and
# normalised to one, and `RATE_NCELL` quadratic cells per axis on `[−RATE_DOMAIN, RATE_DOMAIN]²`. The
# rate is instantaneous: `λ = −(d gap/dt) / gap` at `t = 0`, with `gap = T₁₁ − T₂₂` and its
# derivative from `collisional_vectorfield!`, so the check needs no time step and no seed. It
# prints `λ`, the closed form and their relative difference, and the control, the closed form with
# the constant `4` of `reference_forms.jl` in place of `8`, with its relative difference to `λ`. The
# script exits non-zero, after it has written its files, where the relative difference exceeds
# `--tolerance`, `2 %` by default, or the control's does not. The table goes to `<stem>_rate.md` and
# the values to attributes of the record.
#
# The rate check runs only with `--method pairwise`, `--gamma 0` and `--case bimaxwellian`; any
# other combination with `--rate-check` throws an `ArgumentError` that names the flag and the
# value, before any directory is created. With `γ ≠ 0` there is no closed rate, and a run reports
# conservation and the entropy only. The spline method's vector field has the opposite sign of the
# appendix's `v̇ = −K⁺LJ`, so where `N ≥ 400` and `M ≥ 289` its instantaneous `λ` lies near `−8`
# where the closed one is `8` (K96 of `KNOWN_ISSUES.md`); `scripts/landau_compare.jl` reports it as
# measured.
#
# `--quick` takes 200 particles, 5 steps and 4 cells per axis, `M = 36`. The rate check keeps its
# own lattice and cells.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

module ReferenceForms
include(joinpath(@__DIR__, "reference_forms.jl"))
end

using LinearAlgebra

const NAME = "landau_cases"

const CASES = Dict(
    "bimaxwellian" => seed -> AnisotropicBiMaxwellian(0.5, 0.3; seed),
    "double-maxwellian" =>
        seed -> DoubleMaxwellian2D(;
            variances = ((0.5, 0.5), (0.5, 0.5)), means = ((1.0, 0.0), (-1.0, 0.0)), seed),
    "bkw" => seed -> BKWDistribution(0.2; seed))

const METHODS = Dict(
    "pairwise" => (model = PairwiseLandau, tstep = 1e-5),
    "spline" => (model = Landau, tstep = 1e-8))

# The order of the B-splines of every Landau run: 3, quadratic.
const LANDAU_ORDER = 3

# The half-width of the velocity domain of a run, as a multiple of the largest velocity component
# of the initial particles.
const DOMAIN_MARGIN = 1.5

# The absolute residual to which each step is solved.
const ABSTOL = 1e-10

# The weighted lattice of the rate check: `RATE_LATTICE` points per axis on
# `[−RATE_EXTENT, RATE_EXTENT]`, `RATE_NCELL` cells per axis on `[−RATE_DOMAIN, RATE_DOMAIN]`, and
# the variances of the bi-Maxwellian that weights it.
const RATE_LATTICE = 60
const RATE_EXTENT = 3.4
const RATE_NCELL = 16
const RATE_DOMAIN = 3.5
const LATTICE_VARIANCES = (0.5, 0.3)

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
    landau_spline(ncell, halfwidth)

The quadratic spline distribution of two velocity dimensions on `ncell` cells per axis of
`[−halfwidth, halfwidth]²`, `(ncell + 2)²` basis functions.
"""
function landau_spline(ncell::Integer, halfwidth)
    return SplineDistribution(1, 2, ncell + 1, LANDAU_ORDER, (-halfwidth, halfwidth))
end

"""
    weighted_lattice(n; extent = RATE_EXTENT, variances = LATTICE_VARIANCES)

The `n × n` lattice of velocities on `[−extent, extent]²`, each particle weighted by the
bi-Maxwellian of `variances`, the weights normalised to one.
"""
function weighted_lattice(n::Integer; extent = RATE_EXTENT, variances = LATTICE_VARIANCES)
    dist = ParticleDistribution(1, 2, n^2)
    axis = range(-extent, extent; length = n)
    for (k, (x, y)) in enumerate(Iterators.product(axis, axis))
        dist.particles.v[:, k] .= (x, y)
        dist.particles.w[1, k] = exp(-x^2 / (2 * variances[1]) - y^2 / (2 * variances[2]))
    end
    dist.particles.w ./= sum(dist.particles.w)
    return dist
end

"""
    temperatures(v, w) -> (T₁₁, T₂₂)

The two diagonal second moments `Σ w (v_c − u_c)² / Σ w` of the velocities `v` (`2 × N`) with the
weights `w` about their mean `u`.
"""
function temperatures(v, w)
    c = v .- (v * w) ./ sum(w)
    return (sum(w .* c[1, :] .^ 2), sum(w .* c[2, :] .^ 2)) ./ sum(w)
end

"""
    gap_rate(model, dist) -> (rate, T₁₁, T₂₂)

The instantaneous `(d gap/dt) / gap` of the particles of `dist` under `model`, with
`gap = T₁₁ − T₂₂` and its derivative `2 Σ w (c₁ v̇₁ − c₂ v̇₂) / Σ w` for the velocities `c` about
the mean and `v̇` from `collisional_vectorfield!`, and the two temperatures.
"""
function gap_rate(model, dist)
    v = copy(dist.particles.v)
    v̇ = similar(v)
    w = vec(dist.particles.w)
    VlasovMethods.collisional_vectorfield!(v̇, v, nothing, model)
    c = v .- (v * w) ./ sum(w)
    T₁₁, T₂₂ = temperatures(v, w)
    derivative = 2 * sum(w .* (c[1, :] .* v̇[1, :] .- c[2, :] .* v̇[2, :])) / sum(w)
    return derivative / (T₁₁ - T₂₂), T₁₁, T₂₂
end

"""
    closed_gap_rate(T₁₁, T₂₂)

The closed form of `(d gap/dt) / gap` at `γ = 0`, `ν = 1` and density one: the rate of `T₁₁ − T₂₂`
that `rate_checks` of `scripts/reference_forms.jl` integrates from the moment equations, divided by
`T₁₁ − T₂₂`, which is `−8`.
"""
closed_gap_rate(T₁₁, T₂₂) = ReferenceForms.rate_checks(T₁₁, T₂₂).gap / (T₁₁ - T₂₂)

"""
    control_gap_rate(T₁₁, T₂₂)

The control of the rate check: the closed form of `control_rows` of `scripts/reference_forms.jl`,
with the constant `4` in place of `8`, divided by `T₁₁ − T₂₂`, which is `−4`.
"""
function control_gap_rate(T₁₁, T₂₂)
    return only(ReferenceForms.control_rows(T₁₁, T₂₂)).reference / (T₁₁ - T₂₂)
end

"""
    rate_check(tolerance) -> NamedTuple

The rate check of `PairwiseLandau` at `γ = 0` on the weighted lattice of `RATE_LATTICE²`
particles: `λ = −(d gap/dt) / gap`, the closed form `8`, the control `4`, the relative difference
of `λ` to each, and whether the check passes at `tolerance`.
"""
function rate_check(tolerance)
    dist = weighted_lattice(RATE_LATTICE)
    model = PairwiseLandau(dist,
        CollisionEntropy(landau_spline(RATE_NCELL, RATE_DOMAIN)); γ = 0, ν = 1.0)
    measured, T₁₁, T₂₂ = gap_rate(model, dist)
    λ = -measured
    closed = -closed_gap_rate(T₁₁, T₂₂)
    control = -control_gap_rate(T₁₁, T₂₂)
    difference = abs(λ - closed) / closed
    control_difference = abs(λ - control) / control
    passes = difference ≤ tolerance && control_difference > tolerance
    return (; λ, closed, control, difference, control_difference, passes, T₁₁, T₂₂)
end

function write_rate_table(name, rate, tolerance)
    path = results_path(name)
    open(path, "w") do io
        println(io, "# The rate check of the pairwise method at the Maxwell molecules")
        println(io)
        println(io, "`", COMMAND, "`, commit `", COMMIT, "`: `", RATE_LATTICE, " × ",
            RATE_LATTICE, "` lattice on `[−", RATE_EXTENT, ", ", RATE_EXTENT, "]²`, `",
            RATE_NCELL, " × ", RATE_NCELL, "` quadratic cells on `[−", RATE_DOMAIN, ", ",
            RATE_DOMAIN, "]²`, `γ = 0`, `ν = 1`, tolerance ", tolerance, ".")
        println(io)
        println(io, "| method | `T₁₁` | `T₂₂` | `λ = −(d gap/dt)/gap` | closed form |",
            " relative difference | control | control's relative difference |")
        println(io, "|:--|--:|--:|--:|--:|--:|--:|--:|")
        @printf(io, "| pairwise | %.6f | %.6f | %.6f | %.6f | %.3e | %.6f | %.3e |\n",
            rate.T₁₁, rate.T₂₂, rate.λ, rate.closed, rate.difference,
            rate.control,
            rate.control_difference)
    end
    @printf("  wrote %s\n", path)
    return path
end

"""
    initial_particles(case, npart, seed)

`npart` particles of two velocity dimensions drawn from the example of `case` with `seed`, each
of the weight `1 / npart`.
"""
function initial_particles(case, npart::Integer, seed::Integer)
    dist = ParticleDistribution(1, 2, npart)
    dist.particles.w .= 1 / npart
    return initialize!(dist, case(seed))
end

"""
    stepper(method_name, model, dist, tstep)

The function `k -> …` that takes the step `k` of `model` from the velocities of `dist` and leaves
the new velocities there: `discrete_gradient_step!` for `pairwise`, `Picard_iterate_Landau_nls!`
for `spline`, each to the residual `ABSTOL`.
"""
function stepper(method_name, model, dist, tstep)
    method_name == "pairwise" &&
        return k -> discrete_gradient_step!(model, tstep; abstol = ABSTOL)
    v = copy(dist.particles.v)
    v_prev_2 = copy(v)
    rhs = zeros(size(v)..., 2)
    return function (k)
        VlasovMethods.Picard_iterate_Landau_nls!(
            model, tstep, k, k * tstep, v, v_prev_2, rhs;
            abstol = ABSTOL)
        v_prev_2 .= v
        v .= dist.particles.v
    end
end

"""
    landau_diagnostics(v, dist, sdist) -> NamedTuple

The mass `Σ w`, the momentum `Σ w v` (two components per time), the energy `Σ w |v|² / 2` and the
entropy `S_h = ∫ f_s log|f_s| dv` of the projection onto `sdist`, at each time of `v`
(`2 × N × times`). `dist` gives the weights; the projection overwrites the coefficients of `sdist`.
"""
function landau_diagnostics(v, dist, sdist)
    w = vec(dist.particles.w)
    entropy = CollisionEntropy(sdist)
    times = axes(v, 3)
    S = map(times) do k
        projection(v[:, :, k], dist, sdist)
        entropy()
    end
    return (mass = fill(sum(w), length(times)),
        momentum = reduce(hcat, [v[:, :, k] * w for k in times]),
        energy = [sum(w .* vec(sum(abs2, v[:, :, k]; dims = 1))) / 2 for k in times],
        entropy = S)
end

"""
    march(case, method_name, options; record, attributes) -> NamedTuple

The run of `method_name` from `options.npart` particles of `case` with `options.seed`, at the
exponent `options.gamma`, on `options.ncell` quadratic cells per axis, for `options.nsteps` steps
of `options.tstep`. Returns the times `t`, the velocities `v`, the `diagnostics` and the
`halfwidth` of the velocity domain.

Where a step throws a `DomainError` or an `ErrorException`, the record `record` is written with
`attributes`, the attribute `stopped_at`, the time of that step, and the times and the velocities
of the steps reached, and the error is rethrown.
"""
function march(case, method_name, options; record, attributes)
    dist = initial_particles(case, options.npart, options.seed)
    halfwidth = DOMAIN_MARGIN * maximum(abs, dist.particles.v)
    model = METHODS[method_name].model(dist,
        CollisionEntropy(landau_spline(options.ncell, halfwidth)); γ = options.gamma,
        ν = 1.0)
    step! = stepper(method_name, model, dist, options.tstep)

    n = options.nsteps
    v = zeros(2, options.npart, n + 1)
    v[:, :, 1] = dist.particles.v
    reached = 0
    try
        for k in 1:n
            step!(k)
            v[:, :, k + 1] = dist.particles.v
            reached = k
        end
    catch e
        e isa Union{DomainError, ErrorException} || rethrow()
        stopped_at = (reached + 1) * options.tstep
        @printf("  the operator stops at t = %g, after %d steps\n", stopped_at, reached)
        write_record(record, merge(attributes, (; stopped_at, halfwidth));
            t = collect(0:reached) .* options.tstep, v = v[:, :, 1:(reached + 1)])
        rethrow()
    end
    diagnostics = landau_diagnostics(v, dist, landau_spline(options.ncell, halfwidth))
    return (t = collect(0:n) .* options.tstep, v, diagnostics, halfwidth)
end

"""
    landau_options(args) -> (method_name, case_name, rate, quick, options)

The method, the case and the options of a run. A missing `--method`, a method or a case that is
not one of `METHODS` and `CASES`, `--rate-check` or `--quick` given twice, and `--rate-check` with
another method than `pairwise`, another exponent than `0` or another case than `bimaxwellian`
each throw an `ArgumentError` that names the flag, before `parse_arguments` creates any directory.
"""
function landau_options(args)
    method_name = flag_value(args, "--method", nothing)
    method_name === nothing && throw(ArgumentError(
        "`--method` is required; this script takes " *
        join(sort!(collect(keys(METHODS))), ", ")))
    method = choose(METHODS, "--method", method_name)
    case_name = flag_value(args, "--case", "bimaxwellian")
    choose(CASES, "--case", case_name)

    n = count(==("--rate-check"), args)
    n ≤ 1 || throw(ArgumentError("`--rate-check` is given twice"))
    rate = n == 1
    if rate
        method_name == "pairwise" || throw(ArgumentError(
            "`--rate-check` is the pairwise method's, and `--method` is `$(method_name)`: the " *
            "spline method's vector field has the opposite sign of the appendix's (K96)"))
        gamma = tryparse(Float64, flag_value(args, "--gamma", "-3"))
        gamma == 0 || throw(ArgumentError(
            "`--rate-check` needs the closed rate of the Maxwell molecules, `--gamma 0`, and " *
            "the exponent is `$(flag_value(args, "--gamma", "-3"))`"))
        case_name == "bimaxwellian" || throw(ArgumentError(
            "`--rate-check` is the check of `--case bimaxwellian`, and `--case` is " *
            "`$(case_name)`"))
    end

    quick, options = pd_lb_options(filter(!=("--rate-check"), args);
        paper = (
            method = method_name, case = "bimaxwellian", gamma = -3.0, seed = 20261107,
            npart = 1000, nsteps = 100, tstep = method.tstep, ncell = 8, tolerance = 0.02),
        quick = (npart = 200, nsteps = 5, ncell = 4))
    return method_name, case_name, rate, quick, options
end

function main(args)
    method_name, case_name, rate_requested, quick, options = landau_options(args)
    @printf("%s: %s run, method %s, case %s, γ = %g, %d particles, %d cells, h = %g, %d steps\n",
        NAME, quick ? "quick" : "full", method_name, case_name, options.gamma, options.npart,
        options.ncell, options.tstep, options.nsteps)
    stem = @sprintf("landau_%s_%s_gamma%g", case_name, method_name, options.gamma)

    rate = nothing
    attributes = merge(options, (; rate_check = Int(rate_requested)))
    if rate_requested
        rate = rate_check(options.tolerance)
        @printf("  rate check: λ = −(d gap/dt)/gap = %.6f, closed form %.6f, %s %.3e\n",
            rate.λ, rate.closed, "relative difference", rate.difference)
        @printf("  control: the closed form with the constant 4, %.6f, %s %.3e\n",
            rate.control, "relative difference", rate.control_difference)
        write_rate_table(stem * "_rate.md", rate, options.tolerance)
        attributes = merge(attributes,
            (rate_measured = rate.λ, rate_closed = rate.closed,
                rate_control = rate.control))
    end

    run = march(CASES[case_name], method_name, options; record = stem * ".hdf5", attributes)
    d = run.diagnostics
    write_record(stem * ".hdf5", merge(attributes, (; run.halfwidth)); run.t, run.v, d...)
    @printf("  relative energy drift %.3e, largest momentum drift %.3e, S_h from %.6f to %.6f\n",
        abs(d.energy[end] - d.energy[1]) / d.energy[1],
        maximum(abs, d.momentum .- d.momentum[:, 1]), d.entropy[1], d.entropy[end])
    save_figure(stem * "_conservation.pdf",
        plot_drift(run.t, (d.mass, d.energy); labels = ("mass", "energy")))
    save_figure(stem * "_entropy.pdf", plot_entropy(run.t, d.entropy))

    report_time()

    if rate !== nothing && !rate.passes
        @printf("%s: the rate check fails at the tolerance %g\n", NAME, options.tolerance)
        exit(1)
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
