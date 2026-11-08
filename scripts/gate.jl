# The short run of the figure drivers of the collision manuscripts: each driver with `--quick`, in
# a process of its own, into a temporary directory.
#
#   julia --startup-file=no --project=scripts scripts/gate.jl
#
# A driver passes when it exits 0 within `BUDGET` seconds and writes each of its figures and each
# of its HDF5 files, and each HDF5 file carries
#
# - the attribute `commit`, equal to `git rev-parse HEAD` of this checkout,
# - the attribute `command`, equal to the arguments the gate passed to the driver as
#   `Base.shell_escape` quotes them,
# - the values of `--quick` in its option attributes, so that a driver that ignores `--quick`
#   fails here even where its full run fits the budget.
#
# Every driver that runs to its end must print its wall-clock time. `pd_lb_fig11.jl` must also
# print the relative error of the fourth cumulant at each of the first three steps, each below
# `5e-3`, where the quick run measures `4.4e-4` to `1.5e-3`, and its first step must solve the
# implicit midpoint rule to a relative residual below `1e-10`, where it measures `7e-15`.
# `pd_lb_fig7_8.jl` must start from the velocities that `DoubleMaxwellian` draws after
# `Random.seed!(1)`, and the quick records of `pd_lb_fig1_2.jl` must hold no `stopped_at`, because
# neither quick sweep stops.
#
# Six more runs check the edges of the drivers:
#
# - `pd_lb_fig11.jl` in an empty working directory, with `--runs-dir` and `--results-dir` two
#   levels below a directory that does not have them, the first level named with a space, and
#   with `GIT_DIR` naming a path that does not exist: both directories are created and hold its
#   files, the working directory stays empty, `commit` is `"unknown"`, and `command` quotes the
#   paths with the space;
# - `pd_lb_fig1_2.jl` with the sweeps `M ∈ {24, 25, 26}` at `10⁴` particles and
#   `N ∈ {1000, 5623, 10000, 31623}`, where the projection is not positive at a particle of
#   `M = 25`, of `N = 5623` and of `N = 31623`: it exits 0, prints each stop, records
#   `stopped_at = 25` beside `M ∈ {24, 26}` and `stopped_at = 5623` beside `N ∈ {1000, 10000}`,
#   and writes both figures;
# - `pd_lb_fig1_2.jl` with `M ∈ {24, 3, 41}`, where `M = 3` throws an `ArgumentError`, which is
#   not the model's: it exits non-zero with that error and writes no figure 1;
# - `pd_lb_fig11.jl` with `--tolerance 1e-4`, below the relative errors of the quick run, which
#   are above `4e-4`: it prints the three relative errors and exits non-zero;
# - `pd_lb_fig7_8.jl` with 60 particles, the step `0.03` and `--seed 4`, where the projection is
#   not positive at a particle in the fifth step: it exits non-zero with the `DomainError`, writes
#   no figure, and writes `pd_lb_fig7_8.hdf5` with `commit`, `command`, `stopped_at = 5 * 0.03`,
#   the four steps reached in `t` and `z`, and the velocities of the seed 4 in the first column;
# - `pd_lb_fig11.jl` with `--runs-dir` below a regular file, which cannot be created: it exits
#   non-zero with an error that names the path given.
#
# The quick record of `pd_lb_fig11.jl` must hold the tolerance `0.05`.
#
# `lenard_bernstein_metriplectic.jl` runs with `--quick` six times, once for each of its three
# cases and two brackets, each in a process of its own. Each run must exit 0, print its wall-clock
# time, write its three figures and its record with the quick parameters, its case and its bracket,
# and no `stopped_at`. The two runs of the shifted normal also print the rate check, which relaxes
# from the temperature `1` to `T₀ = 2` at the density `n₀ = 1` of the weights `1 / N`: the gate
# takes the closed form `2ν` for the Dougherty bracket and `2ν n₀ / T₀ = ν` for the Jeyakumar
# bracket, from the bracket that it passed, and checks that the printed closed form and control are
# these two, that the record's `rate_temperature` starts at `1`, that the fitted rate lies within a
# relative `5 %` of the closed form, and that it lies more than `5 %` from the control. The record
# of the double Maxwellian with the Dougherty bracket must start from the velocities that
# `DoubleMaxwellian` draws after `Random.seed!(1)`. Six more runs check its edges:
#
# - `--case maxwellian` and `--bracket landau`: each exits non-zero with the `ArgumentError` that
#   names the flag and the value, and creates neither `--runs-dir` nor `--results-dir`;
# - `--case uniform` with `--tolerance`, an option of the shifted normal alone: it exits non-zero
#   with the `ArgumentError` that names `--tolerance`, and creates neither directory;
# - the shifted normal with the Dougherty bracket and 1000 particles, where the Picard solve of the
#   sixth step does not converge: it exits non-zero, writes no figure, and writes its record with
#   `commit`, `command`, `stopped_at = 6 * 8e-4`, the five steps reached in `t` and `z`, and the 64
#   samples of its rate check in `rate_t` and `rate_temperature`;
# - the shifted normal with `--tolerance 1e-3`, below the relative difference of the quick fit to
#   its closed form, which is `1.5e-2`, and with `--tolerance 1.5`, above the relative difference
#   of the fit to its control, which is `0.97`: each prints its wall-clock time, writes its three
#   figures and its record, says that the rate check fails and exits non-zero.
#
# The gate also includes the driver in a module of its own and checks that `fit_rate` throws an
# `ArgumentError` where a temperature reaches or crosses `T₀`, that it fits the rate of an exact
# exponential, that `flag_value` takes the default where the flag is the last argument, and that
# `march` of the 200 particles of the quick shifted normal over three steps returns the initial
# particles at the first column of its velocities and the final ones at the last.
#
# Beside the runs, the gate includes `pd_lb_common.jl` and checks `--quick`, the step count and the
# spline basis directly.
#
# `landau_cases.jl` runs with `--quick` four times at its default `γ = −3`, without a rate check:
# case (a), the bi-Maxwellian, under both methods, the double Maxwellian under `pairwise` and the
# BKW profile under `spline`. Each run must exit 0, print its wall-clock time, write its two
# figures and its record with the quick parameters, its method, its case, `gamma = −3` and
# `rate_check = 0`, and print no rate check. The rate check runs with `--method pairwise --gamma 0
# --rate-check`: the gate checks that it prints the closed form `8` and the control `4`, that the
# instantaneous rate lies within `2 %` of `8`, where the quick run measures `8.009493`, and more
# than `2 %` from `4`, and that its record carries the measured rate. Two more runs check its
# edges:
#
# - `--method pairwise --tstep 1e-3 --ncell 8`, whose third step does not converge: it exits
#   non-zero, writes no figure, and writes its record with `commit`, `command`, `ncell = 8`,
#   `stopped_at = 3e-3` and the two steps reached. On the 4 cells of `--quick` the same steps
#   converge;
# - the rate check with `--tolerance 1e-4`, below the relative difference `1.2e-3` of the rate: it
#   writes its files, says that the rate check fails and exits non-zero.
#
# The refusals are checked on `landau_options`, which `main` calls first, so that seven processes
# do not load the package only to throw: the gate includes `landau_cases.jl` in a module of its
# own, and `landau_options` must throw the `ArgumentError` that names the flag and the value, and
# create neither `--runs-dir` nor `--results-dir`, for `--rate-check` at the default `γ = −3`,
# with `--method spline`, with `--case bkw` and given twice, for the unknown values
# `--method splines` and `--case maxwellian`, and for no `--method`. `rate_check` must fail at the
# tolerance `1.5`, which its relative difference meets and the control's, `1.0`, does not.
#
# `landau_compare.jl` runs with `--quick`: it must exit 0, print its wall-clock time, write its
# table, its two figures and its record with the quick sweep, and its table must hold one row per
# method, lattice and cell count, each with `dN = 2N` and the verdict of `dN ≥ M − 1` that the gate
# computes from `N` and `M` itself — `N = 144`, `M = 289` holds at the bound and `N = 144`,
# `M = 324` fails —, a particle in every support, the closed form `−8` and the relative error of
# the measured value against it. The measured value carries its sign: at `N = 400`, `M = 324` the
# pairwise method's is negative and the spline method's positive (K96). No word of the output of
# `landau_cases.jl` and `landau_compare.jl`, of their tables, of the Landau section of
# `docs/src/collisions.md` or of the CHANGELOG entries that name them ranks one method above the
# other: the gate fails on `faster`, `better`, `slower`, `cheaper`, `outperforms` and `wins` in
# any case.
#
# The drivers run at once, so the gate takes about as long as the slowest of them. `BUDGET` is
# 400 s, measured from the gate's start, because whole gate runs take up to 229 s and a loaded
# machine slows them. The gate prints the size of each file and its own wall-clock time,
# and exits 1 when a check fails. A driver still running when the gate exits is killed.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

module MPLB
include(joinpath(@__DIR__, "lenard_bernstein_metriplectic.jl"))
end

module LandauCases
include(joinpath(@__DIR__, "landau_cases.jl"))
end

const ROOT = dirname(@__DIR__)
const BUDGET = 400.0
const MIDPOINT_RESIDUAL = 1e-10
const EXPECTED_COMMIT = readchomp(`git -C $(ROOT) rev-parse HEAD`)

const LANDAU_SCRIPT = "landau_cases.jl"
const LANDAU_QUICK = (npart = 200, nsteps = 5, ncell = 4)
const LANDAU_RUNS = (("pairwise", "bimaxwellian"), ("spline", "bimaxwellian"),
    ("pairwise", "double-maxwellian"), ("spline", "bkw"))
const LANDAU_FIGURES = ("conservation", "entropy")
const LANDAU_RATE = 8.0
const LANDAU_CONTROL = 4.0
const LANDAU_RATE_TOLERANCE = 0.02
const LANDAU_RATE_ROW = r"^  rate check: λ = −\(d gap/dt\)/gap = (\S+), closed form (\S+), relative difference \S+$"m
const LANDAU_CONTROL_ROW = r"^  control: the closed form with the constant 4, (\S+), relative difference \S+$"m
# The arguments that `landau_options` refuses, and the start of the message of each refusal.
const LANDAU_REFUSED = [
    "--quick --method pairwise --rate-check" =>
        "`--rate-check` needs the closed rate of the Maxwell molecules, `--gamma 0`, and the " *
        "exponent is `-3`",
    "--quick --method spline --gamma 0 --rate-check" => "`--rate-check` is the pairwise method's, and `--method` is `spline`",
    "--quick --method pairwise --gamma 0 --case bkw --rate-check" => "`--rate-check` is the check of `--case bimaxwellian`, and `--case` is `bkw`",
    "--quick --method pairwise --gamma 0 --rate-check --rate-check" => "`--rate-check` is given twice",
    "--quick --method splines" => "unknown --method `splines`",
    "--quick --method pairwise --case maxwellian" => "unknown --case `maxwellian`",
    "--quick" => "`--method` is required"]
const COMPARE_SCRIPT = "landau_compare.jl"
const COMPARE_QUICK = (lattice = [12, 20], ncell = [8, 15, 16], repeats = 3)
const COMPARE_ROW = r"^\| (pairwise|spline) \| (\d+) \| (\d+) \| (\d+) \| (holds|fails) \| ([^|]+) \| (\S+) \| (\S+) \| (\S+) \| (\S+) \|$"
const RANKING = r"\b(faster|better|slower|cheaper|outperforms|wins)\b"i
const LANDAU_DOCS_HEADING = "## Particle methods for the Landau operator"

const QUICK = Dict(
    "pd_lb_fig1_2.jl" => (figures = ["pd_lb_fig1.pdf", "pd_lb_fig2.pdf"],
        records = ["pd_lb_fig1.hdf5" => (npart = 2000, nbasis = [24, 41, 78]),
            "pd_lb_fig2.hdf5" => (sweep_npart = [100, 178, 316, 562, 1000, 1778],)]),
    "pd_lb_fig3_6.jl" => (
        figures = ["pd_lb_fig3.pdf", "pd_lb_fig4.pdf", "pd_lb_fig5.pdf", "pd_lb_fig6.pdf"],
        records = ["pd_lb_fig3_5.hdf5" => (npart = 200, final_time = 0.05),
            "pd_lb_fig6.hdf5" => (npart = 200, final_time = 0.05)]),
    "pd_lb_fig7_8.jl" => (figures = ["pd_lb_fig7.pdf", "pd_lb_fig8.pdf"],
        records = ["pd_lb_fig7_8.hdf5" => (npart = 200, final_time = 0.2)]),
    "pd_lb_fig9_10.jl" => (figures = ["pd_lb_fig9.pdf", "pd_lb_fig10.pdf"],
        records = ["pd_lb_fig9_10.hdf5" => (tstep = 1e-3, final_time = 0.05)]),
    "pd_lb_fig11.jl" => (figures = ["pd_lb_fig11.pdf"],
        records = ["pd_lb_fig11.hdf5" =>
            (npart = 200, tstep = 1e-2, final_time = 0.1,
                tolerance = 0.05)]))

const CUMULANT_ROW = r"^  κ₄ at t = (\S+): run \S+, closed form \S+, relative error (\S+)$"

const MPLB_SCRIPT = "lenard_bernstein_metriplectic.jl"
const MPLB_QUICK = Dict(
    "shifted-normal" => (npart = 200, final_time = 0.05),
    "double-maxwellian" => (npart = 1000, final_time = 0.05),
    "uniform" => (npart = 200, tstep = 1e-3, final_time = 0.05))
const MPLB_FIGURES = ("initial_final_distribution", "conservation", "entropy_log")
const MPLB_RATE_TOLERANCE = 0.05
const RATE_ROW = r"^  rate check: fitted λ = (\S+), closed form (\S+), relative difference \S+$"m
const CONTROL_ROW = r"^  control: the other form's rate (\S+), relative difference \S+$"m

const FAILURES = String[]

function check(condition::Bool, message)
    condition || push!(FAILURES, message)
    @printf("  [%s] %s\n", condition ? "pass" : "FAIL", message)
    return condition
end

# A driver started in its own process, its output going to `log`. `out` holds `runs` and `results`.
function start(script, args; out = mktempdir(), dir = ROOT, env = nothing)
    args = [args; "--runs-dir"; joinpath(out, "runs"); "--results-dir";
            joinpath(out, "results")]
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(joinpath(ROOT, "scripts"))
        $(joinpath(ROOT, "scripts", script)) $args`
    env === nothing || (command = addenv(command, env))
    log = tempname() * ".log"
    stream = open(log, "w")
    process = run(pipeline(Cmd(command; dir); stdout = stream, stderr = stream); wait = false)
    return (; script, args, process, log, stream, runs = joinpath(out, "runs"),
        results = joinpath(out, "results"))
end

# Waits for a driver until `BUDGET` seconds after the gate started, kills it after that, and
# checks its exit status, 0 or not 0 where `fails`, and, where `timed`, that it printed its
# wall-clock time. Returns its output.
function finish(job; fails = false, timed = true)
    remaining = BUDGET - (time() - STARTED)
    status = timedwait(() -> process_exited(job.process), max(remaining, 0.0); pollint = 1.0)
    if status === :timed_out
        kill(job.process, Base.SIGKILL)
        wait(job.process)
    end
    close(job.stream)
    output = read(job.log, String)
    finished = check(status === :ok, "$(job.script) finishes within $(BUDGET) s")
    exited = check(finished && success(job.process) != fails,
        "$(job.script) exits $(fails ? "non-zero" : "0")")
    timed &&
        check(occursin(r"wall-clock time \d", output), "$(job.script) prints its wall-clock time")
    exited || print(output)
    return output
end

function check_record(path, args, parameters)
    check(isfile(path), "$(basename(path)) is written") || return nothing
    @printf("    %s: %d bytes\n", basename(path), filesize(path))
    h5open(path, "r") do file
        attributes = HDF5.attrs(file)
        check(haskey(attributes, "commit") && attributes["commit"] == EXPECTED_COMMIT,
            "$(basename(path)) carries commit = $(EXPECTED_COMMIT)")
        command = Base.shell_escape(args...)
        check(haskey(attributes, "command") && attributes["command"] == command,
            "$(basename(path)) carries command = $(command)")
        for (key, value) in pairs(parameters)
            name = String(key)
            check(haskey(attributes, name) && attributes[name] == value,
                "$(basename(path)) records the quick $(name) = $(value)")
        end
    end
    return nothing
end

function check_figure(path)
    check(isfile(path) && filesize(path) > 0, "$(basename(path)) is written") &&
        @printf("    %s: %d bytes\n", basename(path), filesize(path))
    return nothing
end

# The times and the relative errors of the rows of the cumulant comparison in `output`.
function cumulant_rows(output)
    rows = filter(!isnothing, map(line -> match(CUMULANT_ROW, line), split(output, '\n')))
    times = [parse(Float64, row[1]) for row in rows]
    errors = [parse(Float64, row[2]) for row in rows]
    return times, errors
end

# The residual of the implicit midpoint rule `z₁ = z₀ + h v̇((z₀ + z₁) / 2)` over the first step
# of the record of `pd_lb_fig11.jl` at `path`, with `v̇` the right-hand side of `LenardBernstein`
# for the particles of the record on 41 basis functions, relative to the largest change
# `|z₁ − z₀|`.
function midpoint_residual(path, tstep)
    z = h5read(path, "z")
    dist = initialize!(ParticleDistribution(1, 1, size(z, 1)),
        DoubleMaxwellian(; domain = (-5.0, 5.0), shift = 2.0))
    model = LenardBernstein(dist, CollisionEntropy(spline_distribution(41)))
    problem = GeometricIntegrator(model, (0.0, tstep), tstep).equation
    v̇ = similar(z[:, 1])
    problem.equation.v(v̇, 0.0, (z[:, 1] .+ z[:, 2]) ./ 2, problem.parameters)
    return maximum(abs, z[:, 2] .- z[:, 1] .- tstep .* v̇) / maximum(abs, z[:, 2] .- z[:, 1])
end

# The Landau section of `docs/src/collisions.md`: from `LANDAU_DOCS_HEADING` to the next heading
# of its level, or to the end of the page; empty where the page has no such heading.
function landau_docs_section()
    lines = readlines(joinpath(ROOT, "docs", "src", "collisions.md"))
    from = findfirst(==(LANDAU_DOCS_HEADING), lines)
    from === nothing && return ""
    to = findnext(line -> startswith(line, "## "), lines, from + 1)
    return join(lines[from:(to === nothing ? end : to - 1)], "\n")
end

# The entries of `CHANGELOG.md` that name `landau_cases.jl` or `landau_compare.jl`: each entry is a
# bullet at the start of a line, with the lines that continue it.
function landau_changelog_entries()
    entries = String[]
    for line in readlines(joinpath(ROOT, "CHANGELOG.md"))
        if startswith(line, "- ")
            push!(entries, line)
        elseif startswith(line, "  ") && !isempty(entries)
            entries[end] *= "\n" * line
        end
    end
    return filter(e -> occursin(r"landau_(cases|compare)\.jl", e), entries)
end

jobs = [start(script, ["--quick"]) for script in sort!(collect(keys(QUICK)))]

const WORK = mktempdir()
const NESTED = joinpath(mktempdir(), "new dir", "deeper")
stray = start("pd_lb_fig11.jl", ["--quick"]; out = NESTED, dir = WORK,
    env = "GIT_DIR" => joinpath(WORK, "no-such-git-dir"))
stop = start("pd_lb_fig1_2.jl",
    ["--quick", "--npart", "10000", "--nbasis", "24,25,26", "--sweep-npart",
        "1000,5623,10000,31623"])
other = start("pd_lb_fig1_2.jl", ["--quick", "--nbasis", "24,3,41"])
miss = start("pd_lb_fig11.jl", ["--quick", "--tolerance", "1e-4"])
halt = start("pd_lb_fig7_8.jl",
    ["--quick", "--npart", "60", "--tstep", "0.03", "--seed", "4"])
const NOT_A_DIRECTORY = joinpath(mktempdir(), "file")
touch(NOT_A_DIRECTORY)
unmade = start("pd_lb_fig11.jl", ["--quick"]; out = joinpath(NOT_A_DIRECTORY, "sub"))
mplb = [start(MPLB_SCRIPT, ["--quick", "--case", case, "--bracket", bracket])
        for case in sort!(collect(keys(MPLB_QUICK)))
        for bracket in ("dougherty", "jeyakumar")]
unknown_case = start(MPLB_SCRIPT, ["--quick", "--case", "maxwellian"])
unknown_bracket = start(MPLB_SCRIPT, ["--quick", "--bracket", "landau"])
untolerant = start(MPLB_SCRIPT, ["--quick", "--case", "uniform", "--tolerance", "0.05"])
unsolved = start(MPLB_SCRIPT, ["--quick", "--bracket", "dougherty", "--npart", "1000"])
tight = start(MPLB_SCRIPT, ["--quick", "--tolerance", "1e-3"])
loose = start(MPLB_SCRIPT, ["--quick", "--tolerance", "1.5"])
landau = [start(LANDAU_SCRIPT, ["--quick", "--method", method, "--case", case])
          for (method, case) in LANDAU_RUNS]
landau_rate = start(LANDAU_SCRIPT,
    ["--quick", "--method", "pairwise", "--gamma", "0", "--rate-check"])
landau_stop = start(LANDAU_SCRIPT, split("--quick --method pairwise --tstep 1e-3 --ncell 8"))
landau_tight = start(LANDAU_SCRIPT,
    split("--quick --method pairwise --gamma 0 --rate-check --tolerance 1e-4"))
compare = start(COMPARE_SCRIPT, ["--quick"])

# A driver still running when the gate exits, after an error of the gate itself, is killed.
atexit(() -> foreach(
    job -> process_running(job.process) && kill(job.process, Base.SIGKILL),
    [jobs; mplb; landau;
     [stray, stop, other, miss, halt, unmade, unknown_case, unknown_bracket,
         untolerant, unsolved, tight, loose, landau_rate, landau_stop, landau_tight, compare]]))

println("the helpers of pd_lb_common.jl")
let dirs = ["--runs-dir", mktempdir(), "--results-dir", mktempdir()],
    paper = (npart = 1000, tstep = 8e-4), quick = (npart = 200,)

    check(
        try
            pd_lb_options(["--quick", "--quick"]; paper, quick)
            false
        catch e
            e isa ArgumentError
        end,
        "`--quick` given twice throws an ArgumentError")
    isquick, options = pd_lb_options(dirs; paper, quick)
    check(!isquick && options.npart == 1000 && options.tstep == 8e-4,
        "without `--quick` the options are the manuscript's")
    isquick, options = pd_lb_options(["--quick"; dirs]; paper, quick)
    check(isquick && options.npart == 200 && options.tstep == 8e-4,
        "`--quick` replaces the options it names and keeps the others")
end
# `0.07 / 5e-3` is `14.000000000000002` in floating point, which the rounding takes to 14 steps
check(
    nsteps(1.0, 8e-4) == 1250 && nsteps(5.0, 8e-4) == 6250 && nsteps(1.0, 1e-4) == 10_000 &&
        nsteps(1.0, 5e-3) == 200 && nsteps(0.05, 8e-4) == 63 && nsteps(0.1, 1e-2) == 10 &&
        nsteps(0.07, 5e-3) == 14,
    "nsteps gives 1250, 6250, 10000, 200, 63, 10 and 14 steps")
check(
    nbasis(spline_distribution(41).basis) == 41 &&
        ncells(spline_distribution(41).basis) == 38,
    "41 cubic basis functions lie on 38 cells")

println("the rate fit of lenard_bernstein_metriplectic.jl")
for (temperatures, where) in (([1.0, 2.0, 1.5], "reaches"), ([1.0, 2.5, 3.0], "crosses"))
    check(
        try
            MPLB.fit_rate([0.0, 1.0, 2.0], temperatures, 2.0)
            false
        catch e
            e isa ArgumentError
        end,
        "fit_rate throws an ArgumentError where the temperature $(where) T₀: $(temperatures)")
end
let t = collect(0:0.1:1), fitted = MPLB.fit_rate(t, 2.0 .- exp.(-3 .* t), 2.0)
    check(isapprox(fitted, 3.0; rtol = 1e-12),
        "fit_rate gives the rate 3 of T₀ − exp(−3 t): $(fitted)")
end
check(
    MPLB.flag_value(["--npart", "5", "--case"], "--case", "uniform") == "uniform" &&
        MPLB.flag_value(["--case", "uniform"], "--case", "shifted-normal") == "uniform",
    "flag_value takes the default where the flag is the last argument, and the value otherwise")
let options = (npart = 200, nbasis = 41, tstep = 8e-4, final_time = 3 * 8e-4, nu = 1.0),
    relaxed = MPLB.march(MPLB.CASES["shifted-normal"], DoughertyBracket(), options;
        record = "unused.hdf5")

    check(
        size(relaxed.z) == (200, 4) &&
            relaxed.final.particles.v[1, :] == relaxed.z[:, end] &&
            relaxed.initial.particles.v[1, :] == relaxed.z[:, 1] &&
            relaxed.z[:, end] != relaxed.z[:, 1],
        "march returns the initial particles at the first column of z and the final ones at " *
        "the last")
end

println("the options and the rate check of landau_cases.jl")
for (args, refused) in LANDAU_REFUSED
    edge_dir = mktempdir()
    runs, results = joinpath(edge_dir, "runs"), joinpath(edge_dir, "results")
    check(
        try
            LandauCases.landau_options([split(args); "--runs-dir"; runs; "--results-dir";
                                        results])
            false
        catch e
            e isa ArgumentError && startswith(e.msg, refused)
        end,
        "landau_options refuses $(args) with the ArgumentError $(refused)")
    check(!ispath(runs) && !ispath(results), "it creates neither --runs-dir nor --results-dir")
end
let rate = LandauCases.rate_check(1.5)
    check(rate.difference ≤ 1.5 && !rate.passes,
        "the rate check fails at the tolerance 1.5, above the relative difference to its " *
        "control: $(rate.control_difference)")
end

for job in jobs
    println(job.script)
    output = finish(job)
    expected = QUICK[job.script]
    foreach(name -> check_figure(joinpath(job.results, name)), expected.figures)
    for (name, parameters) in expected.records
        check_record(joinpath(job.runs, name), job.args, parameters)
    end
    if job.script == "pd_lb_fig11.jl"
        times, errors = cumulant_rows(output)
        check(times == [0.01, 0.02, 0.03],
            "pd_lb_fig11.jl prints the relative error of κ₄ at t = 0.01, 0.02, 0.03: $(times)")
        check(!isempty(errors) && all(<(5e-3), errors),
            "pd_lb_fig11.jl's relative errors lie below 5e-3: $(errors)")
        record = joinpath(job.runs, "pd_lb_fig11.hdf5")
        if isfile(record)
            residual = midpoint_residual(record, 1e-2)
            check(residual < MIDPOINT_RESIDUAL,
                "the first step of pd_lb_fig11.jl solves the implicit midpoint rule: $(residual)")
        end
    elseif job.script == "pd_lb_fig7_8.jl" &&
           isfile(joinpath(job.runs, "pd_lb_fig7_8.hdf5"))
        Random.seed!(1)
        drawn = initialize!(ParticleDistribution(1, 1, 200),
            DoubleMaxwellian(; domain = (-5.0, 5.0), shift = 2.0))
        z = h5read(joinpath(job.runs, "pd_lb_fig7_8.hdf5"), "z")
        check(z[:, 1] == drawn.particles.v[1, :],
            "pd_lb_fig7_8.jl starts from the velocities of the seed 1")
    elseif job.script == "pd_lb_fig1_2.jl"
        for record in ("pd_lb_fig1.hdf5", "pd_lb_fig2.hdf5")
            path = joinpath(job.runs, record)
            isfile(path) || continue
            check(!h5open(file -> haskey(HDF5.attrs(file), "stopped_at"), path, "r"),
                "$(record) records no stopped_at where the sweep does not stop")
        end
    end
end

println("pd_lb_fig11.jl in an empty working directory, outside a git checkout")
finish(stray)
check(isempty(readdir(WORK)), "the working directory stays empty: $(readdir(WORK))")
check(isfile(joinpath(stray.results, "pd_lb_fig11.pdf")),
    "--results-dir is created and holds the figure")
if check(isfile(joinpath(stray.runs, "pd_lb_fig11.hdf5")),
    "--runs-dir is created and holds the record")
    commit = h5open(file -> get(HDF5.attrs(file), "commit", nothing),
        joinpath(stray.runs, "pd_lb_fig11.hdf5"), "r")
    check(commit == "unknown", "commit is \"unknown\" outside a git checkout, got $(commit)")
    command = h5open(file -> get(HDF5.attrs(file), "command", nothing),
        joinpath(stray.runs, "pd_lb_fig11.hdf5"), "r")
    check(command == Base.shell_escape(stray.args...) && occursin("'", command),
        "command quotes the paths with a space: $(command)")
end

println("pd_lb_fig1_2.jl with sweeps that the projection stops")
output = finish(stop)
for (figure, record, dataset, reached, stops) in (
    ("pd_lb_fig1.pdf", "pd_lb_fig1.hdf5", "nbasis", [24, 26], [25]),
    ("pd_lb_fig2.pdf", "pd_lb_fig2.hdf5", "npart", [1000, 10000], [5623, 31623]))
    check(
        all(
            p -> occursin(Regex("^  $(dataset) = +$(p): the operator stops: ", "m"), output),
            stops),
        "it prints the stop at $(dataset) = $(stops)")
    check_figure(joinpath(stop.results, figure))
    check(isfile(joinpath(stop.runs, record)), "$(record) is written") || continue
    h5open(joinpath(stop.runs, record), "r") do file
        attributes = HDF5.attrs(file)
        check(haskey(attributes, "stopped_at") && attributes["stopped_at"] == first(stops),
            "$(record) records stopped_at = $(first(stops))")
        check(read(file[dataset]) == reached, "$(record) records $(dataset) = $(reached) alone")
    end
end

println("pd_lb_fig1_2.jl with a basis too small for the spline")
output = finish(other; fails = true, timed = false)
check(occursin("ArgumentError: at least two breakpoints are needed", output),
    "it stops with the ArgumentError of the basis")
check(!isfile(joinpath(other.results, "pd_lb_fig1.pdf")), "pd_lb_fig1.pdf is not written")

println("pd_lb_fig11.jl with a tolerance below its relative errors")
output = finish(miss; fails = true)
times, _ = cumulant_rows(output)
check(times == [0.01, 0.02, 0.03], "it prints the three relative errors: $(times)")

println("pd_lb_fig7_8.jl with a run that the projection stops")
output = finish(halt; fails = true, timed = false)
check(occursin("DomainError", output), "it stops with the DomainError of the projection")
foreach(name -> check(!isfile(joinpath(halt.results, name)), "$(name) is not written"),
    ["pd_lb_fig7.pdf", "pd_lb_fig8.pdf"])
check_record(joinpath(halt.runs, "pd_lb_fig7_8.hdf5"), halt.args,
    (npart = 60, tstep = 0.03, seed = 4, stopped_at = 5 * 0.03))
if isfile(joinpath(halt.runs, "pd_lb_fig7_8.hdf5"))
    Random.seed!(4)
    drawn = initialize!(ParticleDistribution(1, 1, 60),
        DoubleMaxwellian(; domain = (-5.0, 5.0), shift = 2.0))
    t = h5read(joinpath(halt.runs, "pd_lb_fig7_8.hdf5"), "t")
    z = h5read(joinpath(halt.runs, "pd_lb_fig7_8.hdf5"), "z")
    check(t == collect(0:4) .* 0.03 && size(z) == (60, 5),
        "it records the four steps reached: t = $(t), z of size $(size(z))")
    check(z[:, 1] == drawn.particles.v[1, :], "it records the velocities of the seed 4 first")
end

println("pd_lb_fig11.jl with a --runs-dir that cannot be created")
output = finish(unmade; fails = true, timed = false)
check(occursin("--runs-dir `$(unmade.runs)` cannot be created", output),
    "it names the path given: $(unmade.runs)")

for job in mplb
    case, bracket = job.args[3], job.args[5]
    println("$(job.script) --case $(case) --bracket $(bracket)")
    run_output = finish(job)
    run_stem = "lenard_bernstein_metriplectic_$(case)_$(bracket)"
    foreach(name -> check_figure(joinpath(job.results, "$(run_stem)_$(name).pdf")),
        MPLB_FIGURES)
    run_record = joinpath(job.runs, run_stem * ".hdf5")
    check_record(run_record, job.args, (; MPLB_QUICK[case]..., case, bracket))
    isfile(run_record) || continue
    check(!h5open(file -> haskey(HDF5.attrs(file), "stopped_at"), run_record, "r"),
        "$(run_stem).hdf5 records no stopped_at")
    if case == "double-maxwellian" && bracket == "dougherty"
        Random.seed!(1)
        drawn = initialize!(ParticleDistribution(1, 1, 1000),
            DoubleMaxwellian(; domain = (-5.0, 5.0), shift = 2.0))
        check(h5read(run_record, "z")[:, 1] == drawn.particles.v[1, :],
            "$(run_stem).hdf5 starts from the velocities of the seed 1")
    end
    case == "shifted-normal" || continue
    check_figure(joinpath(job.results, run_stem * "_rate.md"))
    rate = match(RATE_ROW, run_output)
    control = match(CONTROL_ROW, run_output)
    check(rate !== nothing && control !== nothing, "it prints the rate check and its control") ||
        continue
    fitted, closed = parse.(Float64, rate.captures)
    alternative_closed = parse(Float64, control[1])
    expected, alternative = bracket == "dougherty" ? (2.0, 1.0) : (1.0, 2.0)
    rate_temperature = h5read(run_record, "rate_temperature")
    check(isapprox(first(rate_temperature), 1.0; rtol = 1e-12),
        "the rate check starts at the temperature 1: $(first(rate_temperature))")
    check(
        isapprox(closed, expected; rtol = 1e-5) &&
            isapprox(alternative_closed, alternative; rtol = 1e-5),
        "the closed form is $(expected) and the control $(alternative): " *
        "it prints $(closed) and $(alternative_closed)")
    difference = abs(fitted - expected) / expected
    check(difference ≤ MPLB_RATE_TOLERANCE,
        "the fitted rate $(fitted) lies within 5 % of the closed form: $(difference)")
    control_difference = abs(fitted - alternative) / alternative
    check(control_difference > MPLB_RATE_TOLERANCE,
        "the fitted rate lies more than 5 % from the control: $(control_difference)")
end

for (job, flag, value) in ((unknown_case, "--case", "maxwellian"), (
    unknown_bracket, "--bracket", "landau"))
    println("$(job.script) with $(flag) $(value)")
    edge_output = finish(job; fails = true, timed = false)
    check(occursin("ArgumentError: unknown $(flag) `$(value)`", edge_output),
        "it stops with the ArgumentError that names $(flag) and $(value)")
    check(!ispath(job.runs) && !ispath(job.results),
        "it creates neither --runs-dir nor --results-dir")
end

println("$(MPLB_SCRIPT) --case uniform with --tolerance, which only the shifted normal takes")
let edge_output = finish(untolerant; fails = true, timed = false)
    check(occursin("ArgumentError: unknown argument `--tolerance`", edge_output),
        "it stops with the ArgumentError that names --tolerance")
    check(!ispath(untolerant.runs) && !ispath(untolerant.results),
        "it creates neither --runs-dir nor --results-dir")
end

println("$(MPLB_SCRIPT) with a Picard solve that does not converge")
let unsolved_output = finish(unsolved; fails = true, timed = false),
    stem = "lenard_bernstein_metriplectic_shifted-normal_dougherty",
    record = joinpath(unsolved.runs, stem * ".hdf5")

    check(
        occursin(
            "the Picard solve of the metriplectic Lenard–Bernstein step did not converge",
            unsolved_output),
        "it stops with the error of the Picard solve")
    foreach(
        name -> check(!isfile(joinpath(unsolved.results, "$(stem)_$(name).pdf")),
            "$(stem)_$(name).pdf is not written"),
        MPLB_FIGURES)
    check_record(record, unsolved.args, (npart = 1000, stopped_at = 6 * 8e-4))
    if isfile(record)
        t = h5read(record, "t")
        z = h5read(record, "z")
        check(t == collect(0:5) .* 8e-4 && size(z) == (1000, 6),
            "it records the five steps reached: t = $(t), z of size $(size(z))")
        check(
            h5open(
                file -> haskey(file, "rate_t") &&
                        length(read(file, "rate_t")) == 64 &&
                        length(read(file, "rate_temperature")) == 64,
                record,
                "r"),
            "it records the 64 samples of its rate check")
    end
end

for (job, condition) in ((tight, "closed form"), (loose, "control"))
    println("$(MPLB_SCRIPT) with $(join(job.args[1:3], " ")), where the $(condition) fails")
    edge_output = finish(job; fails = true)
    check(occursin("the rate check fails at the tolerance", edge_output),
        "it says that the rate check fails")
    stem = "lenard_bernstein_metriplectic_shifted-normal_dougherty"
    foreach(name -> check_figure(joinpath(job.results, "$(stem)_$(name).pdf")), MPLB_FIGURES)
    check(isfile(joinpath(job.runs, stem * ".hdf5")), "$(stem).hdf5 is written")
end

# The text that `RANKING` must not match: each driver's output and table, the docs section and
# the CHANGELOG entries.
ranked_texts = Pair{String, String}[]

for (job, (method, case)) in zip(landau, LANDAU_RUNS)
    println("$(LANDAU_SCRIPT) --method $(method) --case $(case)")
    run_output = finish(job)
    push!(ranked_texts, "the output of $(LANDAU_SCRIPT) --method $(method)" => run_output)
    stem = "landau_$(case)_$(method)_gamma-3"
    foreach(name -> check_figure(joinpath(job.results, "$(stem)_$(name).pdf")), LANDAU_FIGURES)
    record = joinpath(job.runs, stem * ".hdf5")
    check_record(record, job.args,
        (; LANDAU_QUICK..., method, case, gamma = -3.0, rate_check = 0))
    check(!occursin("rate check", run_output), "it prints no rate check at γ = −3")
    isfile(record) &&
        check(!h5open(file -> haskey(HDF5.attrs(file), "stopped_at"), record, "r"),
            "$(stem).hdf5 records no stopped_at")
end

println("$(LANDAU_SCRIPT) --method pairwise --gamma 0 --rate-check")
let rate_output = finish(landau_rate), stem = "landau_bimaxwellian_pairwise_gamma0",
    record = joinpath(landau_rate.runs, stem * ".hdf5"),
    table = joinpath(landau_rate.results, stem * "_rate.md")

    push!(ranked_texts, "the output of the rate check" => rate_output)
    foreach(name -> check_figure(joinpath(landau_rate.results, "$(stem)_$(name).pdf")),
        LANDAU_FIGURES)
    check_figure(table)
    isfile(table) && push!(ranked_texts, "the rate table" => read(table, String))
    check_record(record, landau_rate.args,
        (; LANDAU_QUICK..., method = "pairwise", gamma = 0.0, rate_check = 1))
    rate = match(LANDAU_RATE_ROW, rate_output)
    control = match(LANDAU_CONTROL_ROW, rate_output)
    if check(rate !== nothing && control !== nothing, "it prints the rate check and its control")
        λ, closed = parse.(Float64, rate.captures)
        alternative = parse(Float64, control[1])
        check(
            isapprox(closed, LANDAU_RATE; rtol = 1e-8) &&
                isapprox(alternative, LANDAU_CONTROL; rtol = 1e-8),
            "the closed form is 8 and the control 4: it prints $(closed) and $(alternative)")
        difference = abs(λ - LANDAU_RATE) / LANDAU_RATE
        check(difference ≤ LANDAU_RATE_TOLERANCE,
            "the instantaneous rate $(λ) lies within 2 % of 8: $(difference)")
        control_difference = abs(λ - LANDAU_CONTROL) / LANDAU_CONTROL
        check(control_difference > LANDAU_RATE_TOLERANCE,
            "the instantaneous rate lies more than 2 % from 4: $(control_difference)")
        isfile(record) && check(
            isapprox(
                h5open(file -> HDF5.attrs(file)["rate_measured"], record, "r"), λ;
                rtol = 1e-6),
            "$(stem).hdf5 records the measured rate")
    end
end

println("$(LANDAU_SCRIPT) with a step that does not converge")
let stop_output = finish(landau_stop; fails = true, timed = false),
    stem = "landau_bimaxwellian_pairwise_gamma-3",
    record = joinpath(landau_stop.runs, stem * ".hdf5")

    check(occursin("the discrete gradient iteration did not converge", stop_output),
        "it stops with the error of the discrete gradient iteration")
    foreach(
        name -> check(!isfile(joinpath(landau_stop.results, "$(stem)_$(name).pdf")),
            "$(stem)_$(name).pdf is not written"),
        LANDAU_FIGURES)
    check_record(record, landau_stop.args,
        (npart = 200, ncell = 8, tstep = 1e-3, stopped_at = 3 * 1e-3))
    if isfile(record)
        t = h5read(record, "t")
        v = h5read(record, "v")
        check(t == collect(0:2) .* 1e-3 && size(v) == (2, 200, 3),
            "it records the two steps reached: t = $(t), v of size $(size(v))")
    end
end

println("$(LANDAU_SCRIPT) with a rate check below its relative difference")
let tight_output = finish(landau_tight; fails = true),
    stem = "landau_bimaxwellian_pairwise_gamma0"

    check(occursin("the rate check fails at the tolerance", tight_output),
        "it says that the rate check fails")
    foreach(name -> check_figure(joinpath(landau_tight.results, "$(stem)_$(name).pdf")),
        LANDAU_FIGURES)
    check(isfile(joinpath(landau_tight.runs, stem * ".hdf5")), "$(stem).hdf5 is written")
end

println("$(COMPARE_SCRIPT)")
let compare_output = finish(compare), table = joinpath(compare.results, "landau_compare.md")
    push!(ranked_texts, "the output of $(COMPARE_SCRIPT)" => compare_output)
    foreach(name -> check_figure(joinpath(compare.results, name)),
        ["landau_compare_pairwise.pdf", "landau_compare_spline.pdf"])
    check_record(joinpath(compare.runs, "landau_compare.hdf5"), compare.args, COMPARE_QUICK)
    check_figure(table)
    if isfile(table)
        text = read(table, String)
        push!(ranked_texts, "the comparison's table" => text)
        rows = filter(!isnothing, map(line -> match(COMPARE_ROW, line), split(text, '\n')))
        sizes = Set((row[1], parse(Int, row[2]), parse(Int, row[3])) for row in rows)
        expected = Set((method, n^2, (ncell + 2)^2) for method in ("pairwise", "spline")
        for n in COMPARE_QUICK.lattice for ncell in COMPARE_QUICK.ncell)
        check(length(rows) == length(expected) && sizes == expected,
            "the table holds one row per method, lattice and cell count: $(length(rows)) rows")
        wrong = String[]
        for row in rows
            N, M, dN = parse.(Int, row.captures[2:4])
            measured, closed, err = parse.(Float64, row.captures[7:9])
            verdict = 2N ≥ M - 1 ? "holds" : "fails"
            dN == 2N && row[5] == verdict || push!(wrong, "$(row.match): rank")
            isapprox(closed, -LANDAU_RATE; atol = 1e-4) ||
                push!(wrong, "$(row.match): closed")
            isapprox(err, abs(measured + LANDAU_RATE) / LANDAU_RATE; rtol = 5e-3, atol = 2e-5) ||
                push!(wrong, "$(row.match): error")
            occursin(r"^[+-]", row[7]) || push!(wrong, "$(row.match): sign")
            row[6] == "yes" || push!(wrong, "$(row.match): support")
        end
        check(isempty(wrong),
            "every row has dN = 2N, the gate's verdict of dN ≥ M − 1, a particle in every " *
            "support, the closed form −8 and its relative error: $(wrong)")
        verdicts = Dict((parse(Int, row[2]), parse(Int, row[3])) => row[5] for row in rows)
        check(
            get(verdicts, (144, 289), "") == "holds" &&
                get(verdicts, (144, 324), "") == "fails",
            "N = 144 holds at M = 289, where dN = M − 1, and fails at M = 324")
        signs = Dict(row[1] => parse(Float64, row[7])
        for row in rows if (row[2], row[3]) == ("400", "324"))
        check(get(signs, "pairwise", 0.0) < 0 && get(signs, "spline", 0.0) > 0,
            "at N = 400, M = 324 the pairwise rate is negative and the spline rate positive: " *
            "$(signs)")
    end
end

println("the text that names the Landau methods")
let section = landau_docs_section(), entries = landau_changelog_entries()
    check(!isempty(section), "docs/src/collisions.md has the section $(LANDAU_DOCS_HEADING)")
    check(!isempty(entries), "CHANGELOG.md has an entry that names the Landau scripts")
    push!(ranked_texts, "the Landau section of docs/src/collisions.md" => section)
    foreach(entry -> push!(ranked_texts, "the CHANGELOG entry" => entry), entries)
    hits = [name * ": " * m.match for (name, text) in ranked_texts
            for m in eachmatch(RANKING, text)]
    check(isempty(hits), "no text ranks the two methods: $(hits)")
end

@printf("gate: %d failures, wall-clock time %.1f s\n", length(FAILURES), time() - STARTED)
exit(isempty(FAILURES) ? 0 : 1)
