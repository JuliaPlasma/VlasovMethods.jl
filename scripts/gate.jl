# The short run of the figure drivers of the collision manuscripts: each driver with `--quick`, in
# a process of its own, into a temporary directory.
#
#   julia --startup-file=no --project=scripts scripts/gate.jl
#
# A driver passes when it exits 0 within `BUDGET` seconds and writes each of its figures and each
# of its HDF5 files, and each HDF5 file carries
#
# - the attribute `commit`, equal to `git rev-parse HEAD` of this checkout,
# - the attribute `command`, equal to the arguments the gate passed to the driver,
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
#   levels below a directory that does not have them, and with `GIT_DIR` naming a path that does
#   not exist: both directories are created and hold its files, the working directory stays
#   empty, and `commit` is `"unknown"`;
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
# `DoubleMaxwellian` draws after `Random.seed!(1)`. Five more runs check its edges:
#
# - `--case maxwellian` and `--bracket landau`: each exits non-zero with the `ArgumentError` that
#   names the flag and the value, and creates neither `--runs-dir` nor `--results-dir`;
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
# `ArgumentError` where the temperature at the first sample is zero or negative, that
# `flag_value` takes the default where the flag is the last argument, and that `march` of the 200
# particles of the quick shifted normal over three steps returns the initial particles at the
# first column of its velocities and the final ones at the last.
#
# Beside the runs, the gate includes `pd_lb_common.jl` and checks `--quick`, the step count and the
# spline basis directly.
#
# The drivers run at once, so the gate takes about as long as the slowest of them. It prints the
# size of each file and its own wall-clock time, and exits 1 when a check fails. A driver still
# running when the gate exits is killed.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

module MPLB
include(joinpath(@__DIR__, "lenard_bernstein_metriplectic.jl"))
end

const ROOT = dirname(@__DIR__)
const BUDGET = 240.0
const MIDPOINT_RESIDUAL = 1e-10
const EXPECTED_COMMIT = readchomp(`git -C $(ROOT) rev-parse HEAD`)

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
        check(haskey(attributes, "command") && attributes["command"] == join(args, " "),
            "$(basename(path)) carries command = $(join(args, " "))")
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
    parameters = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model)
    v̇ = similar(z[:, 1])
    VlasovMethods.LB_rhs_GI!(v̇, 0.0, (z[:, 1] .+ z[:, 2]) ./ 2, parameters)
    return maximum(abs, z[:, 2] .- z[:, 1] .- tstep .* v̇) / maximum(abs, z[:, 2] .- z[:, 1])
end

jobs = [start(script, ["--quick"]) for script in sort!(collect(keys(QUICK)))]

const WORK = mktempdir()
const NESTED = joinpath(mktempdir(), "new", "deeper")
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
unsolved = start(MPLB_SCRIPT, ["--quick", "--bracket", "dougherty", "--npart", "1000"])
tight = start(MPLB_SCRIPT, ["--quick", "--tolerance", "1e-3"])
loose = start(MPLB_SCRIPT, ["--quick", "--tolerance", "1.5"])

# A driver still running when the gate exits, after an error of the gate itself, is killed.
atexit(() -> foreach(
    job -> process_running(job.process) && kill(job.process, Base.SIGKILL),
    [jobs; mplb;
     [stray, stop, other, miss, halt, unmade, unknown_case, unknown_bracket,
         unsolved, tight, loose]]))

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
for first_temperature in (0.0, -1.0)
    check(
        try
            MPLB.fit_rate([0.0, 1.0], [first_temperature, 0.5], 1.0)
            false
        catch e
            e isa ArgumentError
        end,
        "fit_rate throws an ArgumentError where the first temperature is $(first_temperature)")
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

@printf("gate: %d failures, wall-clock time %.1f s\n", length(FAILURES), time() - STARTED)
exit(isempty(FAILURES) ? 0 : 1)
