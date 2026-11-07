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
# Four more runs check the edges of the drivers:
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
#   are above `4e-4`: it prints the three relative errors and exits non-zero.
#
# The quick record of `pd_lb_fig11.jl` must hold the tolerance `0.05`.
#
# Beside the runs, the gate includes `pd_lb_common.jl` and checks `--quick`, the step count and the
# spline basis directly.
#
# The drivers run at once, so the gate takes about as long as the slowest of them. It prints the
# size of each file and its own wall-clock time, and exits 1 when a check fails. A driver still
# running when the gate exits is killed.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

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

# A driver still running when the gate exits, after an error of the gate itself, is killed.
atexit(() -> foreach(
    job -> process_running(job.process) && kill(job.process, Base.SIGKILL),
    [jobs; [stray, stop, other, miss]]))

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

@printf("gate: %d failures, wall-clock time %.1f s\n", length(FAILURES), time() - STARTED)
exit(isempty(FAILURES) ? 0 : 1)
