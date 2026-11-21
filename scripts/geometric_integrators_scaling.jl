# The cost of a GeometricIntegrators step on the state of the grid 1D1V Vlasov–Poisson model, at
# three grid sizes, against a hand-written RK4 loop with the same tableau. Run from the repository
# root:
#
#   julia --project=scripts scripts/geometric_integrators_scaling.jl [--results-dir dir] [--runs-dir dir]
#
# Each size runs the Landau problem of `scripts/landau_damping_1d1v.jl` in `Float64`, with the step
# `dt = 0.5 hₓ / 8.6`. The integrator's RK4 and the loop take the same steps from the same state,
# alternating in samples of `@timed`, each sample through a function barrier; the time per step is
# the fastest sample's and the bytes per step the most that any sample allocated. The two states
# must agree to `8 · nsteps · eps` relative to the loop's largest entry, or the run stops with no
# record. The implicit midpoint rule is measured on its own where its dense Newton Jacobian,
# `8 (nₓ n_v)²` bytes, fits in a quarter of the memory.
#
# The verdict is "accepted" when the integrator's time per step is at most `RATIO_LIMIT` times the
# loop's at every size and its bytes per step do not increase from one size to the next, and
# "rejected" otherwise. The record, with the machine, the load and the versions of the packages a
# step runs, goes to `--results-dir` as `geometric_integrators_scaling.toml`; the script exits 0
# for either verdict. Run it in a fresh process, as the measurement is of compiled code.

include(joinpath(@__DIR__, "options.jl"))

using GeometricIntegrators
using LinearAlgebra
using PoissonSolvers
using TOML
using VlasovMethods

const GI = GeometricIntegrators
const GIBase = GeometricIntegrators.GeometricIntegratorsBase
const GE = GeometricIntegrators.GeometricEquations

"The options of the driver from the command line `args`: the two output directories alone."
driver_options(args) = parse_arguments(args)

"The grid sizes `(nₓ, n_v)`, in increasing `nₓ n_v`."
const SIZES = [(64, 128), (256, 512), (1024, 1024)]

"The largest ratio of the integrator's time per step to the loop's that the verdict accepts."
const RATIO_LIMIT = 1.2

"The file name of the record in `--results-dir`."
const RECORD = "geometric_integrators_scaling.toml"

# The packages whose code a step runs, by the names the record gives their versions under.
const PACKAGES = ["VlasovMethods", "GeometricIntegrators", "GeometricIntegratorsBase",
    "GeometricBase", "GeometricEquations", "GeometricBrackets", "PoissonSolvers",
    "SimpleSolvers", "RungeKutta"]

# The Landau problem of `scripts/landau_damping_1d1v.jl`.
const K = 0.5
const ALPHA = 0.01
const L = 4π
const VDOMAIN = (-8.6, 8.6)

"""
    scaling_verdict(ratios, bytes; limit = RATIO_LIMIT) -> String

`"accepted"` if every ratio is at most `limit` and the bytes per step never increase from one size
to the next, `bytes[k + 1] ≤ bytes[k]`, with both inputs in the order of increasing `nₓ n_v`;
`"rejected"` otherwise, a `NaN` ratio included. Empty inputs and inputs of different lengths throw
an `ArgumentError`.
"""
function scaling_verdict(ratios, bytes; limit = RATIO_LIMIT)
    isempty(ratios) && throw(ArgumentError("scaling_verdict needs at least one size"))
    length(ratios) == length(bytes) || throw(ArgumentError(
        "scaling_verdict takes one ratio and one byte count per size; got $(length(ratios)) " *
        "ratios and $(length(bytes)) byte counts"))
    ratios_hold = all(r -> r <= limit, ratios)
    bytes_hold = all(bytes[k + 1] <= bytes[k] for k in 1:(length(bytes) - 1))
    return ratios_hold && bytes_hold ? "accepted" : "rejected"
end

"The verdict of the RK4 entries of a record: their `ratio` and `integrator_bytes_per_step`."
function rk4_verdict(entries)
    return scaling_verdict([e["ratio"] for e in entries],
        [e["integrator_bytes_per_step"] for e in entries])
end

"The bytes of the dense `Float64` Jacobian of the state of an `nx × nv` grid."
jacobian_bytes(nx, nv) = 8 * (nx * nv)^2

"""
    implicit_fits(nx, nv, memory) -> Bool

Whether the implicit midpoint rule is measured at `nx × nv`: its Jacobian, which the Newton solve
holds about three times over with the LU factor and the dual-number buffers, takes at most a
quarter of `memory` bytes.
"""
implicit_fits(nx, nv, memory) = jacobian_bytes(nx, nv) <= memory ÷ 4

"""
    per_step(times, bytes, nsteps) -> (; seconds, bytes)

The time per step of the fastest sample and the bytes per step of the sample that allocated the
most, for samples of `nsteps` steps each.
"""
function per_step(times, bytes, nsteps)
    (seconds = minimum(times) / nsteps,
        bytes = maximum(bytes) / nsteps)
end

"""
    check_agreement(q, qref, nsteps) -> Float64

The difference `d = maximum(abs, q - qref) / maximum(abs, qref)` of two states advanced by the
same `nsteps` steps from the same state. Throws an `ErrorException` when `d` exceeds
`8 · nsteps · eps(Float64)` or is `NaN`: a timing of two codes that compute different states
compares nothing.
"""
function check_agreement(q, qref, nsteps)
    d = maximum(abs, q - qref) / maximum(abs, qref)
    tol = 8 * nsteps * eps(Float64)
    d <= tol ||
        error("the integrator's and the loop's RK4 states differ by $(d) relative " *
              "after $(nsteps) steps, above the bound $(tol)")
    return d
end

"The output of the program `command`, stripped; `\"unavailable\"` for `nothing`."
uptime_text(::Nothing) = "unavailable"
uptime_text(command) = String(strip(read(`$command`, String)))

"""
    write_record(dir, record) -> String

Write the dictionary `record` as TOML to `joinpath(dir, RECORD)`, replacing an existing file, and
return the path.
"""
function write_record(dir, record)
    path = joinpath(dir, RECORD)
    open(path, "w") do io
        TOML.print(io, record)
    end
    return path
end

# The grid model of the Landau problem at `nx × nv`, and its step `0.5 hₓ / 8.6`.
function landau_model(nx, nv)
    basis = FFTWBasis((0.0, L), nx)
    x = basis.xgrid[1:nx]
    v = range(VDOMAIN[1], VDOMAIN[2]; length = nv)
    values = zeros(nx, nv)
    values .= (1 .+ ALPHA .* cos.(K .* x)) .* exp.(-v' .^ 2 ./ 2) ./ sqrt(2π)
    dist = GridDistribution(values, (0.0, L), VDOMAIN)
    model = VlasovPoisson(dist, Potential(basis))
    return model, 0.5 * step(dist.x) / VDOMAIN[2]
end

# The solution step of `int`, at its initial state and time.
function solution_step(int)
    prob = int.equation
    return GIBase.solutionstep(int.integrator, GE.initialstate(prob)), GE.initialtime(prob)
end

# `nsteps` steps of the integrator, the first one ending at `t₀ + (n + 1) dt`: the clock grows by
# `dt` at every step, because the implicit midpoint rule's initial guess extrapolates from the
# times of the steps before and throws on two equal times.
function integrator_steps!(solstep, integrator, t₀, dt, n, nsteps)
    for s in (n + 1):(n + nsteps)
        GI.reset!(solstep, t₀ + s * dt)
        GI.integrate!(solstep, integrator)
    end
    return nothing
end

# `nsteps` steps of the classic RK4 tableau on `q`, in place, with the buffers `work` and the
# right-hand side of the grid model on `params`. The grid model's right-hand side ignores the time.
function rk4_steps!(q, work, params, dt, nsteps)
    (; k₁, k₂, k₃, k₄, y) = work
    rhs! = VlasovMethods._grid_1d1v_rhs!
    for _ in 1:nsteps
        rhs!(k₁, nothing, q, params)
        @. y = q + dt / 2 * k₁
        rhs!(k₂, nothing, y, params)
        @. y = q + dt / 2 * k₂
        rhs!(k₃, nothing, y, params)
        @. y = q + dt * k₃
        rhs!(k₄, nothing, y, params)
        @. q += dt / 6 * (k₁ + 2 * k₂ + 2 * k₃ + k₄)
    end
    return nothing
end

# The RK4 measurement of one size: `warmup` steps of both codes, then `samples` samples of `steps`
# steps each, alternating the integrator and the loop. Returns the per-step statistics of both and
# the agreement of their final states.
function measure_rk4(
        solstep, integrator, q, work, params, t₀, dt, (; warmup, samples, steps))
    integrator_steps!(solstep, integrator, t₀, dt, 0, warmup)
    rk4_steps!(q, work, params, dt, warmup)
    times = (zeros(samples), zeros(samples))
    bytes = (zeros(samples), zeros(samples))
    for k in 1:samples
        a = @timed integrator_steps!(solstep, integrator, t₀, dt, warmup + (k - 1) * steps,
            steps)
        b = @timed rk4_steps!(q, work, params, dt, steps)
        times[1][k], bytes[1][k] = a.time, a.bytes
        times[2][k], bytes[2][k] = b.time, b.bytes
    end
    nsteps = warmup + samples * steps
    d = check_agreement(solstep.q, q, nsteps)
    return per_step(times[1], bytes[1], steps), per_step(times[2], bytes[2], steps), d
end

# The measurement of an integrator alone: `warmup` steps, then `samples` samples of `steps` steps.
function measure_integrator(solstep, integrator, t₀, dt, (; warmup, samples, steps))
    integrator_steps!(solstep, integrator, t₀, dt, 0, warmup)
    times = zeros(samples)
    bytes = zeros(samples)
    for k in 1:samples
        a = @timed integrator_steps!(solstep, integrator, t₀, dt, warmup + (k - 1) * steps,
            steps)
        times[k], bytes[k] = a.time, a.bytes
    end
    return per_step(times, bytes, steps)
end

# The version of each package of `PACKAGES`, read from the loaded modules.
function package_versions()
    modules = Dict(id.name => m for (id, m) in Base.loaded_modules)
    return Dict{String, Any}(name => string(pkgversion(modules[name])) for name in PACKAGES)
end

"""
    scaling_record(sizes; memory = Sys.total_memory(), uptime_command = Sys.which("uptime"),
                   rk4 = (warmup = 2, samples = 5, steps = 10),
                   midpoint = (warmup = 1, samples = 3, steps = 1)) -> Dict

The measurement at each `(nₓ, n_v)` of `sizes`, with the machine, the versions and the verdict, as
a dictionary that `write_record` writes. `rk4` and `midpoint` give the warm-up steps, the samples
and the steps per sample of each method; the implicit midpoint rule is measured where
`implicit_fits(nₓ, n_v, memory)`. Throws an `ArgumentError`, before it measures anything, for no
sizes and for sizes whose `nₓ n_v` does not strictly increase.
"""
function scaling_record(
        sizes; memory = Sys.total_memory(), uptime_command = Sys.which("uptime"),
        rk4 = (warmup = 2, samples = 5, steps = 10),
        midpoint = (warmup = 1, samples = 3, steps = 1))
    isempty(sizes) && throw(ArgumentError("scaling_record needs at least one size"))
    counts = [nx * nv for (nx, nv) in sizes]
    all(counts[k] < counts[k + 1] for k in 1:(length(counts) - 1)) || throw(ArgumentError(
        "scaling_record takes sizes in strictly increasing nₓ n_v; got $(collect(sizes))"))

    rk4_entries = Dict{String, Any}[]
    midpoint_entries = Dict{String, Any}[]
    for (nx, nv) in sizes
        model, dt = landau_model(nx, nv)

        nsteps = rk4.warmup + rk4.samples * rk4.steps
        int = VlasovMethods.GeometricIntegrator(model, (0.0, (nsteps + 1) * dt), dt, RK4())
        solstep, t₀ = solution_step(int)
        q = copy(vec(model.distribution.values))
        work = (k₁ = similar(q), k₂ = similar(q), k₃ = similar(q), k₄ = similar(q),
            y = similar(q))
        integrator, loop, d = measure_rk4(solstep, int.integrator, q, work,
            int.equation.parameters, t₀, dt, rk4)
        push!(rk4_entries,
            Dict{String, Any}(
                "nx" => nx, "nv" => nv, "dt" => dt, "warmup_steps" => rk4.warmup,
                "samples" => rk4.samples, "steps_per_sample" => rk4.steps,
                "integrator_seconds_per_step" => integrator.seconds,
                "loop_seconds_per_step" => loop.seconds,
                "ratio" => integrator.seconds / loop.seconds,
                "integrator_bytes_per_step" => integrator.bytes,
                "loop_bytes_per_step" => loop.bytes, "state_difference" => d))

        entry = Dict{String, Any}("nx" => nx, "nv" => nv,
            "jacobian_bytes" => jacobian_bytes(nx, nv),
            "measured" => implicit_fits(nx, nv, memory))
        if entry["measured"]
            nsteps = midpoint.warmup + midpoint.samples * midpoint.steps
            int = VlasovMethods.GeometricIntegrator(model, (0.0, (nsteps + 1) * dt), dt,
                ImplicitMidpoint())
            solstep, t₀ = solution_step(int)
            stats = measure_integrator(solstep, int.integrator, t₀, dt, midpoint)
            entry["integrator_seconds_per_step"] = stats.seconds
            entry["integrator_bytes_per_step"] = stats.bytes
        end
        push!(midpoint_entries, entry)
    end

    return Dict{String, Any}(
        "verdict" => rk4_verdict(rk4_entries), "ratio_limit" => RATIO_LIMIT,
        "sizes" => [[nx, nv] for (nx, nv) in sizes],
        "julia_version" => string(VERSION), "threads" => Threads.nthreads(),
        "blas_threads" => LinearAlgebra.BLAS.get_num_threads(),
        "uptime" => uptime_text(uptime_command), "load_average" => collect(Sys.loadavg()),
        "total_memory" => Int(memory), "packages" => package_versions(),
        "rk4" => rk4_entries, "implicit_midpoint" => midpoint_entries)
end

function main(options)
    record = scaling_record(SIZES)
    path = write_record(options.results_dir, record)
    for e in record["rk4"]
        println("RK4 ", e["nx"], " × ", e["nv"], ": integrator ",
            e["integrator_seconds_per_step"], " s, loop ", e["loop_seconds_per_step"],
            " s per step, ratio ", e["ratio"], "; bytes per step ",
            e["integrator_bytes_per_step"], " and ", e["loop_bytes_per_step"])
    end
    println("verdict: ", record["verdict"], " (ratio limit ", RATIO_LIMIT, ")")
    println("  wrote ", path)
    return nothing
end

# The test of this script includes it, so `main` runs only where the file is the program.
if abspath(PROGRAM_FILE) == @__FILE__
    main(driver_options(ARGS))
end
