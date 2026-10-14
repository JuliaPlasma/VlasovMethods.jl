# The allocation profile of one collision-model path.
#
#   julia --project=scripts scripts/profile.jl [--model clb|landau] [--npart 1000] [--nknot 41] \
#       [--order 4] [--tstep 0.001] [--final-time 0.01] [--runs-dir dir] [--results-dir dir]
#
# Run it in a fresh process; allocation profiling needs no `--track-allocation` or
# `--code-coverage` instrumentation, and either would measure the instrumentation. Read the
# report it writes to `--results-dir` against the raw dump it writes to `--runs-dir`. Everything
# the run produces goes under those two directories.
#
#   clb      one `ConservativeLenardBernstein` run through `GeometricIntegrator`, so the profile
#            is the whole collision step: the projection, the two coefficient solves and the
#            time step. The trajectory goes to `--runs-dir` beside the dump.
#   landau   one `Picard_iterate_Landau_nls!` step of the two-dimensional Landau operator. No
#            Landau driver marches in time — every one of them meets the positivity guard of
#            `KNOWN_ISSUES.md` K8 — so the solver step is profiled directly, on the cell-centred
#            lattice of `test/integration/landau_solver.jl`, where the projection stays positive.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using Profile
using VlasovMethods

const MODELS = ("clb", "landau")

const PARSED = let (positional, options) = parse_arguments(ARGS,
        ("--model", "--npart", "--nknot", "--order", "--tstep", "--final-time"))
    no_positional(positional, "profile.jl")
    options
end

const MODEL = option_choice(PARSED, "--model", "clb", MODELS)
const NPART = option_int(PARSED, "--npart", 1000)
const NKNOT = option_int(PARSED, "--nknot", 41)
const ORDER = option_int(PARSED, "--order", 4)
const TSTEP = option_float(PARSED, "--tstep", 1e-3)
const TSPAN = (0.0, option_float(PARSED, "--final-time", 1e-2))
const DOMAINV = (-10.0, 10.0)

# The particle lattice of `test/integration/landau_solver.jl`'s two-dimensional Landau step: 16
# cell centres per velocity axis inside `[-2, 2]`, equal weights. The basis is quadratic over 16
# cells where the test's is cubic over 8, and `Free` on the same interval, so that the projected
# density is positive everywhere the operator reads it.
const N1 = 16
const LANDAU_DOMAIN = (-2.0, 2.0)
const TOL = 1e-10
const FTOL = 1e-10
const BETA = 0.5
const DEPTH = 3
const CHUNKSIZE = 1

"""
    step_clb(h5file, tspan = TSPAN)

One `ConservativeLenardBernstein` run over `tspan`, through the model's `GeometricIntegrator`.
Returns what `run!(integrator, h5file)` returns; the trajectory is at `h5file`. The caller wraps in
`Profile.Allocs.@profile` the call whose allocations it wants.
"""
function step_clb(h5file, tspan = TSPAN)
    dist = initialize!(ParticleDistribution(1, 1, NPART), NormalDistribution())
    sdist = SplineDistribution(1, 1, NKNOT, ORDER, DOMAINV, 0, :Dirichlet)
    entropy = CollisionEntropy(sdist)
    model = ConservativeLenardBernstein(dist, entropy)
    integrator = GeometricIntegrator(model, tspan, TSTEP)

    return VlasovMethods.run!(integrator, h5file)
end

"One Landau Picard step of the two-dimensional setup, warmed by the call in `main`."
function step_landau(landau, v_prev, ti)
    rhs_prev = zeros(2, size(v_prev, 2), 2)
    # `Picard_iterate_Landau_nls!` prints and reports its own convergence; K43 records that.
    redirect_stdout(devnull) do
        VlasovMethods.Picard_iterate_Landau_nls!(
            landau, TOL, FTOL, BETA, TSTEP, ti, 0.0, v_prev, copy(v_prev), rhs_prev,
            DEPTH, CHUNKSIZE)
    end

    return nothing
end

function landau_setup()
    npart = N1^2
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pdist = ParticleDistribution(1, 2, VlasovMethods.ParticleList(zeros(4, npart);
        variables = vars))
    lattice = [LANDAU_DOMAIN[1] + (i - 0.5) * (LANDAU_DOMAIN[2] - LANDAU_DOMAIN[1]) / N1
               for i in 1:N1]
    pdist.particles.v[1, :] .= repeat(lattice; outer = N1)
    pdist.particles.v[2, :] .= repeat(lattice; inner = N1)
    pdist.particles.w[1, :] .= 1 / npart

    # `:nothing` is `Free()`, the condition `test/integration/landau_solver.jl` builds by hand.
    sdist = SplineDistribution(1, 2, N1 + 1, 3, LANDAU_DOMAIN, 0, :nothing)

    return Landau(pdist, CollisionEntropy(sdist)), copy(pdist.particles.v)
end

# The innermost Julia frame of an allocation — the runtime and the allocation profiler below and
# above it are C, so a summary row names the code under profile and not the sampling machinery.
function allocation_frame(allocation)
    for frame in allocation.stacktrace
        frame.from_c && continue
        return string(frame.func, " at ", basename(string(frame.file)), ":", frame.line)
    end

    return "unknown"
end

"""
    allocation_summary(allocations; top = 15) -> (total, rows)

The `top` frames that own the most sampled bytes, each as a `(frame, count, bytes)` tuple, with the
sampled bytes over every frame.
"""
function allocation_summary(allocations; top = 15)
    frames = Dict{String, Tuple{Int, Int}}()
    total = 0
    for allocation in allocations
        frame = allocation_frame(allocation)
        count, bytes = get(frames, frame, (0, 0))
        frames[frame] = (count + 1, bytes + Int(allocation.size))
        total += Int(allocation.size)
    end
    rows = sort([(frame, count, bytes) for (frame, (count, bytes)) in frames];
        by = row -> -row[3])

    return total, rows[1:min(top, length(rows))]
end

function main()
    started = time()
    h5file = MODEL == "clb" ? runs_path("profile_clb.hdf5") : ""

    if MODEL == "clb"
        # a warm one-step run of the same path, so that the profile is of arithmetic and not of
        # compilation; it writes the trajectory first, and the profiled run overwrites it
        step_clb(h5file, (TSPAN[1], TSPAN[1] + TSTEP))
        Profile.Allocs.clear()
        Profile.Allocs.@profile step_clb(h5file)
    else
        landau, v_prev = landau_setup()
        # a warm call of each branch of the solver, so that the profile is of arithmetic and not
        # of compilation
        step_landau(landau, v_prev, 1)
        step_landau(landau, v_prev, 4)
        Profile.Allocs.clear()
        Profile.Allocs.@profile step_landau(landau, v_prev, 4)
    end

    elapsed = time() - started

    allocations = Profile.Allocs.fetch()
    dump = IOBuffer()
    Profile.Allocs.print(dump, allocations)
    profile_text = String(take!(dump))

    numbers = runs_path("profile_" * MODEL * "_allocations.txt")
    report = results_path("profile_" * MODEL * ".md")
    write(numbers, profile_text)

    open(report, "w") do io
        println(io, "# The allocation profile of `scripts/profile.jl --model ", MODEL, "`")
        println(io)
        if MODEL == "clb"
            println(io, "`ConservativeLenardBernstein` through the model's ",
                "`GeometricIntegrator`, at --npart ", NPART, ", --nknot ", NKNOT,
                ", --order ", ORDER, ", --tstep ", TSTEP, ", --final-time ", TSPAN[2],
                ": ", round(Int, TSPAN[2] / TSTEP), " steps.")
            println(io, "The trajectory is `", basename(h5file), "`, beside the raw dump.")
        else
            println(
                io, "One `Picard_iterate_Landau_nls!` step of the two-dimensional Landau ",
                "operator, on a ", N1, "^2 cell-centred lattice inside ", LANDAU_DOMAIN,
                " with a quadratic `Free` basis, `ti = 4`.")
        end
        @printf(io, "\nWall time %.1f s; %d sampled allocations.\n",
            elapsed, length(allocations.allocs))
        println(io, "The allocation profiler samples at its default rate of 0.1, so the count")
        println(io, "above is the sample and not the total.")
        println(io)
        total, rows = allocation_summary(allocations.allocs)
        @printf(io, "The frames below own %d of the %d sampled bytes:\n",
            sum(row -> row[3], rows; init = 0), total)
        println(io)
        println(io, "| allocations | bytes | frame |")
        println(io, "|--:|--:|:--|")
        for (frame, count, bytes) in rows
            @printf(io, "| %d | %d | `%s` |\n", count, bytes, frame)
        end
        println(io)
        println(io, "The full dump is `", basename(numbers), "` in the runs directory.")
    end

    @printf("profile: model %s, %.1f s, %d sampled allocations\n",
        MODEL, elapsed, length(allocations.allocs))
    @printf("  wrote %s%s and %s\n", h5file == "" ? "" : h5file * ", ", numbers, report)
end

main()
