# The allocation profile of one collision-model path.
#
#   julia --project=scripts scripts/profile.jl [--model clb|landau] [--npart 1000] [--nknot 41] \
#       [--order 4] [--tstep 0.001] [--final-time 0.001] [--runs-dir dir] [--results-dir dir]
#
# Run it in a fresh process; allocation profiling needs no `--track-allocation` or
# `--code-coverage` instrumentation, and either would measure the instrumentation. Read the
# report it writes to `--results-dir` against the raw dump it writes to `--runs-dir`. Everything
# the run produces goes under those two directories.
#
#   clb      a `ConservativeLenardBernstein` run through `GeometricIntegrator`, one step at the
#            defaults, so the profile is the whole collision step: the projection, the two
#            coefficient solves and the time step. The trajectory goes to `--runs-dir` beside the
#            dump.
#   landau   one `Picard_iterate_Landau_nls!` step of the two-dimensional Landau operator. No
#            Landau driver marches in time — every one of them meets the positivity guard of
#            `KNOWN_ISSUES.md` K8 — so the solver step is profiled directly, on the cell-centred
#            lattice of `test/integration/landau_solver.jl`, where the projection stays positive.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using Profile
using VlasovMethods

const OPT = parse_arguments(ARGS; model = "clb", npart = 1000, nknot = 41, order = 4,
    tstep = 1e-3, final_time = 1e-3)
const TSTEP = OPT.tstep
const TSPAN = (0.0, OPT.final_time)
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

One `ConservativeLenardBernstein` run over `tspan`, through the model's `GeometricIntegrator`,
with the trajectory written to `h5file`.
"""
function step_clb(h5file, tspan = TSPAN)
    dist = initialize!(ParticleDistribution(1, 1, OPT.npart), NormalDistribution())
    sdist = SplineDistribution(1, 1, OPT.nknot, OPT.order, DOMAINV, 0, :nothing)
    model = ConservativeLenardBernstein(dist, CollisionEntropy(sdist))
    integrator = GeometricIntegrator(model, tspan, TSTEP)

    return VlasovMethods.run!(integrator, h5file)
end

"One Landau Picard step of the two-dimensional setup."
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

# Each profile warms the path it profiles first, so that the profile is of arithmetic and not of
# compilation, and returns the line that describes it and the files it wrote.
function profile_clb()
    h5file = runs_path("profile_clb.hdf5")
    # the warm run writes the trajectory first, and the profiled run overwrites it
    step_clb(h5file, (TSPAN[1], TSPAN[1] + TSTEP))
    Profile.Allocs.clear()
    Profile.Allocs.@profile step_clb(h5file)

    return string(
        "`ConservativeLenardBernstein` through the model's `GeometricIntegrator`, ",
        "at --npart ", OPT.npart, ", --nknot ", OPT.nknot, ", --order ", OPT.order,
        ", --tstep ", TSTEP, ", --final-time ", TSPAN[2], ": ",
        round(Int, TSPAN[2] / TSTEP), " steps. The trajectory is `", basename(h5file),
        "`, beside the raw dump."),
    [h5file]
end

function profile_landau()
    landau, v_prev = landau_setup()
    # one warm call of each branch of the solver
    step_landau(landau, v_prev, 1)
    step_landau(landau, v_prev, 4)
    Profile.Allocs.clear()
    Profile.Allocs.@profile step_landau(landau, v_prev, 4)

    return string("One `Picard_iterate_Landau_nls!` step of the two-dimensional Landau ",
        "operator, on a ", N1, "^2 cell-centred lattice inside ", LANDAU_DOMAIN,
        " with a quadratic `Free` basis, `ti = 4`."),
    String[]
end

const PROFILES = Dict("clb" => profile_clb, "landau" => profile_landau)
const PROFILE = choose(PROFILES, "--model", OPT.model)

# The innermost frame in the repository's own source, `src/` or this script, so that a summary row
# names the code that allocated rather than the `Base` constructor or the runtime below it. An
# allocation with no such frame on its stack is counted under `outside VlasovMethods`.
const PACKAGE = pkgdir(VlasovMethods)

function allocation_frame(allocation)
    i = findfirst(frame -> startswith(string(frame.file), PACKAGE), allocation.stacktrace)
    i === nothing && return "outside VlasovMethods"
    frame = allocation.stacktrace[i]

    return string(frame.func, " at ", relpath(string(frame.file), PACKAGE), ":", frame.line)
end

"""
    allocation_summary(allocations; top = 15) -> (total, rows)

The `top` frames that own the most sampled bytes, each as a `(frame, count, bytes)` tuple, with the
sampled bytes over every frame.
"""
function allocation_summary(allocations; top = 15)
    frames = Dict{String, Tuple{Int, Int}}()
    for allocation in allocations
        frame = allocation_frame(allocation)
        count, bytes = get(frames, frame, (0, 0))
        frames[frame] = (count + 1, bytes + Int(allocation.size))
    end
    rows = sort!([(frame, count, bytes) for (frame, (count, bytes)) in frames];
        by = last, rev = true)

    return sum(allocation -> Int(allocation.size), allocations; init = 0), first(rows, top)
end

function main()
    started = time()
    description, written = PROFILE()
    elapsed = time() - started

    allocations = Profile.Allocs.fetch()
    dump = IOBuffer()
    # a wide display, so that the dump does not cut the file names of its frames
    Profile.Allocs.print(IOContext(dump, :displaysize => (10^6, 1000)), allocations)

    numbers = runs_path("profile_" * OPT.model * "_allocations.txt")
    report = results_path("profile_" * OPT.model * ".md")
    write(numbers, take!(dump))

    open(report, "w") do io
        println(io, "# The allocation profile of `scripts/profile.jl --model ", OPT.model, "`")
        println(io)
        println(io, description)
        @printf(io, "\nWall time %.1f s; %d sampled allocations.\n",
            elapsed, length(allocations.allocs))
        println(io, "The allocation profiler samples at its default rate of 0.1, so the count")
        println(io, "above is the sample and not the total.")
        println(io)
        total, rows = allocation_summary(allocations.allocs)
        @printf(io, "The frames below own %d of the %d sampled bytes:\n",
            sum(last, rows; init = 0), total)
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
        OPT.model, elapsed, length(allocations.allocs))
    @printf("  wrote %s\n", join([written; numbers; report], ", "))
end

main()
