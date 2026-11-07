# The part that the drivers of the Lenard–Bernstein manuscript's figures share: `pd_lb_fig1_2.jl`,
# `pd_lb_fig3_6.jl`, `pd_lb_fig7_8.jl`, `pd_lb_fig9_10.jl` and `pd_lb_fig11.jl`. Each driver
# includes this file first. It holds the command line with `--quick`, the spline basis of
# `M` basis functions, the time stepping, the HDF5 record and the saving of a figure.
#
# Every HDF5 file a driver writes carries the attributes `commit`, the `git rev-parse HEAD` of
# the checkout that holds this file or `"unknown"` outside a git checkout, and `command`, the
# driver's own arguments joined by spaces, plus `script` and one attribute per option, so that the
# file names the command that regenerates it.

const STARTED = time()

include(joinpath(@__DIR__, "options.jl"))

using Printf
using Random

using CairoMakie
using HDF5
using VlasovMethods

import GeometricIntegrators

const VELOCITY_DOMAIN = (-10.0, 10.0)

# The order of the B-splines, 4 for the cubic splines of every figure.
const ORDER = 4

const COMMAND = join(ARGS, " ")
const SCRIPT = basename(PROGRAM_FILE)

function git_commit()
    try
        return readchomp(pipeline(`git -C $(@__DIR__) rev-parse HEAD`; stderr = devnull))
    catch
        return "unknown"
    end
end

const COMMIT = git_commit()

"""
    pd_lb_options(args; paper, quick) -> (quick::Bool, options::NamedTuple)

The options of a driver: `paper` holds the defaults of the manuscript, and `quick` the defaults
that the flag `--quick` replaces them with, for the short run of `scripts/gate.jl`. `--quick`
takes no value. Every other flag goes to `parse_arguments`, which also creates `--runs-dir` and
`--results-dir`.
"""
function pd_lb_options(args; paper::NamedTuple, quick::NamedTuple)
    n = count(==("--quick"), args)
    n ≤ 1 || throw(ArgumentError("`--quick` is given twice"))
    defaults = n == 1 ? merge(paper, quick) : paper
    return n == 1, parse_arguments(filter(!=("--quick"), args); defaults...)
end

"""
    spline_distribution(nbasis)

The cubic spline distribution on `VELOCITY_DOMAIN` with `nbasis` basis functions. A basis of
order `ORDER` on `ncell` cells has `ncell + ORDER - 1` functions, so the mesh has
`nbasis - ORDER + 2` breakpoints. The manuscript's "41 elements" are 41 basis functions.
"""
function spline_distribution(nbasis::Integer)
    sdist = SplineDistribution(1, 1, nbasis - ORDER + 2, ORDER, VELOCITY_DOMAIN)
    @assert VlasovMethods.nbasis(sdist.basis) == nbasis
    return sdist
end

"""
    initial_distribution(example, npart; seed = nothing)

`npart` particles of one position and one velocity dimension, drawn from `example`. With a
`seed`, the global random number generator is seeded with it first. `NormalDistribution` and
`ShiftedNormalV` seed it with `1234` themselves before they draw the velocities.
"""
function initial_distribution(example, npart::Integer; seed = nothing)
    seed === nothing || Random.seed!(seed)
    return initialize!(ParticleDistribution(1, 1, npart), example)
end

"""
    nsteps(final_time, tstep)

The number of steps of a run to `final_time`: `final_time / tstep` rounded up, after rounding the
quotient to eight digits, so that `0.07 / 5e-3`, which is `14.000000000000002` in floating point,
gives 14 steps and `0.05 / 8e-4` gives 63.
"""
nsteps(final_time, tstep) = ceil(Int, round(final_time / tstep; digits = 8))

"""
    integrate!(stopped, method::GeometricIntegrator, final_time, tstep) -> (t, z)

`nsteps(final_time, tstep)` steps of the integrator of `method` from its initial condition, and
the velocities of the particles at each step, written back into the model's particles at the end.
`t` holds the `n + 1` times `k * tstep`, and `z` the velocities, one column per time.

Where a step throws a `DomainError`, the model's own error, `stopped(t, z, time)` is called with
the times and the velocities of the steps reached and the time `k * tstep` of the step `k` that
threw, and the error is rethrown. The first step includes the extrapolation that starts the
solve, so a stop there gives `t = [0]`. Any other error passes through.

The steps are taken here and not by `run!`, because `run!` takes the steps of the range
`t₀:Δt:t_end`, which can hold one step fewer: 62 at `0.05` and `8e-4`, where `nsteps` gives 63.
"""
function integrate!(stopped, method::GeometricIntegrator, final_time, tstep)
    problem = method.equation
    n = nsteps(final_time, tstep)
    z = zeros(length(problem.ics.q), n + 1)
    z[:, 1] = problem.ics.q
    reached = 0
    try
        step = GeometricIntegrators.GeometricIntegratorsBase.solutionstep(
            method.integrator, VlasovMethods.GeometricEquations.initialstate(problem))
        for k in 1:n
            GeometricIntegrators.reset!(step, k * tstep)
            GeometricIntegrators.integrate!(step, method.integrator)
            z[:, k + 1] = step.q
            reached = k
        end
    catch e
        times = collect(0:reached) .* tstep
        e isa DomainError && stopped(times, z[:, 1:(reached + 1)], (reached + 1) * tstep)
        rethrow()
    end
    method.model.dist.particles.v[1, :] .= z[:, end]
    return collect(0:n) .* tstep, z
end

"""
    diagnostics(z, dist, sdist) -> NamedTuple

The mass `Σ w`, the momentum `Σ w v`, the energy `Σ w v² / 2` and the entropy
`S_h = ∫ f_s log|f_s| dv` of the projection onto `sdist`, at each column of `z`. `dist` gives the
weights, which no step changes; the projection overwrites the coefficients of `sdist`.
"""
function diagnostics(z, dist, sdist)
    w = dist.particles.w[1, :]
    entropy = CollisionEntropy(sdist)
    S = map(eachcol(z)) do v
        projection(collect(v), dist, sdist)
        entropy()
    end
    return (mass = fill(sum(w), size(z, 2)), momentum = vec(sum(w .* z; dims = 1)),
        energy = vec(sum(w .* z .^ 2; dims = 1)) ./ 2, entropy = S)
end

"""
    relax(example, npart, options; record, seed = nothing, model = ConservativeLenardBernstein,
        integrator = GeometricIntegrator) -> NamedTuple

One run of a collision operator from `npart` particles of `example`: the cubic spline of
`options.nbasis` functions, the collision frequency `options.nu`, and the integrator
`integrator(model, (0, options.final_time), options.tstep)`, by default the model's own, the
implicit midpoint rule for `ConservativeLenardBernstein`. Returns the initial and the final
particles, the times `t`, the velocities `z` and the `diagnostics`; the driver writes the record
`record` of a run that completes.

A projection that throws, on `f_s ≤ 0` at a particle or on a particle that left the velocity
domain, stops the run with its error. The record `record` is then written with the attribute
`stopped_at`, the time of the step that threw, and the datasets `t` and `z` of the steps reached,
before the error is rethrown, so that the driver writes no figure and exits non-zero.
"""
function relax(example, npart::Integer, options; record, seed = nothing,
        model = ConservativeLenardBernstein, integrator = GeometricIntegrator)
    dist = initial_distribution(example, npart; seed)
    initial = deepcopy(dist)
    sdist = spline_distribution(options.nbasis)
    collisions = model(dist, CollisionEntropy(sdist); ν = options.nu)
    method = integrator(collisions, (0.0, options.final_time), options.tstep)
    t, z = integrate!(method, options.final_time, options.tstep) do t, z, time
        @printf("  the operator stops at t = %g, after %d steps\n", time, length(t) - 1)
        write_record(record, merge(options, (; stopped_at = time)); t, z)
    end
    return (initial, final = dist, t, z, diagnostics = diagnostics(z, dist, sdist))
end

"""
    write_record(name, options; datasets...)

The HDF5 file `name` in `--runs-dir`, with the attributes `commit`, `command` and `script`, one
attribute per entry of `options`, and one dataset per keyword.
"""
function write_record(name, options; datasets...)
    path = runs_path(name)
    h5open(path, "w") do file
        attributes = HDF5.attrs(file)
        attributes["commit"] = COMMIT
        attributes["command"] = COMMAND
        attributes["script"] = SCRIPT
        for (key, value) in pairs(options)
            attributes[String(key)] = value
        end
        for (key, value) in pairs(datasets)
            file[String(key)] = value
        end
    end
    @printf("  wrote %s\n", path)
    return path
end

"""
    write_run(name, options, relaxed)

The record of a run of `relax`: its times `t`, its velocities `z` and each of its diagnostics as
a dataset.
"""
function write_run(name, options, relaxed)
    return write_record(name, options; relaxed.t, relaxed.z, relaxed.diagnostics...)
end

"""
    save_figure(name, fig)

`fig` as the file `name` in `--results-dir`.
"""
function save_figure(name, fig)
    path = results_path(name)
    save(path, fig)
    @printf("  wrote %s\n", path)
    return path
end

"""
    save_distribution(name, relaxed, xdomain)

The initial and the final particles of a run of `relax` as the figure of
`plot_collision_distribution` on `xdomain` and `VELOCITY_DOMAIN`, saved as `name`.
"""
function save_distribution(name, relaxed, xdomain)
    return save_figure(name,
        plot_collision_distribution(relaxed.initial, relaxed.final; xdomain,
            vdomain = VELOCITY_DOMAIN))
end

"""
    report_time()

Prints the wall-clock time since this file was included, the package loads included.
"""
report_time() = @printf("%s: wall-clock time %.1f s\n", SCRIPT, time() - STARTED)
