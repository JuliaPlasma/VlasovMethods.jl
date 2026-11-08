# The full-order Vlasov–Poisson driver: a particle cloud on a periodic spline basis, advanced by
# the splitting method, with the phase-space scatter of the cloud drawn at four times of the run.
#
#   julia --project=scripts scripts/vlasov_poisson.jl [--npart 10000] [--ncells 16] [--order 3] \
#       [--tstep 0.1] [--final-time 20.0] [--seed 1234] [--runs-dir dir] [--results-dir dir]
#
# The HDF5 trajectory goes to `--runs-dir`, the figure to `--results-dir`. `--seed` seeds the
# draw of the particle positions, so one command writes the same trajectory and figure.

include(joinpath(@__DIR__, "options.jl"))

using Printf
using Random

using CairoMakie
using HDF5
using PoissonSolvers
using VlasovMethods

"The options of the driver from the command line `args`, defaulting to the values above."
function driver_options(args)
    return parse_arguments(args; npart = 10000, ncells = 16, order = 3, tstep = 0.1,
        final_time = 20.0, seed = 1234)
end

const DOMAIN = (0.0, 1.0)

function main(o)
    tspan = (0.0, o.final_time)

    # create and initialize particle distribution function; `NormalDistribution` draws the
    # positions before its own `Random.seed!`, so the generator is seeded here
    Random.seed!(o.seed)
    dist = initialize!(ParticleDistribution(1, 1, o.npart), NormalDistribution())

    # create electrostatic potential
    potential = Potential(PeriodicBasisSpline(DOMAIN, o.order, o.ncells))

    # create Vlasov-Poisson model
    model = VlasovPoisson(dist, potential)

    # create integrator
    integrator = SplittingMethod(model, tspan, o.tstep)

    # integrate, writing the trajectory to the runs directory
    h5file = runs_path("vlasov_poisson.hdf5")
    run!(integrator, h5file)

    @printf("vlasov_poisson: %d particles, %d cells, order %d, dt = %g, t ∈ %s\n",
        o.npart, o.ncells, o.order, o.tstep, tspan)

    # `z` is `(nd, npart, ntime)`; the two rows are `x` and `v`. `run!(::SplittingMethod, h5file)`
    # writes the trajectory alone, so the time of a slice is the step count times the step size.
    z = h5read(h5file, "z")
    @printf("  wrote %s (%d time slices)\n", h5file, size(z, 3))
    slice_time(slice) = tspan[1] + (slice - 1) * o.tstep

    # the figure shows the same cloud at four evenly spaced times, so that the filamentation the
    # scheme is meant to resolve is visible rather than asserted
    slices = unique(round.(Int, range(1, size(z, 3); length = 4)))
    fig = plot_phase_space(z, slices; xdomain = DOMAIN,
        titles = [@sprintf("t = %.1f", slice_time(slice)) for slice in slices])

    figure = results_path("vlasov_poisson.png")
    save(figure, fig)
    @printf("  wrote %s\n", figure)
    return nothing
end

# `scripts/gate.jl` includes this file to read `driver_options`, so `main` runs only where the file
# is the program
abspath(PROGRAM_FILE) == @__FILE__() && main(driver_options(ARGS))
