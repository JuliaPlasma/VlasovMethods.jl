# The metriplectic Lenard–Bernstein driver: the particle distribution marched by the model's own
# Picard step, over a list of particle numbers, with the conservation error and the distribution
# at both ends of the run written out for each.
#
#   julia --project=scripts scripts/lenard_bernstein_metriplectic.jl \
#       [--npart 1000,3000,5000,7000,9000] [--nknot 31] [--order 4] [--tstep 0.0005] \
#       [--final-time 2.0] [--runs-dir dir] [--results-dir dir]
#
# For each particle number the HDF5 diagnostics go to `--runs-dir` (`<name>_N<n>.hdf5`) and the
# figure to `--results-dir` (`<name>_N<n>.png`), with a convergence figure and a report over the
# whole list at the end. The list is the scaling study: one run per particle number.
#
# **The initial condition is the cell-centred lattice on `[-2, 2]`.** A shifted double Maxwellian
# on `[-10, 10]` has tails that reach outside any bounded velocity domain; the projection then
# finds particles outside its support and throws, which is `KNOWN_ISSUES.md` K36. A lattice whose
# support *is* the velocity domain keeps the projected density positive everywhere the operator
# reads it, which is the device `test/integration/metriplectic_solve.jl` uses.
#
# The tolerance of the solve is the one the tests pin: `abstol = 3e-16·√N`, `reltol = 1e-50`.
# A fixed `1e-15` compares against `‖F‖₂` and is below the round-off floor of that norm: one step
# of the uniform cloud of `test/integration/metriplectic_solve.jl` throws at N = 64 with it.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using CairoMakie
using HDF5
using ProgressMeter
using VlasovMethods

const NAME = "lenard_bernstein_metriplectic"

const OPT = parse_arguments(ARGS; npart = collect(1000:2000:10000), nknot = 31, order = 4,
    tstep = 5e-4, final_time = 2.0)
const NPART_LIST = OPT.npart
const NKNOT = OPT.nknot
const ORDER = OPT.order
const TSTEP = OPT.tstep
const TSPAN = (0.0, OPT.final_time)
const DOMAINV = (-2.0, 2.0)

# the Picard iteration: the damping and depth the model's own step ignores, and the two
# tolerances the solve reads
const BETA = 0.5
const DEPTH = 3
const RELTOL = 1e-50

# The cell centres of the velocity domain, one per particle, with equal weights.
function lattice(npart)
    return [DOMAINV[1] + (i - 0.5) * (DOMAINV[2] - DOMAINV[1]) / npart for i in 1:npart]
end

"""
One run of `npart` particles over `TSPAN`; returns the times, the initial and the final state,
and the momentum and energy errors at every time.
"""
function march(npart)
    pdist = ParticleDistribution(1, 1, npart)
    pdist.particles.v[1, :] .= lattice(npart)
    pdist.particles.w[1, :] .= 1 / npart

    sdist = SplineDistribution(1, 1, NKNOT, ORDER, DOMAINV, 0, :nothing)
    model = MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist))

    trange = TSPAN[1]:TSTEP:TSPAN[2]

    # the unused `dv` and the extrapolation history, in the argument order of the solve
    dv = zeros(npart)
    dv_history = zeros(npart, 2)
    abstol = 3.0e-16 * sqrt(npart)     # the tolerance test/integration/metriplectic_solve.jl pins

    # Only the moments of each step are kept. The solve returns a new vector, so `v` and `v_old`
    # are the last two states without a copy.
    v₀ = pdist.particles.v[1, :]
    v, v_old = v₀, v₀
    mom, enr = [sum(v)], [sum(abs2, v)]
    @showprogress desc="N = $npart" for (i, t) in pairs(trange[1:(end - 1)])
        v, v_old = Picard_iterate_over_particles(
            dv, v, v_old, dv_history, i, t, TSTEP, DEPTH,
            BETA, abstol, RELTOL, model),
        v
        push!(mom, sum(v))
        push!(enr, sum(abs2, v))
    end

    # The lattice is symmetric about zero, so its momentum is zero to round-off and a *relative*
    # error in it is 0/0. The momentum error is therefore written as a difference, and the energy
    # error — the energy is not zero — relative.
    return collect(trange), v₀, v, (mom .- mom[1]) ./ npart, (enr .- enr[1]) ./ enr[1]
end

function main()
    started = time()
    final_momentum = Float64[]
    final_energy = Float64[]
    report = results_path(NAME * ".md")

    for npart in NPART_LIST
        trange, v_initial, v_final, momentum_error, energy_error = march(npart)

        h5file = runs_path("$(NAME)_N$(npart).hdf5")
        h5open(h5file, "w") do h5
            h5["t"] = trange
            h5["momentum_error"] = momentum_error
            h5["energy_error"] = energy_error
            h5["v"] = v_final
        end

        fig = Figure(size = (1200, 450))
        # an error of exactly zero has no logarithm; it is drawn at `eps()`, and a value below
        # `eps()` is drawn at that floor too. A floor at `floatmin()` would put the axis' lower
        # limit at the denormal bound, where Makie's log-limit padding underflows
        ax1 = Axis(fig[1, 1], xlabel = "t", ylabel = "error", title = "conservation",
            yscale = log10)
        lines!(ax1, trange, max.(abs.(momentum_error), eps()), color = :steelblue,
            label = "|momentum|")
        lines!(ax1, trange, max.(abs.(energy_error), eps()), color = :firebrick,
            label = "|energy|")
        axislegend(ax1; position = :rb)

        edges = collect(range(DOMAINV[1], DOMAINV[2]; length = 81))
        ax2 = Axis(fig[1, 2], xlabel = "v", ylabel = "f",
            title = "initial and final distribution")
        for (v, colour, label) in ((v_initial, :steelblue, "initial"),
            (v_final, :firebrick, "final"))
            hist!(ax2, v, bins = edges, normalization = :pdf, color = colour, label = label)
        end
        axislegend(ax2; position = :lt)

        figure = results_path("$(NAME)_N$(npart).png")
        save(figure, fig)

        push!(final_momentum, momentum_error[end])
        push!(final_energy, energy_error[end])
        @printf("%s: N = %6d, momentum error %.3e, energy error %.3e\n",
            NAME, npart, momentum_error[end], energy_error[end])
    end

    fig = Figure(size = (1000, 450))
    ax = Axis(fig[1, 1], xlabel = "N", ylabel = "error at t = $(TSPAN[2])",
        xscale = log10, yscale = log10)
    # the same `eps()` floor as the conservation figure, for the same reason
    scatterlines!(ax, NPART_LIST, max.(abs.(final_momentum), eps()), color = :steelblue,
        label = "|momentum|")
    scatterlines!(ax, NPART_LIST, max.(abs.(final_energy), eps()), color = :firebrick,
        label = "|energy|")
    axislegend(ax; position = :rt)

    convergence = results_path(NAME * "_convergence.png")
    save(convergence, fig)

    open(report, "w") do io
        println(io, "# The metriplectic Lenard–Bernstein run at ", length(NPART_LIST),
            " particle numbers")
        println(io)
        println(io, "`--nknot ", NKNOT, " --order ", ORDER, " --tstep ", TSTEP,
            " --final-time ", TSPAN[2], "`, the cell-centred lattice on ", DOMAINV,
            ", `abstol = 3e-16·√N`.")
        println(io)
        println(io, "| N | momentum error at t = ", TSPAN[2], " | energy error |")
        println(io, "|--:|--:|--:|")
        for (npart, momentum, energy) in zip(NPART_LIST, final_momentum, final_energy)
            @printf(io, "| %d | %+.3e | %+.3e |\n", npart, momentum, energy)
        end
        @printf(io, "\nWall time %.1f s.\n", time() - started)
    end

    @printf("%s: %d runs, %.1f s\n", NAME, length(NPART_LIST), time() - started)
    @printf("  wrote %s and %s\n", convergence, report)
end

main()
