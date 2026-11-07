# Figures 1 and 2 of the Lenard–Bernstein manuscript (PD-LB): the convergence of the
# semi-discrete conservative operator, under projection and under particle sampling
# (`fig:spline-convergence`, `fig:particle-convergence`).
#
#   julia --project=scripts scripts/pd_lb_fig1_2.jl [--quick] [--npart 100000] \
#       [--nbasis 24,25,…,78] [--sweep-nbasis 41] [--sweep-npart 100,178,…,100000] \
#       [--runs-dir dir] [--results-dir dir]
#
# Figure 1: `npart` velocities drawn from the standard normal distribution, which are never
# projected. For each number of cubic basis functions `M` of `--nbasis` on [−10, 10], the exact
# Maxwellian is projected onto the basis (`eq:numerical-example-1-coefficients`), and the right-hand
# side `eq:velocity_maxwellian_ode` is evaluated at the particles, with `A₁` and `A₂` of
# `eq:discrete_A_solution`. The figure is its Euclidean norm `‖v̇‖₂` against the cell width
# `20 / (M − 3)`, with the fitted exponent; the manuscript's reference curve is `M⁻⁴`.
#
# Figure 2: for each `N` of `--sweep-npart`, `N` normal velocities are projected onto
# `--sweep-nbasis` basis functions (`eq:projection_f`) and `eq:velocity_ode` is evaluated at them.
# The figure is the sample variance `Σ v̇² / N` against `1 / N`, with the fitted exponent; the
# manuscript's reference curve is `1 / N`.
#
# Each figure is one sweep, and an error stops it: the HDF5 file of the figure then records the
# `M` or the `N` that stopped as the attribute `stopped_at`, beside the values before it, the figure
# is not written, and the script exits non-zero with the error after both sweeps have run. A
# projection that is not positive at a particle is such an error, the `DomainError` of
# `compute_coefficients`. With the defaults, both sweeps stop on it: figure 1 at `M = 25`, where
# the projected Maxwellian is negative at two of the particles, and figure 2 at `N = 5623`.
#
# `--quick` takes 2000 particles, `M ∈ {24, 41, 78}` and the first six `N` of the sweep, 100 to
# 1778, where neither sweep stops.
#
# `NormalDistribution` seeds the random number generator itself, so each `N` has the same
# velocities in every run.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

using LinearAlgebra: norm

const QUICK, OPT = pd_lb_options(ARGS;
    paper = (npart = 100_000, nbasis = collect(24:78), sweep_nbasis = 41,
        sweep_npart = round.(Int, 10 .^ (2:0.25:5))),
    quick = (npart = 2000, nbasis = [24, 41, 78],
        sweep_npart = round.(Int, 10 .^ (2:0.25:3.25))))

# The right-hand side of `eq:velocity_ode` at the velocities `v` of `dist`, for the spline of
# `sdist` as it stands, with `ν = 1`.
function velocity_field(sdist, dist, v)
    A₁, A₂ = VlasovMethods.compute_coefficients(sdist, dist, v)
    fs = sdist.spline
    dfs = derivative(fs)
    return -(dfs.(v) ./ fs.(v) .+ A₁ .+ A₂ .* v)
end

# cell width of a cubic basis of `M` functions on the velocity domain
cell_width(M) = (last(VELOCITY_DOMAIN) - first(VELOCITY_DOMAIN)) / (M - ORDER + 1)

# The values of `f` at each of `points` until `f` throws: the points reached, their values, and the
# error, or `nothing`. The record `name` holds the points and the values as the datasets that
# `labels` names, and the point that threw as the attribute `stopped_at`.
function sweep(f, name, points, labels)
    values = Float64[]
    failure = nothing
    for p in points
        try
            push!(values, f(p))
        catch e
            failure = e
            break
        end
        @printf("  %s = %6d: %s = %.4e\n", first(labels), p, last(labels), last(values))
    end
    reached = points[1:length(values)]
    stopped = failure === nothing ? (;) : (; stopped_at = points[length(values) + 1])
    isempty(stopped) ||
        @printf("  %s = %6d: the sweep stops\n", first(labels), stopped.stopped_at)
    write_record(name, merge(OPT, stopped); (Symbol.(labels) .=> (reached, values))...)
    return reached, values, failure
end

@printf("pd_lb_fig1_2: %s run, %d particles, M = %s\n", QUICK ? "quick" : "paper",
    OPT.npart, OPT.nbasis)

dist = initial_distribution(NormalDistribution(), OPT.npart)
v = dist.particles.v[1, :]
M, norms, error1 = sweep("pd_lb_fig1.hdf5", OPT.nbasis, ("nbasis", "norm")) do M
    sdist = spline_distribution(M)
    project_Maxwellian(sdist)
    return norm(velocity_field(sdist, dist, v))
end
error1 === nothing &&
    save_figure("pd_lb_fig1.pdf", plot_convergence(cell_width.(M), norms))

sdist = spline_distribution(OPT.sweep_nbasis)
N, variances, error2 = sweep("pd_lb_fig2.hdf5", OPT.sweep_npart, ("npart", "variance")) do N
    sample = initial_distribution(NormalDistribution(), N)
    u = sample.particles.v[1, :]
    projection(u, sample, sdist)
    return sum(abs2, velocity_field(sdist, sample, u)) / N
end
error2 === nothing && save_figure("pd_lb_fig2.pdf", plot_convergence(1 ./ N, variances))

report_time()

for e in (error1, error2)
    e === nothing || throw(e)
end
