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
# Each figure is one sweep. A point where the projection is not positive at a particle, the
# `DomainError` of `compute_coefficients`, is a point where the operator stops: the sweep prints
# it, leaves it out and goes on, the HDF5 file of the figure records the first such `M` or `N` as
# the attribute `stopped_at` beside the points it computed, and the figure draws those points. The
# script then exits 0. Any other error stops the script, which exits non-zero, and so does the
# `ArgumentError` of `plot_convergence` for a sweep that computes fewer than two points. With the
# defaults, figure 1 stops at `M = 25`, where the projected Maxwellian is negative at two of the
# particles, and figure 2 at `N = 5623` and `N = 31623`.
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

# The values of `f` at each of `points`: the points reached and their values. A point where `f`
# throws a `DomainError`, the model's own error, is skipped and its error printed; any other error
# passes through. The record `name` holds the points reached and their values as the datasets that
# `labels` names, and the first point that threw, where one did, as the attribute `stopped_at`.
function sweep(f, name, points, labels)
    reached = eltype(points)[]
    values = Float64[]
    stopped = (;)
    for p in points
        value = try
            f(p)
        catch e
            e isa DomainError || rethrow()
            @printf("  %s = %6d: the operator stops: %s\n", first(labels), p, e.msg)
            isempty(stopped) && (stopped = (; stopped_at = p))
            continue
        end
        push!(reached, p)
        push!(values, value)
        @printf("  %s = %6d: %s = %.4e\n", first(labels), p, last(labels), value)
    end
    write_record(name, merge(OPT, stopped); (Symbol.(labels) .=> (reached, values))...)
    return reached, values
end

@printf("pd_lb_fig1_2: %s run, %d particles, M = %s\n", QUICK ? "quick" : "paper",
    OPT.npart, OPT.nbasis)

dist = initial_distribution(NormalDistribution(), OPT.npart)
v = dist.particles.v[1, :]
M, norms = sweep("pd_lb_fig1.hdf5", OPT.nbasis, ("nbasis", "norm")) do M
    sdist = spline_distribution(M)
    project_Maxwellian(sdist)
    return norm(velocity_field(sdist, dist, v))
end
save_figure("pd_lb_fig1.pdf", plot_convergence(cell_width.(M), norms))

sdist = spline_distribution(OPT.sweep_nbasis)
N, variances = sweep("pd_lb_fig2.hdf5", OPT.sweep_npart, ("npart", "variance")) do N
    sample = initial_distribution(NormalDistribution(), N)
    u = sample.particles.v[1, :]
    projection(u, sample, sdist)
    return sum(abs2, velocity_field(sdist, sample, u)) / N
end
save_figure("pd_lb_fig2.pdf", plot_convergence(1 ./ N, variances))

report_time()
