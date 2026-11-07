# Figure 11 of the Lenard–Bernstein manuscript (PD-LB): the decay of the fourth cumulant under the
# operator with fixed coefficients, against the closed form of the cumulant appendix
# (`cumulant-scaling.tex`, `fig:analytic_cumulant_decay`).
#
#   julia --project=scripts scripts/pd_lb_fig11.jl [--quick] [--npart 1000] [--nbasis 41] \
#       [--tstep 5e-3] [--final-time 1.0] [--nu 1.0] [--seed 1] [--tolerance 0.05] \
#       [--runs-dir dir] [--results-dir dir]
#
# The appendix runs `eq:velocity_ode` with the coefficients fixed at `A₁ = 0` and `A₂ = 1`, in place
# of `eq:coefficient_lin_system`, "in our method": the implicit midpoint scheme of
# `sec:time_discretisation`. `LenardBernstein` is that right-hand side; its own `GeometricIntegrator`
# takes `RK438`, so this script integrates it with `ImplicitMidpoint` instead, as
# `implicit_midpoint` below builds it. The velocities are those of `pd_lb_fig7_8.jl`, half drawn
# from the normal distribution of mean 2 and half from that of mean −2, after seeding the random
# number generator with `--seed`, and the spline has `--nbasis` cubic basis functions on [−10, 10].
#
# The closed form: with the coefficients fixed, the cumulants `κ_k` of order `k > 2` decay as
# `κ_k(t) = κ_k(0) exp(−k ν t / σ²)`, where `σ² = 1 / A₂ = 1` is the variance of the equilibrium.
# The script compares the fourth cumulant of the particles, `κ₄ = μ₄ − 3 μ₂²` of their central
# moments, with `κ₄(0) exp(−4 ν t)` at the first three steps, prints the relative error of each,
# and exits non-zero where one exceeds `--tolerance`, 5 %, after it has written its record and its
# figure.
#
# Figure 11 is the error `|κ₄(t) − κ₄(0) exp(−4 ν t)| / |κ₄(0)|` and the closed form
# `exp(−4 ν t)`, both relative to `κ₄(0)`, against `t`, as `plot_error` draws them. It departs
# from the manuscript's figure, which draws `κ₄` and `κ₅` themselves, each divided by its initial
# value, against their closed forms as dashed lines. This script draws no `κ₅`: its initial value
# is sampling noise, because the bi-Maxwellian is symmetric.
#
# `--quick` takes 200 particles, the step 1e-2 and the final time 0.1.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

const QUICK, OPT = pd_lb_options(ARGS;
    paper = (npart = 1000, nbasis = 41, tstep = 5e-3, final_time = 1.0, nu = 1.0, seed = 1,
        tolerance = 0.05),
    quick = (npart = 200, tstep = 1e-2, final_time = 0.1))

const EXAMPLE = DoubleMaxwellian(; domain = (-5.0, 5.0), shift = 2.0)

# The implicit midpoint rule on the right-hand side of `LenardBernstein`, with the parameters of
# the model's own `GeometricIntegrator`.
function implicit_midpoint(model::LenardBernstein, tspan, tstep)
    parameters = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    equation = VlasovMethods.GeometricEquations.ODEProblem(VlasovMethods.LB_rhs_GI!, tspan,
        tstep, copy(model.dist.particles.v[1, :]); parameters)
    integrator = GeometricIntegrators.Integrators.GeometricIntegrator(equation,
        GeometricIntegrators.Integrators.ImplicitMidpoint())
    return GeometricIntegrator(model, equation, integrator)
end

# The fourth cumulant `μ₄ − 3 μ₂²` of the velocities `v` with the weights `w`.
function fourth_cumulant(v, w)
    n = sum(w)
    c = v .- sum(w .* v) / n
    μ₂ = sum(w .* c .^ 2) / n
    μ₄ = sum(w .* c .^ 4) / n
    return μ₄ - 3 * μ₂^2
end

@printf("pd_lb_fig11: %s run, %d particles, M = %d, h = %g, t = %g, seed %d\n",
    QUICK ? "quick" : "paper", OPT.npart, OPT.nbasis, OPT.tstep, OPT.final_time, OPT.seed)

relaxed = relax(EXAMPLE, OPT.npart, OPT; record = "pd_lb_fig11.hdf5", seed = OPT.seed,
    model = LenardBernstein, integrator = implicit_midpoint)
w = relaxed.initial.particles.w[1, :]
κ₄ = [fourth_cumulant(v, w) for v in eachcol(relaxed.z)]
closed = first(κ₄) .* exp.(-4 * OPT.nu .* relaxed.t)
relative = abs.(κ₄ .- closed) ./ abs.(closed)

write_record("pd_lb_fig11.hdf5", OPT; relaxed.t, relaxed.z, relaxed.diagnostics...,
    kappa4 = κ₄, kappa4_closed = closed)
save_figure("pd_lb_fig11.pdf",
    plot_error(relaxed.t, abs.(κ₄ .- closed) ./ abs(first(κ₄)), closed ./ first(κ₄)))

@printf("  κ₄(0) = %.6f\n", first(κ₄))
# the first three steps after `t = 0`, the columns 2 to 4
const COMPARED = 2:4
for k in COMPARED
    @printf("  κ₄ at t = %.4f: run %.6f, closed form %.6f, relative error %.3e\n",
        relaxed.t[k], κ₄[k], closed[k], relative[k])
end

report_time()

if maximum(relative[COMPARED]) > OPT.tolerance
    @printf("pd_lb_fig11: a relative error exceeds %g\n", OPT.tolerance)
    exit(1)
end
