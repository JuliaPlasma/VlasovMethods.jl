# Figures 9 and 10 of the Lenard–Bernstein manuscript (PD-LB): the relaxation of a uniform
# distribution under the conservative operator.
#
#   julia --project=scripts scripts/pd_lb_fig9_10.jl [--quick] [--npart 200] [--nbasis 41] \
#       [--tstep 1e-4] [--final-time 1.0] [--nu 1.0] [--seed 1] \
#       [--runs-dir dir] [--results-dir dir]
#
# The velocities are drawn from the uniform distribution on [−2, 2] (`UniformDistribution`), after
# seeding the random number generator with `--seed`. The spline has `--nbasis` cubic basis
# functions on [−10, 10], and `ConservativeLenardBernstein` is integrated with its implicit
# midpoint rule (`sec:time_discretisation`) at the step `--tstep` to `--final-time`.
#
# Figure 9 is the initial and the final distribution (`fig:initial_final_dists_uniform`), and
# figure 10 the entropy `S_h` (`fig:entropy_ex3`). Figure 10 departs from the manuscript's: the
# manuscript draws `S_h / |S_h(0)|`, and this script draws `S_h` itself, as `plot_entropy` does.
#
# `--quick` takes the step 1e-3 and the final time 0.05.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

const QUICK, OPT = pd_lb_options(ARGS;
    paper = (npart = 200, nbasis = 41, tstep = 1e-4, final_time = 1.0, nu = 1.0, seed = 1),
    quick = (tstep = 1e-3, final_time = 0.05))

const EXAMPLE = UniformDistribution((0.0, 1.0), (-2.0, 2.0))

@printf("pd_lb_fig9_10: %s run, %d particles, M = %d, h = %g, t = %g, seed %d\n",
    QUICK ? "quick" : "paper", OPT.npart, OPT.nbasis, OPT.tstep, OPT.final_time, OPT.seed)

relaxed = relax(EXAMPLE, OPT.npart, OPT; record = "pd_lb_fig9_10.hdf5", seed = OPT.seed)
write_run("pd_lb_fig9_10.hdf5", OPT, relaxed)
save_distribution("pd_lb_fig9.pdf", relaxed, EXAMPLE.xdomain)
save_figure("pd_lb_fig10.pdf", plot_entropy(relaxed.t, relaxed.diagnostics.entropy))

report_time()
