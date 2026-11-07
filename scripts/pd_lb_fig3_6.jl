# Figures 3 to 6 of the Lenard–Bernstein manuscript (PD-LB): the relaxation of a shifted normal
# distribution under the conservative operator (`sec:first_example`).
#
#   julia --project=scripts scripts/pd_lb_fig3_6.jl [--quick] [--npart 1000] [--npart-fig6 200] \
#       [--nbasis 41] [--tstep 8e-4] [--final-time 1.0] [--nu 1.0] \
#       [--runs-dir dir] [--results-dir dir]
#
# The velocities are drawn from the normal distribution of mean 2 and variance 1
# (`ShiftedNormalV`), the spline has `--nbasis` cubic basis functions on [−10, 10], and
# `ConservativeLenardBernstein` is integrated with its implicit midpoint rule
# (`sec:time_discretisation`) at the step `--tstep` to `--final-time`.
#
# Figures 3, 4 and 5 come from the run of `--npart` particles: the initial and the final
# distribution (`fig:shifted_normal_example`), the relative drift of the mass, the momentum and
# the energy (`fig:conservation`), and the entropy `S_h` (`fig:entropy_ex1`). Figure 6 is the
# initial and the final distribution of the same run with `--npart-fig6` particles
# (`fig:shifted_normal_ex2`).
#
# Figure 5 departs from the manuscript's: the manuscript draws `S_h / |S_h(0)|`, and this script
# draws `S_h` itself, as `plot_entropy` does.
#
# `--quick` takes 200 particles and the final time 0.05. `ShiftedNormalV` seeds the random number
# generator itself, so the velocities are the same in every run.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

const QUICK, OPT = pd_lb_options(ARGS;
    paper = (npart = 1000, npart_fig6 = 200, nbasis = 41, tstep = 8e-4, final_time = 1.0,
        nu = 1.0),
    quick = (npart = 200, final_time = 0.05))

const EXAMPLE = ShiftedNormalV((-5.0, 5.0), 2.0)

@printf("pd_lb_fig3_6: %s run, %d and %d particles, M = %d, h = %g, t = %g\n",
    QUICK ? "quick" : "paper", OPT.npart, OPT.npart_fig6, OPT.nbasis, OPT.tstep,
    OPT.final_time)

relaxed = relax(EXAMPLE, OPT.npart, OPT; record = "pd_lb_fig3_5.hdf5")
d = relaxed.diagnostics
write_run("pd_lb_fig3_5.hdf5", OPT, relaxed)
save_distribution("pd_lb_fig3.pdf", relaxed, EXAMPLE.domain)
save_figure("pd_lb_fig4.pdf", plot_drift(relaxed.t, (d.mass, d.momentum, d.energy)))
save_figure("pd_lb_fig5.pdf", plot_entropy(relaxed.t, d.entropy))

const OPT6 = merge(OPT, (; npart = OPT.npart_fig6))
relaxed6 = relax(EXAMPLE, OPT6.npart, OPT6; record = "pd_lb_fig6.hdf5")
write_run("pd_lb_fig6.hdf5", OPT6, relaxed6)
save_distribution("pd_lb_fig6.pdf", relaxed6, EXAMPLE.domain)

report_time()
