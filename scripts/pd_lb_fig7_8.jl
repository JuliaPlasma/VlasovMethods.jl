# Figures 7 and 8 of the Lenard–Bernstein manuscript (PD-LB): the relaxation of a bi-Maxwellian
# under the conservative operator.
#
#   julia --project=scripts scripts/pd_lb_fig7_8.jl [--quick] [--npart 1000] [--nbasis 41] \
#       [--tstep 8e-4] [--final-time 5.0] [--nu 1.0] [--seed 1] \
#       [--runs-dir dir] [--results-dir dir]
#
# Half of the velocities are drawn from the normal distribution of mean 2 and variance 1, and half
# from that of mean −2 (`DoubleMaxwellian`), after seeding the random number generator with
# `--seed`. The spline has `--nbasis` cubic basis functions on [−10, 10], and
# `ConservativeLenardBernstein` is integrated with its implicit midpoint rule
# (`sec:time_discretisation`) at the step `--tstep` to `--final-time`.
#
# Figure 7 is the initial and the final distribution (`fig:double_maxwellian_example`), and
# figure 8 the entropy `S_h` (`fig:entropy_ex2`). The page `docs/src/collisions.md` states how
# each figure departs from the manuscript's.
#
# Where the projection throws, the run stops with its `DomainError`: `relax` writes the record
# `pd_lb_fig7_8.hdf5` with the attribute `stopped_at`, the time of the step that threw, and the
# steps reached, and the script writes no figure and exits non-zero. With the manuscript's
# parameters the run stops so at `t = 1.4144` (K94).
#
# `--quick` takes 200 particles and the final time 0.2.

include(joinpath(@__DIR__, "pd_lb_common.jl"))

const QUICK, OPT = pd_lb_options(ARGS;
    paper = (npart = 1000, nbasis = 41, tstep = 8e-4, final_time = 5.0, nu = 1.0, seed = 1),
    quick = (npart = 200, final_time = 0.2))

const EXAMPLE = DoubleMaxwellian(; domain = (-5.0, 5.0), shift = 2.0)

@printf("pd_lb_fig7_8: %s run, %d particles, M = %d, h = %g, t = %g, seed %d\n",
    QUICK ? "quick" : "paper", OPT.npart, OPT.nbasis, OPT.tstep, OPT.final_time, OPT.seed)

relaxed = relax(EXAMPLE, OPT.npart, OPT; record = "pd_lb_fig7_8.hdf5", seed = OPT.seed)
write_run("pd_lb_fig7_8.hdf5", OPT, relaxed)
save_distribution("pd_lb_fig7.pdf", relaxed, EXAMPLE.domain)
save_figure("pd_lb_fig8.pdf", plot_entropy(relaxed.t, relaxed.diagnostics.entropy))

report_time()
