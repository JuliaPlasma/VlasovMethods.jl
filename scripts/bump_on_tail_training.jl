# The training set of the paper's bump-on-tail example: the bump-on-tail of `κ = 0.3` marched full
# order at the ten time scales `χ = κₛ/0.3` of `κₛ ∈ [0.1, 0.5]`, and written to the runs
# directory. `bump_on_tail_projections.jl` builds the reduced bases from it, and
# `bump_on_tail_testing.jl` marches the reduced model of those bases.
#
#   julia --project=scripts scripts/bump_on_tail_training.jl [--np 10000] [--nh 16] [--order 3] \
#       [--dt 0.1] [--final-time 25.0] [--seed 1234] [--runs-dir dir] [--results-dir dir]
#
# The run id of the file comes from the parameters, so two runs of different parameters never
# share one. It draws no figure.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

const OPT = parse_arguments(ARGS; np = BOT_NP, nh = BOT_NH, order = BOT_ORDER, dt = BOT_DT,
    final_time = BOT_T, seed = BOT_SEED)

const PARAMS = (
    np = OPT.np, nh = OPT.nh, order = OPT.order, dt = OPT.dt, T = OPT.final_time,
    κlo = KAPPA_RANGE[1], κhi = KAPPA_RANGE[2], seed = OPT.seed)

sweep = training_sweep(PARAMS)
fpath = runs_path("$(run_id(PARAMS)).h5")

h5open(fpath, "w") do io
    h5save(io, sweep.snapshots)
end

@printf("bump_on_tail_training: %d samples, %d particles, %d cells, degree %d\n",
    length(sweep.kappa), PARAMS.np, PARAMS.nh, PARAMS.order)
@printf("  dt = %g, t ∈ [0, %g]\n", PARAMS.dt, PARAMS.T)
@printf("  χ of the samples: %s\n", join(string.(round.(sweep.chi; digits = 4)), ", "))
@printf("  wrote %s\n", fpath)
