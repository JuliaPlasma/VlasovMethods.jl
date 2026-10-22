# The reduced bases of the paper's bump-on-tail example, and P39's two input files: the particle
# basis of the training set at the notebook's tolerance, the field basis at its own, and the
# reduced positions and accelerations of every training sample. It also draws the particle and the
# field spectrum (item 11).
#
# The paper's command, whose basis tolerances are `VP-ROM-Projections.ipynb`'s cells 14 and 15:
#
#   julia --project=scripts scripts/bump_on_tail_projections.jl --np 10000 --nh 16 --order 3 \
#       --dt 0.1 --final-time 25.0 --particle-tol 1e-9 --field-tol 1e-5
#
# `--particle-tol 1e-9` is the tolerance of the particle basis and `--field-tol 1e-5` that of the
# field basis; the old `scripts/bump_on_tail_2_projections.jl` of the reference used `1e-4` and
# `1e-2`. The script marches the ten training samples itself, as `bump_on_tail_training.jl` does,
# so it needs no earlier run; P39 reads the two files it writes.
#
#   julia --project=scripts scripts/bump_on_tail_projections.jl [--np 10000] [--nh 16] \
#       [--order 3] [--dt 0.1] [--final-time 25.0] [--seed 1234] [--particle-tol 1e-9] \
#       [--field-tol 1e-5] [--runs-dir dir] [--results-dir dir]

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie

const OPT = parse_arguments(ARGS; np = BOT_NP, nh = BOT_NH, order = BOT_ORDER, dt = BOT_DT,
    final_time = BOT_T, seed = BOT_SEED, particle_tol = 1e-9, field_tol = 1e-5)

const PARAMS = (
    np = OPT.np, nh = OPT.nh, order = OPT.order, dt = OPT.dt, T = OPT.final_time,
    κlo = KAPPA_RANGE[1], κhi = KAPPA_RANGE[2], seed = OPT.seed)

sweep = training_sweep(PARAMS)
SS = sweep.snapshots

# the two bases of the paper: the cotangent lift of the positions and velocities at the particle
# tolerance, and the POD basis of the accelerations at the field tolerance
particle = particle_basis(CotangentLiftEVD(), SS;
    parameters = (κ = BOT_KAPPA, α = 0.1, σ = 0.5, v₀ = 4.5), tolerance = OPT.particle_tol)
field = field_basis(PODEVD(), SS;
    parameters = (κ = BOT_KAPPA, α = 0.1, σ = 0.5, v₀ = 4.5), tolerance = OPT.field_tol)

# the reduced data of P39: every saved state projected onto the particle basis
k, nₛ, nparam = particle.k, size(SS.X, 3), length(sweep.kappa)
X = zeros(k, nₛ, nparam)
A = zeros(k, nₛ, nparam)
for p in axes(X, 3), ts in axes(X, 2)

    mul!(@view(X[:, ts, p]), particle.Φ', @view(SS.X[1, :, ts, p]))
    mul!(@view(A[:, ts, p]), particle.Φ', @view(SS.A[1, :, ts, p]))
end

ppath, vpath = write_p39_files(RUNS_DIR[], PARAMS, particle, field, X, A;
    commit = git_commit(normpath(joinpath(@__DIR__, ".."))), command = join(ARGS, " "),
    seed = OPT.seed)

@printf("bump_on_tail_projections: %d particles, %d cells, dt = %g, t ∈ [0, %g]\n",
    PARAMS.np, PARAMS.nh, PARAMS.dt, PARAMS.T)
@printf("  particle basis: k = %d at tolerance %g\n", particle.k, OPT.particle_tol)
@printf("  field basis:    k = %d at tolerance %g\n", field.k, OPT.field_tol)
@printf("  wrote %s\n", ppath)
@printf("  wrote %s\n", vpath)

# item 11: the particle and the field spectrum, each against its own largest eigenvalue
_normalised(v) = max.(v ./ first(v), eps(Float64))
_pad(v, n) = vcat(v, fill(NaN, n - length(v)))
n = max(length(particle.Λ), length(field.Λ))
figure = plot_series(1:n,
    [log10.(_pad(_normalised(particle.Λ), n)), log10.(_pad(_normalised(field.Λ), n))];
    labels = ["particles", "field"], xlabel = "i", ylabel = "log₁₀ λᵢ/λ₁")
fpath = results_path(item_figure(11, "spectra"))
save(fpath, figure)
@printf("  wrote %s\n", fpath)
