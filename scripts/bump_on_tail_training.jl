# The training sweep of the self-consistent bump-on-tail set: the ten χ of the paper's sweep, each
# a full-order march of the same initial condition, saved into one HDF5 file. The figure draws the
# field energy of each sample, and the fitted growth rate of each is printed.
#
#   julia --project=scripts scripts/bump_on_tail_training.jl [--particles 10000] [--ncells 16] \
#       [--degree 3] [--tstep 0.1] [--final-time 25.0] [--runs-dir dir] [--results-dir dir]
#
# The default is the paper's sweep. The snapshots go to `--runs-dir` as `<runid>.h5`, which
# `bump_on_tail_projections.jl` reads from the same directory, and the figure to `--results-dir`.

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie
using HDF5
using PoissonSolvers: PeriodicBasisSpline, Potential
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods

const OPT = parse_arguments(ARGS; particles = NPARTICLES_BOT, ncells = NCELLS_BOT,
    degree = DEGREE_BOT, tstep = DT_BOT, final_time = T_BOT)

const T = Float64
const NPARTICLES = OPT.particles
const DT = OPT.tstep
const NT = round(Int, OPT.final_time / DT)
const NS = NT + 1
const CHI = training_chi()

basis = PeriodicBasisSpline((zero(T), domain_length()), OPT.degree, OPT.ncells)
potential = Potential(basis)
field = vp_field(potential)
nh = nbasis(basis)

dist = initial_distribution(; n = NPARTICLES)
ps = ParameterSpace(
    Parameter(:χ, first(CHI), last(CHI), CHI),
    Parameter(:ε, ε_BOT, ε_BOT, 1),
    Parameter(:a, α_BOT, α_BOT, 1),
    Parameter(:v₀, v₀_BOT, v₀_BOT, 1),
    Parameter(:σ, σ_BOT, σ_BOT, 1))

SS = Snapshots(T, 1, NPARTICLES, nh, NS, ps)
IP = IntegratorParameters(DT, NT, NS, nh, NPARTICLES, length(ps))

@printf("bump_on_tail_training: %d particles, %d cells, degree %d, dt = %g, t ∈ [0, %g]\n",
    NPARTICLES, OPT.ncells, OPT.degree, DT, OPT.final_time)

for p in eachindex(ps)
    integrate_vp(dist.particles, march_parameters(ps(p).χ), field, SS, IP, p)
end

# The acceleration of each saved state, which the field basis of the projections is built from:
# the plain field at the saved positions, with the weights the march deposited. No march writes
# `Snapshots.A`, so it is taken here from the positions the march saved.
a = zeros(T, NPARTICLES)
w = vec(dist.particles.w)
for p in eachindex(ps), ts in axes(SS.X, 3)

    x = @view SS.X[1, :, ts, p]
    field(a, x, w, IP.t[ts])
    SS.A[1, :, ts, p] .= a
end

# the snapshots, with the commit, the command and the seed the files of the set carry
fpath = runs_path("$(RUNID_BOT).h5")
h5open(fpath, "w") do file
    VlasovMethods.h5save(file, SS)
    VlasovMethods.h5save(file, IP)
    HDF5.attributes(file)["commit"] = commit_id()
    HDF5.attributes(file)["command"] = join(ARGS, " ")
    HDF5.attributes(file)["seed"] = SEED_BOT
end
@printf("  wrote %s\n", fpath)

# the growth rate of every sample, half the slope of log W at the local maxima of the field energy.
# The fit takes the first two maxima: the slowest sample of the sweep, at the smallest χ, reaches
# only two within T = 25, and one range serves every sample.
γ = growth_rate(IP.t, SS.W; maxima = 1:2)
for p in eachindex(ps)
    @printf("  χ = %.4f  W(T) = %.4e  γ = %.4f\n", ps(p).χ, SS.W[end, p], γ[p])
end

fig = plot_series(IP.t, eachcol(SS.W);
    labels = [@sprintf("χ = %.3f", χ) for χ in CHI], xlabel = "t", ylabel = "W")
figure = results_path("bump_on_tail_training_energy.png")
save(figure, fig)
@printf("  wrote %s\n", figure)
