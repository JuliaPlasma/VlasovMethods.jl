# The analytic quadratic-potential example: `n` particles in the external field `E(x) = β²x`,
# whose exact solution is known, so the full model has no discretisation error. The bases are
# trained on the closed form at six `β ∈ [5.95, 6.05]` and tested off-sample at `β = 6.0`: the POD
# basis of `k = 5` with the reduced ODE form and RK4, and the cotangent-lift basis of `k = 20` with
# the reduced PODE form and Lobatto IIIA–IIIB(2), both at the reduced step `rtstep`.
#
#   julia --project=scripts scripts/analytic_quadratic.jl [--runid analytic_quadratic] \
#       [--npart 1000] [--beta 6.0] [--tstep 0.001] [--train-step 1000] [--save-step 100] \
#       [--rtstep 0.0001] [--final-time 100.0] [--kode 5] [--kpode 20] \
#       [--runs-dir dir] [--results-dir dir]
#
# `tstep` is the step of the full model, which the closed form needs only to place its saves:
# the training snapshots are `train_step` steps apart, and the states that the reduced runs are
# compared at are `save_step` steps apart. The paper trains on every 10th step up to `T = 1000`, a
# snapshot of 2000 × 600006 entries, 9.6 GB, and its PODE run alone takes 10⁷ reduced steps. The
# default trains on every 1000th step and runs to `T = 100`, and is otherwise the paper's run. The
# paper's parameters are
#
#   julia --project=scripts scripts/analytic_quadratic.jl --train-step 10 --final-time 1000.0
#
# The run writes `<runid>.h5` to `--runs-dir`, and to `--results-dir` the errors of the positions
# and the velocities (item 24), the drift of the energy (27) and the spectra of the two bases (11).

include(joinpath(@__DIR__, "options.jl"))
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using Printf

using CairoMakie
using HDF5
using ReducedBasisMethods: ReducedBasis, PODSVD, CotangentLiftSVD, project, lift,
                           reduced_ode,
                           reduced_pode, spectrumplot
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods
# an only-list: both GeometricIntegrators and VlasovMethods export `GeometricIntegrator`
using GeometricIntegrators: LobattoIIIAIIIB, RK4, integrate

const OPT = parse_arguments(ARGS; runid = "analytic_quadratic", npart = 1000, beta = 6.0,
    tstep = 0.001, train_step = 1000, save_step = 100, rtstep = 1e-4, final_time = 100.0,
    kode = 5, kpode = 20)

# the six training `β` of the paper
const TRAINING_BETA = [5.95, 5.97, 5.99, 6.01, 6.03, 6.05]

"""
The initial condition of the paper: `q₀ = 10 randn(n)` from the seed 1; `p₀ = randn(n)` from the
seed 10, and `p₀ ← 0.5 p₀ + 4` for each particle whose draw from the seed 100 is at least
`1 / (1 + 0.3)`.
"""
function initial_condition(n)
    q₀ = 10 .* randn(Xoshiro(1), n)
    p₀ = randn(Xoshiro(10), n)
    rng = Xoshiro(100)
    for i in 1:n
        rand(rng) >= 1 / (1 + 0.3) && (p₀[i] = 0.5 * p₀[i] + 4.0)
    end
    return q₀, p₀
end

"The closed form at the time `t`: `q₀ cos βt + p₀ sin(βt)/β` and `p₀ cos βt − β q₀ sin βt`."
function exact(β, q₀, p₀, t)
    return q₀ .* cos(β * t) .+ (p₀ ./ β) .* sin(β * t),
    p₀ .* cos(β * t) .- (β .* q₀) .* sin(β * t)
end

n = OPT.npart
q₀, p₀ = initial_condition(n)

# the training snapshot `[Q; P]` of the six `β`, each saved every `train_step` steps
ttrain = 0:(OPT.train_step * OPT.tstep):(OPT.final_time)
S = reduce(hcat, [vcat(exact(β, q₀, p₀, t)...) for β in TRAINING_BETA for t in ttrain])
paramspace = ParameterSpace(Parameter(:β, TRAINING_BETA))
pod = ReducedBasis(PODSVD(), S, NamedTuple(), paramspace; k = OPT.kode)
psd = ReducedBasis(CotangentLiftSVD(), S, NamedTuple(), paramspace; k = OPT.kpode)
@printf("analytic_quadratic: %d training snapshots, k = %d (ODE), k = %d (PODE)\n",
    size(S, 2), pod.k, psd.k)

# the saves of the reduced runs: `stride` reduced steps of `rtstep` each, `save_step` full steps
stride = save_stride(OPT.save_step, OPT.tstep, OPT.rtstep)
Δs = OPT.save_step * OPT.tstep
nsaves = round(Int, OPT.final_time / Δs)
t = (0:nsaves) .* Δs
span = (0.0, stride * OPT.rtstep)
h = AnalyticHamiltonian(OPT.beta)

"""
The errors of the reduced run in the basis `b` against the closed form at the times `t`: the step
`advance(z)` takes the reduced state from one save to the next. A run whose state stops being
finite ends there, and its later entries stay `NaN`.
"""
function compare(b, z, advance)
    errors = (x2 = fill(NaN, nsaves + 1), xmax = fill(NaN, nsaves + 1),
        v2 = fill(NaN, nsaves + 1), vmax = fill(NaN, nsaves + 1))
    H = fill(NaN, nsaves + 1)
    for j in eachindex(t)
        j > 1 && (z = advance(z))
        all(isfinite, z) || break
        x = lift(b, z)
        qe, pe = exact(OPT.beta, q₀, p₀, t[j])
        e = trajectory_errors(reshape(qe, n, 1), reshape(pe, n, 1),
            reshape(x[1:n], n, 1), reshape(x[(n + 1):(2n)], n, 1))
        for k in keys(errors)
            errors[k][j] = only(e[k])
        end
        H[j] = h(0.0, x, NamedTuple())
    end
    return errors, H
end

function advance_ode(z)
    sol = integrate(reduced_ode(pod, h, span, OPT.rtstep, z), RK4())
    return sol.q[end]
end

function advance_pode(z)
    sol = integrate(reduced_pode(psd, h, span, OPT.rtstep, z), LobattoIIIAIIIB(2);
        f_abstol = 1e-14)
    return vcat(sol.q[end], sol.p[end])
end

time_ode = @elapsed errors_ode, H_ode = compare(pod, project(pod, [q₀; p₀]), advance_ode)
@printf("  ODE, RK4: %.1f s\n", time_ode)
time_pode = @elapsed errors_pode, H_pode = compare(psd, project(psd, [q₀; p₀]), advance_pode)
@printf("  PODE, Lobatto IIIA–IIIB(2): %.1f s\n", time_pode)

drift_ode = energy_drift(H_ode)
drift_pode = energy_drift(H_pode)

fpath = runs_path(OPT.runid * ".h5")
h5open(fpath, "w") do file
    file["t"] = collect(t)
    for (name, errors, H) in (("ode", errors_ode, H_ode), ("pode", errors_pode, H_pode))
        for k in keys(errors)
            file["$(name)/$(k)"] = errors[k]
        end
        file["$(name)/H"] = H
    end
    HDF5.attributes(file)["commit"] = repository_commit()
    HDF5.attributes(file)["command"] = command_line("analytic_quadratic.jl", ARGS)
end
@printf("  wrote %s\n", fpath)

function write_figure(fig, name; logscale = false)
    logscale && (fig.content[1].yscale = log10)
    figure = results_path(OPT.runid * "_" * name * ".png")
    save(figure, fig)
    @printf("  wrote %s\n", figure)
end

# item 24: the errors of the positions and the velocities against the closed form, from the first
# save on, because the errors at `t = 0` are those of the projection alone
later = 2:(nsaves + 1)
for (q, name) in ((:x, "positions"), (:v, "velocities"))
    series = [errors[k][later] for errors in (errors_ode, errors_pode)
              for k in (Symbol(q, 2), Symbol(q, :max))]
    labels = ["$(m), $(k)" for m in ("ODE k = $(pod.k)", "PODE k = $(psd.k)")
              for k in ("RMS", "max")]
    write_figure(
        plot_series(t[later], series; labels = labels,
            ylabel = "error of the $(name)"),
        "errors_$(q)"; logscale = true)
end

# item 27: the relative drift of the energy, which the closed form conserves
write_figure(
    plot_series(t[later], (abs.(drift_ode[later]), abs.(drift_pode[later]));
        labels = ["ODE k = $(pod.k)", "PODE k = $(psd.k)"], ylabel = "|ΔH/H|"),
    "energy_drift"; logscale = true)

# item 11: the spectra of the two bases
for (b, name) in ((pod, "pod"), (psd, "cotangent_lift"))
    write_figure(
        spectrumplot(b.Λ; title = "the spectrum of the $(name) basis",
            ylabel = "Λᵢ"),
        "spectrum_" * name)
end
