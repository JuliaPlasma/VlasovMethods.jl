# The analytic quadratic model of the paper: six training β, the closed-form trajectories of the
# quadratic potential `E(x) = β²x`, the cotangent-lift and POD bases of them, and the reduced ODE
# and PODE marches at the off-sample β = 6.00.
#
#   julia --project=scripts scripts/analytic_quadratic.jl [--particles 1000] [--final-time 100.0] \
#       [--tstep 1e-4] [--k-ode 5] [--k-pode 20] [--runs-dir dir] [--results-dir dir]
#
# The default stops at `T = 100` of the paper's `T = 1000`, at the same `dt = 1e-4` and with the
# same six β, so that a run is minutes rather than hours; the paper's run is
#
#   julia --project=scripts scripts/analytic_quadratic.jl --final-time 1000.0
#
# The closed form is evaluated at the saved times only, never at every step, so no trajectory of
# `2n × 10⁶` states is stored. The figures go to `--results-dir`: the spectrum of the training
# basis, the errors against the closed form and the drift of the Hamiltonian of the two forms.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using CairoMakie
using Random
using ReducedBasisMethods: CotangentLiftSVD, PODSVD, ReducedBasis, lift, project,
                           reduced_ode, reduced_pode, spectrumplot
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods
using GeometricIntegrators: LobattoIIIAIIIB, RK4, integrate

const OPT = parse_arguments(ARGS; particles = 1000, final_time = 100.0, tstep = 1e-4,
    k_ode = 5, k_pode = 20)

const T = Float64
const DT = OPT.tstep
const TSPAN = OPT.final_time

# the six training β of §4 and the off-sample test β, on which the reduced model is tested
const TRAINING_β = (5.95, 5.97, 5.99, 6.01, 6.03, 6.05)
const TEST_β = 6.0

# the spacing of the training snapshot and of the reduced saves: the paper saves every 10th step of
# the closed form, and the reduced model every 1000th step of its `10⁷`
const SAVE_TIME = 1.0
const NSAVE = 1000

const NT = round(Int, TSPAN / DT)
const NS = div(NT, NSAVE) + 1

"""
    initial_condition(::Type{T}; n) -> (q₀, p₀)

The initial condition of §4: `q₀ = 10·randn(n)` under seed 1, `p₀ ~ N(0,1)` under seed 10, then
`p₀ ← 0.5·p₀ + 4.0` for the particles whose seeded draw of seed 100 is `≥ 1/(1 + 0.3)`.
"""
function initial_condition(::Type{T}; n) where {T}
    q0 = T(10) .* randn(Xoshiro(1), T, n)
    p0 = randn(Xoshiro(10), T, n)
    rng = Xoshiro(100)
    for i in 1:n
        rand(rng) >= 1 / (1 + T(0.3)) && (p0[i] = T(0.5) * p0[i] + T(4))
    end
    return q0, p0
end

"`q(t) = q₀cos(βt) + (p₀/β)sin(βt)` and `p(t) = p₀cos(βt) - βq₀sin(βt)`, one column per time."
function closed_form(::Type{T}, β, q0, p0, ts) where {T}
    Q = reduce(hcat, [(p0 ./ β) .* sin(β * t) .+ q0 .* cos(β * t) for t in ts])
    P = reduce(hcat, [p0 .* cos(β * t) .- (β .* q0) .* sin(β * t) for t in ts])
    return Q, P
end

"The stacked `[Q; P]` training snapshot of the six β, one column per saved time and β."
function training_snapshot(::Type{T}, q0, p0; step = SAVE_TIME) where {T}
    n = length(q0)
    ts = T(0):T(step):T(TSPAN)
    S = Matrix{T}(undef, 2n, length(TRAINING_β) * length(ts))
    for (j, β) in enumerate(TRAINING_β)
        Q, P = closed_form(T, T(β), q0, p0, ts)
        cols = ((j - 1) * length(ts) + 1):(j * length(ts))
        S[1:n, cols] .= Q
        S[(n + 1):(2n), cols] .= P
    end
    return S
end

q0, p0 = initial_condition(T; n = OPT.particles)
h = AnalyticHamiltonian(T(TEST_β))
tps = ParameterSpace(Parameter(:β, T(TRAINING_β[1]), T(TRAINING_β[end]), length(TRAINING_β)))

S = training_snapshot(T, q0, p0)
bc = ReducedBasis(CotangentLiftSVD(), S, NamedTuple(), tps; k = OPT.k_pode)
bp = ReducedBasis(PODSVD(), S, NamedTuple(), tps; k = OPT.k_ode)
z0 = project(bc, [q0; p0])
zp = project(bp, [q0; p0])

@printf("analytic_quadratic: %d particles, 6 training β, dt = %g, t ∈ [0, %g]\n",
    OPT.particles, DT, TSPAN)
@printf("  %d saved states, k = %d (PODE) and k = %d (ODE)\n", NS, bc.k, bp.k)

# the PODE form of the cotangent-lift basis, one chunk of `NSAVE` steps per saved state, at the
# solver tolerance of the verification run
Zc = Matrix{T}(undef, 2bc.k, NS)
Zc[:, 1] = z0
for j in 1:(NS - 1)
    prob = reduced_pode(bc, h, (zero(T), T(NSAVE) * DT), DT, Zc[:, j])
    sol = integrate(prob, LobattoIIIAIIIB(2); f_abstol = 1e-14)
    Zc[:, j + 1] = vcat(sol.q[end], sol.p[end])
end

# the ODE form of the POD basis, the Galerkin form of `k = 5`
Zp = Matrix{T}(undef, bp.k, NS)
Zp[:, 1] = zp
for j in 1:(NS - 1)
    prob = reduced_ode(bp, h, (zero(T), T(NSAVE) * DT), DT, Zp[:, j])
    Zp[:, j + 1] = integrate(prob, RK4()).q[end]
end

# the lifted reduced trajectories, and the closed form at the saved times
Xc = lift(bc, Zc)
Xp = lift(bp, Zp)
ts = [T(NSAVE * j) * DT for j in 0:(NS - 1)]
Q, P = closed_form(T, T(TEST_β), q0, p0, ts)

# item 24: the errors of the lifted PODE trajectory against the closed form
n = length(q0)
err = trajectory_errors(Q, P, Xc[1:n, :], Xc[(n + 1):(2n), :])

# item 27: the Hamiltonian of the lifted state of each form, and its drift
energy_at(X) = [h(zero(T), X[:, j], NamedTuple()) for j in axes(X, 2)]
Hc = energy_at(Xc)
Hp = energy_at(Xp)

figure = results_path("analytic_quadratic_spectrum.png")
save(figure, spectrumplot(bc.Λ; ylabel = "λᵢ", title = "training basis, k = $(bc.k)"))
@printf("  wrote %s\n", figure)

figure = results_path("analytic_quadratic_errors.png")
save(figure,
    plot_series(ts, (err.x2, err.xmax, err.v2, err.vmax);
        labels = ["RMS x", "L∞ x", "RMS v", "L∞ v"], xlabel = "t", ylabel = "error"))
@printf("  wrote %s\n", figure)

figure = results_path("analytic_quadratic_energy_drift.png")
save(figure,
    plot_series(ts, (energy_drift(Hp), energy_drift(Hc));
        labels = ["ODE, k = $(bp.k)", "PODE, k = $(bc.k)"], xlabel = "t", ylabel = "ΔH/H"))
@printf("  wrote %s\n", figure)
