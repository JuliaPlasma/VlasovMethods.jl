# The analytic quadratic-potential example: the exact solution of the model of the acceptance
# list, the reduced basis of its training trajectories, and the reduced march at the test
# parameter.
#
#   julia --project=scripts scripts/analytic_quadratic.jl [--npart 200] [--dt 1e-3] \
#       [--rdt 1e-3] [--final-time 10.0] [--save-step 10] [--k 20] [--runs-dir dir] \
#       [--results-dir dir]
#
# The data go to `--runs-dir`, the figures to `--results-dir`. The paper's parameters are
# `--npart 1000 --final-time 1000 --dt 1e-3 --rdt 1e-4`; the defaults are smaller, so that the
# run stays inside a few minutes.

include(joinpath(@__DIR__, "options.jl"))
# the shared file of the scripts, for the commit of the checkout
include(joinpath(@__DIR__, "bump_on_tail_common.jl"))

using CairoMakie
using HDF5
using LinearAlgebra
using Printf
using Random: Xoshiro, randn, rand
using ReducedBasisMethods: CotangentLiftSVD, ReducedBasis, lift, project, reduced_pode
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods
using VlasovMethods: AnalyticHamiltonian, plot_series

using GeometricIntegrators: LobattoIIIAIIIB, integrate

const OPT = parse_arguments(ARGS; npart = 200, dt = 1e-3, rdt = 1e-3, final_time = 10.0,
    save_step = 10, k = 20)

const NPART = OPT.npart
const DT = OPT.dt
const RDT = OPT.rdt
const FINAL_TIME = OPT.final_time
const SAVE_STEP = OPT.save_step
const K = OPT.k

# the parameters of the model: the six training values of β and the test one, which lies off the
# training range
const TRAINING_β = (5.95, 5.97, 5.99, 6.01, 6.03, 6.05)
const TEST_β = 6.00

const RUNID = "analytic_beta600_n$(NPART)_T$(replace(string(FINAL_TIME), "." => "p"))"

"""
    initial_condition(::Type{T}, n)

The initial condition of the acceptance list: `q₀ = 10·randn(n)` from the generator of seed 1,
`p₀ ~ N(0, 1)` from the generator of seed 10, and then, for every particle whose draw of the
generator of seed 100 is `≥ 1/(1 + a)`, the bump `p₀ ← σᵥ·p₀ + v₀` with `a = 0.3`, `σᵥ = 0.5` and
`v₀ = 4.0`.
"""
function initial_condition(::Type{T}, n) where {T}
    q0 = T(10) .* randn(Xoshiro(1), T, n)
    p0 = randn(Xoshiro(10), T, n)

    a, σv, v0 = T(0.3), T(0.5), T(4)
    rng = Xoshiro(100)
    for i in 1:n
        rand(rng) >= 1 / (1 + a) && (p0[i] = σv * p0[i] + v0)
    end

    return q0, p0
end

"""
    analytic_solution(::Type{T}, β, q₀, p₀, ts)

The closed form of the model at the times `ts`, as two `n × length(ts)` matrices: the positions
`q(t) = q₀ cos(βt) + (p₀/β) sin(βt)` and the velocities `p(t) = p₀ cos(βt) - βq₀ sin(βt)`. The
model is `q̇ = p`, `ṗ = -β²q`, and its solution has no discretisation error.
"""
function analytic_solution(::Type{T}, β, q0, p0, ts) where {T}
    Q = reduce(hcat, [(p0 ./ β) .* sin(β * t) .+ q0 .* cos(β * t) for t in ts])
    P = reduce(hcat, [p0 .* cos(β * t) .- (β .* q0) .* sin(β * t) for t in ts])
    return Q, P
end

const Q0, P0 = initial_condition(Float64, NPART)

# the saved times: the closed form is evaluated at them only, and never at every step
const TS = collect(range(0.0, FINAL_TIME; step = SAVE_STEP * DT))

# the training snapshot: the stacked `[Q; P]` of the six training wave numbers
const S = Matrix{Float64}(undef, 2 * NPART, length(TRAINING_β) * length(TS))
for (j, β) in enumerate(TRAINING_β)
    Q, P = analytic_solution(Float64, β, Q0, P0, TS)
    columns = ((j - 1) * length(TS) + 1):(j * length(TS))
    S[1:NPART, columns] .= Q
    S[(NPART + 1):(2 * NPART), columns] .= P
end

const PSC = ParameterSpace(Parameter(:β, 5.95, 6.05, length(TRAINING_β)))
const B = ReducedBasis(CotangentLiftSVD(), S, NamedTuple(), PSC; k = K)

@printf("the basis holds %d of %d columns\n", B.k, size(B.Φ, 1))

# the reduced march at the test value of β, which lies off the training range
const H = AnalyticHamiltonian(TEST_β)
const X0 = vcat(Q0, P0)
const PROB = reduced_pode(B, H, (0.0, FINAL_TIME), RDT, project(B, X0))
const SOL = integrate(PROB, LobattoIIIAIIIB(2); f_abstol = 1e-14)

# The columns of the reduced trajectory that hold the saved times: the march saves every step of
# `RDT`, and the closed form is evaluated every `SAVE_STEP` steps of `DT`. Only these columns are
# lifted, so that the memory of the run does not grow with the number of steps.
const QRED = Array(SOL.q)
const PRED = Array(SOL.p)
const STRIDE = round(Int, SAVE_STEP * DT / RDT)
const COLUMNS = 1:STRIDE:size(QRED, 2)
const LIFTED = lift(B, vcat(QRED[:, COLUMNS], PRED[:, COLUMNS]))
const QL = LIFTED[1:NPART, :]
const PL = LIFTED[(NPART + 1):(2 * NPART), :]

@printf("the reduced march took %d steps\n", size(QRED, 2) - 1)

# the comparison against the closed form, at the saved times only
const QEXACT, PEXACT = analytic_solution(Float64, TEST_β, Q0, P0, TS)
length(COLUMNS) == length(TS) || error(
    "the reduced march holds $(length(COLUMNS)) states at the saved times; the closed form is " *
    "evaluated at $(length(TS))")

const ERRORS = trajectory_errors(QEXACT, PEXACT, QL, PL)

# the total energy of the two trajectories: `H = ½(p² + β²q²)`, summed over the particles
const HEXACT = [(sum(abs2, PEXACT[:, t]) .+ TEST_β^2 * sum(abs2, QEXACT[:, t])) / 2
                for t in axes(QEXACT, 2)]
const HREDUCED = [(sum(abs2, PL[:, t]) .+ TEST_β^2 * sum(abs2, QL[:, t])) / 2
                  for t in axes(QL, 2)]

const COMMIT = repository_commit()

h5open(runs_path("$(RUNID).h5"), "w") do file
    file["t"] = TS
    file["Q"] = QEXACT
    file["P"] = PEXACT
    file["Q_reduced"] = QL
    file["P_reduced"] = PL
    HDF5.attributes(file)["runid"] = RUNID
    HDF5.attributes(file)["commit"] = COMMIT
    HDF5.attributes(file)["command"] = join(ARGS, " ")
    HDF5.attributes(file)["seed"] = 1
end

# the figures: the errors of the reduced trajectory against the closed form (item 24) and the
# drift of the total energy (item 27)
fig = plot_series(TS, (ERRORS.x2, ERRORS.xmax, ERRORS.v2, ERRORS.vmax);
    labels = ["x2", "xmax", "v2", "vmax"], xlabel = "t", ylabel = "error")
save(results_path("$(RUNID)_errors.png"), fig)

fig = plot_series(TS, (energy_drift(HEXACT), energy_drift(HREDUCED));
    labels = ["closed form", "reduced"], xlabel = "t", ylabel = "ΔH/H")
save(results_path("$(RUNID)_energy_drift.png"), fig)
