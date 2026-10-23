# The analytic quadratic model of the paper's second example: the closed-form motion of `n`
# particles in the harmonic potential `φ = β²x²/2`, the two bases of its snapshots, and the two
# reduced forms of its equation — the ODE form in the POD basis integrated by RK4, and the PODE
# form in the cotangent-lift basis integrated by Lobatto IIIA–IIIB(2).
#
# The paper's command, from `reference/Vlasov_Analytic.jl`, `Reduced_Vlasov_ODE.jl` and
# `Reduced_Vlasov_PODE_CotangentLift.jl`: the full model at `dt = 0.001` to `T = 1000`, and the
# reduced model at `dt = 0.0001` to `T = 1000` saved every 1000th step, which is the 10 000 saved
# states of the paper, with `k = 5` for the ODE/RK4 form and `k = 20` for the PODE form. The
# default is smaller — `T = 20` at `dt = 0.001` — so that a reproduction takes about a minute
# rather than hours:
#
#   julia --project=scripts scripts/analytic_quadratic.jl --final-time 1000 --dt 0.001 \
#       --dt-reduced 0.0001 --nsave 10 --nsave-reduced 1000
#
#   julia --project=scripts scripts/analytic_quadratic.jl [--final-time 20] [--dt 0.001] \
#       [--dt-reduced 0.001] [--nsave 10] [--nsave-reduced 10] [--runs-dir dir] \
#       [--results-dir dir]
#
# The three files it writes go to the runs directory: `analytic_full.h5` with the closed-form
# states and the root attribute `dt`, and `analytic_ode_k5.h5` and `analytic_pode_k20.h5` with the
# reduced states of the two forms. It also draws the spectrum of the snapshots (item 11).

include(joinpath(@__DIR__, "options.jl"))

using Printf

using CairoMakie
using GeometricIntegrators: LobattoIIIAIIIB, RK4, integrate
using HDF5
using LinearAlgebra
using Random
using ReducedBasisMethods: CotangentLiftEVD, PODSVD, ReducedBasis, project, reduced_ode,
                           reduced_pode, spectrumplot
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods

# The parameters of §4's table, as `reference/Vlasov_Analytic.jl` draws them.
const BETA = 6.0            # the frequency of the harmonic potential
const SIGMA_X = 10.0        # spread of the initial positions
const SIGMA_V = 0.5         # spread of the bump's velocities
const V0 = 4.0              # offset velocity of the bump
const A_TAIL = 0.3          # fraction of the particles in the bump
const NPART = 1000          # particles of the full model
const SEED_Q, SEED_P, SEED_TAIL = 1, 10, 100
const K_ODE = 5             # the basis size of the ODE/RK4 form
const K_PODE = 20           # the basis size of the PODE form

const OPT = parse_arguments(
    ARGS; final_time = 20.0, dt = 1e-3, dt_reduced = 1e-3, nsave = 10,
    nsave_reduced = 10)

"The initial positions and velocities of §4, drawn from the three seeds of the reference."
function analytic_initial_condition(n)
    q0 = SIGMA_X .* randn(Xoshiro(SEED_Q), n)
    p0 = randn(Xoshiro(SEED_P), n)

    rng = Xoshiro(SEED_TAIL)
    for i in 1:n
        rand(rng) >= 1 / (1 + A_TAIL) && (p0[i] = SIGMA_V * p0[i] + V0)
    end

    return q0, p0
end

"The closed form of the reference: `q(t) = q₀cos(βt) + (p₀/β)sin(βt)` and
`p(t) = p₀cos(βt) - βq₀sin(βt)`, as the two blocks of one `2n × length(ts)` matrix."
function analytic_states(q0, p0, ts)
    n = length(q0)
    S = Matrix{Float64}(undef, 2n, length(ts))
    for (j, t) in enumerate(ts)
        S[1:n, j] .= q0 .* cos(BETA * t) .+ (p0 ./ BETA) .* sin(BETA * t)
        S[(n + 1):2n, j] .= p0 .* cos(BETA * t) .- (BETA .* q0) .* sin(BETA * t)
    end
    return S
end

"The reduced states `z` of `b` over `nchunks` chunks of `nsave` steps of `dt`, one state per
chunk, the initial state first: the form `step` is `reduced_ode` with `RK4`, or `reduced_pode` with
Lobatto IIIA–IIIB(2). The model is autonomous, so each chunk restarts the clock at zero and the
physical time is tracked as the chunk index times the chunk's span."
function reduced_states(form, b, h, z0, nchunks, nsave, dt)
    n = length(z0)
    Z = Matrix{Float64}(undef, n, nchunks + 1)
    Z[:, 1] .= z0
    z = z0
    for j in 1:nchunks
        span = (0.0, nsave * dt)
        if form === :ode
            sol = integrate(reduced_ode(b, h, span, dt, z), RK4())
            z = sol.q[end]
        else
            sol = integrate(reduced_pode(b, h, span, dt, z), LobattoIIIAIIIB(2);
                f_abstol = 1e-14)
            z = vcat(sol.q[end], sol.p[end])
        end
        Z[:, j + 1] .= z
    end
    return Z
end

q0, p0 = analytic_initial_condition(NPART)

# the full model: the closed form at every saved time
nt = round(Int, OPT.final_time / OPT.dt)
nchunks = div(nt, OPT.nsave)
t = OPT.nsave * OPT.dt .* (0:nchunks)
S = analytic_states(q0, p0, t)

fpath = runs_path("analytic_full.h5")
h5open(fpath, "w") do io
    io["t"] = collect(t)
    io["q"] = S
    HDF5.attributes(io)["dt"] = OPT.dt
end
@printf("analytic_quadratic: %d particles, β = %g, dt = %g, t ∈ [0, %g]\n", NPART, BETA,
    OPT.dt, OPT.final_time)
@printf("  wrote %s (%d saved states)\n", fpath, length(t))

# the bases of the snapshots: the POD basis of the ODE form and the cotangent-lift basis of the
# PODE form, both from the stacked `[q; p]` snapshot matrix
space = ParameterSpace(Parameter(:β, BETA, BETA, 1))
parameters = (β = BETA, σx = SIGMA_X, σv = SIGMA_V, v₀ = V0, a = A_TAIL)
b_ode = ReducedBasis(PODSVD(), S, parameters, space; k = K_ODE)
b_pode = ReducedBasis(CotangentLiftEVD(), S, parameters, space; k = K_PODE)
h = AnalyticHamiltonian(BETA)

# item 11: the spectrum of the snapshots
spectrum = spectrumplot(b_ode.Λ; ylabel = "λᵢ")
ipath = results_path("item11_analytic_spectrum.png")
save(ipath, current_figure())
@printf("  wrote %s\n", ipath)

# the reduced model, in the two forms
nt_reduced = round(Int, OPT.final_time / OPT.dt_reduced)
nchunks_reduced = div(nt_reduced, OPT.nsave_reduced)
t_reduced = OPT.nsave_reduced * OPT.dt_reduced .* (0:nchunks_reduced)

Z_ode = reduced_states(:ode, b_ode, h, project(b_ode, vcat(q0, p0)), nchunks_reduced,
    OPT.nsave_reduced, OPT.dt_reduced)
fpath = runs_path("analytic_ode_k$(K_ODE).h5")
h5open(fpath, "w") do io
    io["t"] = collect(t_reduced)
    io["z"] = Z_ode
    HDF5.attributes(io)["dt"] = OPT.dt_reduced
end
@printf("  wrote %s (%d saved states, k = %d)\n", fpath, length(t_reduced), K_ODE)

Z_pode = reduced_states(:pode, b_pode, h, project(b_pode, vcat(q0, p0)), nchunks_reduced,
    OPT.nsave_reduced, OPT.dt_reduced)
fpath = runs_path("analytic_pode_k$(K_PODE).h5")
h5open(fpath, "w") do io
    io["t"] = collect(t_reduced)
    io["z"] = Z_pode
    HDF5.attributes(io)["dt"] = OPT.dt_reduced
end
@printf("  wrote %s (%d saved states, k = %d)\n", fpath, length(t_reduced), K_PODE)
