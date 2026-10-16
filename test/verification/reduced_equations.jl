# Check 3: the paper's central claim. The PODE form of the reduced equation — the PSD cotangent basis
# at `k = 20` integrated by Lobatto IIIA–IIIB(2) — conserves the analytic model's Hamiltonian to
# `O(dt²)`. The ODE form — the POD basis at `k = 5` integrated by RK4 — runs beside it, and the save
# at which it leaves the bound is logged, not asserted. It is the `T = 1000` run and belongs to the
# `slow` group alone: the PODE's `10^7` steps are about thirteen minutes, and the ODE's run stops at
# the chunk where it first leaves the bound, if it does.
#
# The run advances in chunks of `1000` steps and keeps one state per chunk — the `0.1` grid of the
# saves, `10001` of them, where every step would hold `10^7` states (3.2 GB for the `k = 20` PODE).
# Each chunk is one emitted problem over the timespan `chunk_span`, and **what is asserted is the
# number of steps its solution takes, never the problem's own count**: the plain span
# `(0, T(1000)·dt)` takes exactly `1000` steps, while `ntime(prob)` counts `1001` and the span
# rounded down by one `eps` takes `999`. A chunk built from the run's own clock — the float span
# `(t_{j-1}, t_j)` — takes `1000` at some `j` and `999` at others, which puts every later save at the
# wrong time. The chunk's equation is autonomous — neither form reads `t` — so the clock starts at
# zero each time and the test tracks the physical time as `0.1·j` itself.

using LinearAlgebra
using ReducedBasisMethods: ReducedBasis, PODSVD, CotangentLiftSVD, project, lift,
                           reduced_ode,
                           reduced_pode
using Test
using VlasovMethods
# an explicit only-list, and not `using GeometricIntegrators`: both it and `VlasovMethods` export
# `GeometricIntegrator`, and the only-list binds it to this one
using GeometricIntegrators: LobattoIIIAIIIB, RK4, integrate, ntime

include("../helpers/analytic_vlasov.jl")
include("../helpers/reduced_runs.jl")

@testset "check 3: the PODE form conserves H" begin
    T = Float64
    β, dt = T(ANALYTIC_β), T(1e-4)
    nsave, nₜ, nₛ = 1000, 10^7, 10^4 + 1
    n, kp, ko = 1000, 20, 5
    nchunks = nₛ - 1

    # the leading Verlet error for frequency `β`: `(β·dt)² = 3.6e-7`
    bound = (β * dt)^2

    q0, p0 = analytic_initial_condition(T; n = n)
    S = training_snapshots(T, q0, p0)
    h = AnalyticHamiltonian(β)

    # `H` of the lifted state at `j = 0`, which the drifts are relative to; it is positive for §4's
    # draws, so there is no guard on the denominator
    H0(z, b) = h(zero(T), lift(b, z), NamedTuple())

    # the PODE form: the PSD cotangent basis at `k = 20`, integrated by Lobatto IIIA–IIIB(2)
    bc = ReducedBasis(
        CotangentLiftSVD(), S, NamedTuple(), training_parameter_space(T); k = kp)
    zc = project(bc, [q0; p0])
    reference = H0(zc, bc)
    @test reference > zero(T)

    drift, counts = Float64[], Int[]
    for j in 1:nchunks
        prob = reduced_pode(bc, h, (zero(T), chunk_span(T, nsave, dt)), dt, zc)

        sol = integrate(prob, LobattoIIIAIIIB(2); f_abstol = lobatto_abstol(T, kp))
        push!(counts, ntime(sol))
        zc = vcat(sol.q[end], sol.p[end])
        push!(drift, abs((H0(zc, bc) - reference) / reference))
    end

    # the steps taken: every chunk of the plain span takes exactly the `1000` steps it means, and the
    # `10^4` chunks are the `10^7` steps of the horizon. `ntime(prob)` is not this count — it is
    # `div(tend - tbegin, tstep, RoundUp)`, which is `1001` here — and the span rounded down by one
    # `eps` takes `999`, as the float span `(t_{j-1}, t_j)` takes at some `j`, so this is where a
    # chunk that falls a step short is caught and the saves after it are known to be at the right
    # time
    @test all(==(nsave), counts)
    @test length(counts) == nchunks
    @test sum(counts) == nₜ
    @test all(isfinite, drift)
    @test maximum(drift) ≤ bound

    # the ODE form: the POD basis at `k = 5`, integrated by RK4. The run stops at the chunk where it
    # first exceeds the bound, if it does, because a run that has left the bound can overflow to `Inf`
    # by `T = 1000`.
    bp = ReducedBasis(PODSVD(), S, NamedTuple(), training_parameter_space(T); k = ko)
    zp = project(bp, [q0; p0])
    reference_ode = H0(zp, bp)
    @test reference_ode > zero(T)

    # `j_sep` is the first `j` with `!(|ΔH/H| ≤ bound)`, so that a `NaN` counts as an exceedance, and
    # `nothing` when the drift stays inside the bound over the whole horizon. It is logged and not
    # asserted, because the ODE half measures round-off: the drift comes from POD vectors 4 and 5,
    # which are round-off, so whether and when the ODE form leaves the bound depends on the platform.
    # One run measured `j_sep = 568`, `t_sep = 56.8`, with a drift of `4.11e-7` against the bound
    # `3.6e-7` (Julia 1.13.1); `t_sep` moves from `61.0` to `53.1` when only the spacing of the
    # training set changes; and on windows with Julia 1.13.2 the drift stays at `2.3e-12` over all
    # `10^4` chunks.
    drift_ode, counts_ode = Float64[], Int[]
    for j in 1:nchunks
        prob = reduced_ode(bp, h, (zero(T), chunk_span(T, nsave, dt)), dt, zp)

        sol = integrate(prob, RK4())
        push!(counts_ode, ntime(sol))
        zp = sol.q[end]
        push!(drift_ode, abs((H0(zp, bp) - reference_ode) / reference_ode))
        # the run stops at the chunk where the bound is first exceeded
        (isfinite(drift_ode[end]) && drift_ode[end] ≤ bound) || break
    end

    # the chunks the run did take, each of them the plain span's `1000` steps
    @test all(==(nsave), counts_ode)

    # the separation of the ODE form, on the platform the run is on
    j_sep = first_exceedance(drift_ode, bound)
    t_sep = j_sep === nothing ? nothing : 0.1 * j_sep
    @info "check 3: the ODE form" j_sep t_sep max_drift=maximum(drift_ode) bound
end
