# Check 3: the paper's central claim. The PODE form of the reduced equation — the PSD cotangent basis
# at `k = 20` integrated by Lobatto IIIA–IIIB(2) — conserves the analytic model's Hamiltonian to
# `O(dt²)`. This file runs it to `T = 100` and belongs to the `slow` group: the PODE's `10^6` steps
# take over a minute. The `T = 1000` run of both reduced forms, the PODE form and the ODE
# form in the POD basis at `k = 5` integrated by RK4, is the paper command of
# `scripts/analytic_quadratic.jl`, which writes the bound, the largest drift and the first save
# that exceeds the bound of each form; it trains on the closed form at its own saved times, not on
# this file's training set.
#
# The run advances in chunks of `1000` steps and keeps one state per chunk — the `0.1` grid of the
# saves, `1001` of them, where every step would hold `10^6` states. Each chunk is one emitted
# problem over the timespan `chunk_span`, and **what is asserted is the number of steps its
# solution takes, never the problem's own count**: the plain span `(0, T(1000)·dt)` takes exactly
# `1000` steps, while `ntime(prob)` counts `1001` and the span rounded down by one `eps` takes
# `999`. A chunk built from the run's own clock — the float span `(t_{j-1}, t_j)` — takes `1000` at
# some `j` and `999` at others, which puts every later save at the wrong time. The chunk's equation
# is autonomous — the PODE form does not read `t` — so the clock starts at zero each time.

using LinearAlgebra
using ReducedBasisMethods: ReducedBasis, CotangentLiftSVD, project, lift, reduced_pode
using Test
using VlasovMethods
# an explicit only-list, and not `using GeometricIntegrators`: both it and `VlasovMethods` export
# `GeometricIntegrator`, and the only-list binds it to this one
using GeometricIntegrators: LobattoIIIAIIIB, integrate, ntime

include("../helpers/analytic_vlasov.jl")
include("../helpers/reduced_runs.jl")

@testset "check 3: the PODE form conserves H" begin
    T = Float64
    β, dt = T(ANALYTIC_β), T(1e-4)
    nsave, nₜ, nₛ = 1000, 10^6, 10^3 + 1
    n, kp = 1000, 20
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
    # `10^3` chunks are the `10^6` steps of the horizon. `ntime(prob)` is not this count — it is
    # `div(tend - tbegin, tstep, RoundUp)`, which is `1001` here — and the span rounded down by one
    # `eps` takes `999`, as the float span `(t_{j-1}, t_j)` takes at some `j`, so this is where a
    # chunk that falls a step short is caught and the saves after it are known to be at the right
    # time
    @test all(==(nsave), counts)
    @test length(counts) == nchunks
    @test sum(counts) == nₜ
    @test all(isfinite, drift)
    @test maximum(drift) ≤ bound
end
