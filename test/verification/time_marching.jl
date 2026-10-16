# Check 1 of §5V.P37: the reduced leapfrog of the analytic model at `k = 20` against the closed form,
# over `T = 1000`. It is the `T = 1000` run and belongs to the `slow` group alone: the two marches of
# `10^7` steps are about three minutes, and the check is not run beside the core files.

using LinearAlgebra
using ParticleMethods: ParticleList
using Random
using ReducedBasisMethods: ReducedBasis, CotangentLiftSVD
using Test
using VlasovMethods
using VlasovMethods: Snapshots

include("../helpers/analytic_field.jl")
include("../helpers/analytic_vlasov.jl")

@testset "check 1: the reduced march at k = 20 against the closed form" begin
    T = Float64
    β, dt = T(ANALYTIC_β), T(1e-4)

    # the leapfrog takes `nₜ = 10^7` steps and `nₛ = 10001` saved states, so `nsave` — the stride —
    # is 1000 and the reduced model saves every `0.1`
    nₜ, nₛ = 10^7, 10^4 + 1
    nsave = div(nₜ, nₛ - 1)
    @test nsave == 1000
    @test (nₛ - 1) * nsave == nₜ

    n, k = 1000, 20
    q0, p0 = analytic_initial_condition(T; n = n)

    # the draws are `Xoshiro`, as §5V.P37 states: the seeded global RNG gives the same values, so
    # the test's initial condition is §4's
    @test q0 == T(10) .* (Random.seed!(1); randn(n))

    # the basis is trained on the six training `β`, each saved every `1.0` over `[0, 1000]`
    S = training_snapshots(T, q0, p0)
    b = ReducedBasis(
        CotangentLiftSVD(), S, NamedTuple(), training_parameter_space(T); k = k)
    Φ = b.Φ
    @test size(Φ) == (n, k)

    # the model's field, in the one basis that lifts both the positions and the velocities
    efield = ReducedElectricField(QuadraticField(β), Φ, Φ' * q0)
    plt = ParticleList(reshape(q0, 1, n), reshape(p0, 1, n), ones(T, 1, n))
    SS = Snapshots(T, 1, n, 0, nₛ, test_parameter_space(T))
    IP = IntegratorParameters(dt, nₜ, nₛ, 0, n, 1)
    IC = ReducedIntegratorCache(IP, k)

    @test reduced_integrate_vp(plt, Φ, (χ = T(1),), efield, SS, IP, IC, 1) === nothing

    # the state of the last save is a computed one: a save decided on the save counter writes the
    # initial state and nothing after it, and leaves every later state at zero
    @test norm(SS.V[1, :, end, 1]) > 0

    """
    The relative `L²` error of the saved states of `snapshots` against the closed form, as one global
    ratio `√(Σⱼ‖xⱼ - x̄ⱼ‖²) / √(Σⱼ‖x̄ⱼ‖²)` over the 10001 states, with `x = [q; p]` unweighted. The
    reduced save `j` is paired with the closed form at the time of step `1000j`, with the integer
    product formed first — the ratio 10 of the full model's spacing of `0.01` to the reduced `0.1` is
    that integer step count, and no floating time is ever matched.
    """
    function relative_error(snapshots)
        num, den = zero(T), zero(T)
        for j in 0:(nₛ - 1)
            t = T(nsave * j) * dt
            Q, P = analytic_solution(T, β, q0, p0, (t,))
            x = [snapshots.X[1, :, j + 1, 1]; snapshots.V[1, :, j + 1, 1]]
            num += sum(abs2, x - [Q[:, 1]; P[:, 1]])
            den += sum(abs2, [Q[:, 1]; P[:, 1]])
        end
        return sqrt(num) / sqrt(den)
    end

    err20 = relative_error(SS)

    # the bound is one term per error source: the `k = 20` basis truncation, the leapfrog's phase
    # error and its step truncation. The phase error is `β·T·(βdt)²/(24√3) = 5.196e-5`, and it is the
    # whole of the measured error at every `k` from 2 to 20; the truncation term `σ₂₁/σ₁` is about
    # `1e-17`, because the analytic data has rank 2 under the cotangent lift.
    @test b.Λ[1] > 0
    truncation = sqrt(b.Λ[k + 1] / b.Λ[1])
    bound = truncation +
            β * (T(nₜ) * dt) * (β * dt)^2 / (24 * sqrt(T(3))) * (1 + 1 / 10^4) +
            (β * dt)^2
    @test err20 ≤ bound

    # the basis truncation alone is far below the error, which is why the bound is composite
    @test truncation < err20

    # the error at `k = 20` is below the error at `k = 1`, where the basis is below the data's rank
    # and the truncation error `σ₂/σ₁ ≈ 3.46e-2` dominates. A `k = 5` comparison is not discriminating:
    # both `k = 5` and `k = 20` resolve the rank-2 data, so round-off decides their error ordering.
    fill!(SS.X, zero(T))
    fill!(SS.V, zero(T))
    b1 = ReducedBasis(
        CotangentLiftSVD(), S, NamedTuple(), training_parameter_space(T); k = 1)
    IC1 = ReducedIntegratorCache(IP, 1)
    efield1 = ReducedElectricField(QuadraticField(β), b1.Φ, b1.Φ' * q0)
    reduced_integrate_vp(plt, b1.Φ, (χ = T(1),), efield1, SS, IP, IC1, 1)
    err1 = relative_error(SS)

    @test err20 < err1
end
