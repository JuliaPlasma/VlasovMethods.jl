using ForwardDiff
using LinearAlgebra
using Random
using ReducedBasisMethods: ReducedBasis, PODSVD, CotangentLiftSVD, project,
                           lift, reduced_ode, reduced_pode
using Test
using VlasovMethods
# an explicit only-list, so that no name it exports can collide with one of `VlasovMethods`'
using GeometricIntegrators: LobattoIIIAIIIB, equation, integrate, ntime, parameters

include("../helpers/allocation_barriers.jl")
include("../helpers/analytic_vlasov.jl")
include("../helpers/reduced_runs.jl")
include("../helpers/vlasov_ode_reference.jl")

@testset "check 4: the POD form is the reference's Galerkin form, $T" for T in (Float64, Float32)
    n, k, nsteps = 200, 5, 1000
    dt = T(1e-4)
    q0, p0 = analytic_initial_condition(T; n = n)
    b = ReducedBasis(PODSVD(), training_snapshots(T, q0, p0; step = 10), NamedTuple(),
        training_parameter_space(T); k = k)
    h = AnalyticHamiltonian(T(ANALYTIC_β))

    # `k = 5` is odd, so `J_5` does not exist: the POD form is the Galerkin form `ż = UᵀJ_N∇H(Uz)`,
    # with `J_N` of `2N` rows, and `U = Φ` holds `2N` of them
    @test b.k == k
    @test size(b.Φ) == (2n, k)

    z0 = project(b, [q0; p0])
    prob = reduced_ode(b, h, (zero(T), chunk_span(T, nsteps, dt)), dt, z0)

    # **The step count is asserted on an integrated solution, never on the problem.** The problem
    # here is built and not integrated, so it asserts no count: GeometricEquations' `ntime(prob)` is
    # `div(tend - tbegin, tstep, RoundUp)`, which is not the number of steps the integrator takes
    # over the solution's `tbegin:Δt:tend` grid — it is one larger for this span in `Float64`. Check
    # 3's chunks and check 5's problems carry the assertion on their solutions.

    # the emitted vector field is the reference's to `√eps`, at 100 points from
    # `randn(Xoshiro(3), T, 5)`
    vf = equation(prob).v
    params = parameters(prob)
    Z = randn(Xoshiro(3), T, k, 100)
    ref, v = zeros(T, k), zeros(T, k)
    for j in axes(Z, 2)
        z = Z[:, j]
        vlasov_ode_reference!(ref, zero(T), z, T(ANALYTIC_β), b.Φ)
        @test norm(ref, Inf) > zero(T)
        @test vf(v, zero(T), z, params) === v
        @test norm(v - ref, Inf) ≤ √eps(T) * norm(ref, Inf)
    end
end

@testset "check 5: the lift, $T" for T in (Float64, Float32)
    n, k, nsteps = 1000, 20, 1000
    dt = T(1e-4)
    q0, p0 = analytic_initial_condition(T; n = n)
    # **The basis is check 3's**: the cotangent-lift SVD of §4's training snapshot at `n = 1000`,
    # saved every `1.0`, at `k = 20`, built in `T` from `T` snapshots
    b = ReducedBasis(CotangentLiftSVD(), training_snapshots(T, q0, p0), NamedTuple(),
        training_parameter_space(T); k = k)
    h = AnalyticHamiltonian(T(ANALYTIC_β))

    @test size(b.Φ) == (n, k)
    # `Φ` is orthonormal to round-off, and what is left of that is the first term of the round-trip
    # bound below: the decomposition's own non-orthogonality, `‖ΦᵀΦ - I‖₂`, measured `16.4 eps` in
    # `Float64` and `14.1 eps` in `Float32`, both inside `sqrt(eps(T))`
    nonorth = opnorm(b.Φ' * b.Φ - Matrix{T}(I, k, k))
    @test nonorth ≤ sqrt(eps(T))

    # the two states of the check: §4's `x₀`, which lies in the training span, and a seeded state
    # that does not. The projection of a stacked state is the lift of its projection, which is what
    # the two `Φ Φᵀ` blocks of `lift` and `project` are.
    offspan = randn(Xoshiro(2), T, 2n)
    @test norm(offspan - lift(b, project(b, offspan))) ≥ T(0.5) * norm(offspan)

    for (what, X0) in (("§4's x₀", [q0; p0]), ("the off-span x₀", offspan))
        q, p = X0[1:n], X0[(n + 1):(2n)]
        prob = reduced_pode(b, h, (zero(T), chunk_span(T, nsteps, dt)), dt, project(b, X0))
        sol = integrate(prob, LobattoIIIAIIIB(2); f_abstol = lobatto_abstol(T, k))

        # the plain chunk span takes the chunk's own `1000` steps, which is what the solution counts
        # and what `ntime(prob)` does not report
        @test ntime(sol) == nsteps

        ql, pl = lift_solution(b, sol)
        # one column per stored time, which is every step of the 1000-step problem
        @test size(ql) == (n, length(sol.q))
        @test size(pl) == (n, length(sol.q))

        # the lifted initial state is `(ΦΦᵀq₀, ΦΦᵀp₀)`, within the worst case of two products with an
        # orthonormal `Φ`. On the off-span `x₀` this refuses a lifted state taken as `x₀` itself,
        # which §4's `x₀` — it lies in the span — cannot see.
        tol = T(n + k) * sqrt(T(k)) * eps(T) * norm(X0)
        @test norm(ql[:, 1] - b.Φ * (b.Φ' * q)) ≤ tol
        @test norm(pl[:, 1] - b.Φ * (b.Φ' * p)) ≤ tol
    end

    # the round trip of the reduced coordinates, at 100 seeded `z`, against the basis's own
    # non-orthogonality plus the rounding of the two products,
    # `‖project(lift(z)) - z‖ ≤ (‖ΦᵀΦ - I‖₂ + 10·eps(T))·‖z‖`: measured `10.67 ≤ 26.4 eps` in
    # `Float64` and `7.33 ≤ 24.1 eps` in `Float32`, and the round trip sits below `‖ΦᵀΦ - I‖₂` alone
    Z = randn(Xoshiro(3), T, 2k, 100)
    err = maximum(j -> norm(project(b, lift(b, Z[:, j])) - Z[:, j]) / norm(Z[:, j]),
        axes(Z, 2))
    @test err ≤ nonorth + 10 * eps(T)
end

@testset "the gradient is inferred and allocates nothing, $T" for T in (Float64, Float32)
    h = AnalyticHamiltonian(T(ANALYTIC_β))
    x = randn(Xoshiro(5), T, 1000)
    g = zeros(T, 1000)

    @test @inferred(h(g, zero(T), x, NamedTuple())) === g

    # the value half of the instance, `h(t, x, params)`, and not only the gradient: on a small
    # integer state both blocks sum exactly, so the assertion is `==` against the closed form
    # `½(β²|q|² + |p|²)` of the model, `½(36·30 + 174) = 627`, and not a tolerance
    xh = T[1, 2, 3, 4, 5, 6, 7, 8]
    @test h(zero(T), xh, NamedTuple()) ==
          (T(ANALYTIC_β)^2 * (1 + 4 + 9 + 16) + (25 + 36 + 49 + 64)) / 2

    if Base.JLOptions().code_coverage == 0
        @test _alloc4(h, g, zero(T), x, NamedTuple()) == 0
    end
    @test _alloc4((g, t, x, params) -> similar(g), g, zero(T), x, NamedTuple()) > 0

    # the `Dual` states of Newton's `JacobianAutodiff`, which the Lobatto run above passes
    xd = ForwardDiff.Dual.(x, one(T))
    gd = similar(xd)
    @test @inferred(h(gd, zero(T), xd, NamedTuple())) === gd

    if Base.JLOptions().code_coverage == 0
        @test _alloc4(h, gd, zero(T), xd, NamedTuple()) == 0
    end
    @test _alloc4((g, t, x, params) -> similar(g), gd, zero(T), xd, NamedTuple()) > 0
end
