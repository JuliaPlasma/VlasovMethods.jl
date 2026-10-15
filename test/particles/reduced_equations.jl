using ForwardDiff
using LinearAlgebra
using Random
using ReducedBasisMethods: ReducedBasis, PODSVD, CotangentLiftEVD, CotangentLiftSVD,
                           project,
                           lift, reduced_ode, reduced_pode
using Test
using VlasovMethods
# an explicit only-list, so that no name it exports can collide with one of `VlasovMethods`'
using GeometricIntegrators: LobattoIIIAIIIB, equation, integrate, ntime, parameters

include("../helpers/allocation_barriers.jl")
include("../helpers/analytic_vlasov.jl")
include("../helpers/reduced_runs.jl")
include("../helpers/vlasov_ode_reference.jl")

@testset "the t_sep comparator counts a NaN as an exceedance" begin
    # `!(value ≤ bound)`, in the order the helper writes it: `NaN` is neither `≤` nor `>`, so the
    # mutant that tests `value > bound` drops it and returns 3
    @test first_exceedance([0.0, NaN, 1.0], 0.5) == 2
    @test first_exceedance([0.0, 0.25], 0.5) === nothing
    @test first_exceedance([0.0, 0.5, 0.75], 0.5) == 3
end

# The solver tolerance of the Lobatto IIIA–IIIB(2) runs, as §5V.P32 sets it: `1e-14` in `Float64`
# and `max(8, 2k)·eps(T)` in `Float32`, so neither run depends on the default of whichever
# GeometricIntegratorsBase resolves.
lobatto_abstol(::Type{Float64}, k) = 1e-14
lobatto_abstol(::Type{T}, k) where {T} = max(8, 2k) * eps(T)

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

    # **GeometricIntegrators' count is checked**, and it is the contract the chunks of check 3 rest
    # on: `div(tend - tbegin, tstep, RoundUp)` of the chunk's span is exactly `1000` steps, where the
    # run's own clock, the float span `(t_{j-1}, t_j)` of `t_j = j·1000·dt`, counts another number at
    # some `j` and puts every later save at the wrong time.
    @test ntime(prob) == chunk_steps(T, nsteps, dt) == nsteps
    @test any(
        j -> div(T(j + 1) * (T(nsteps) * dt) - T(j) * (T(nsteps) * dt), dt, RoundUp) !=
             nsteps,
        0:20)

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
    n, k, nsteps = 200, 20, 1000
    dt = T(1e-4)
    q0, p0 = analytic_initial_condition(T; n = n)
    b = ReducedBasis(
        CotangentLiftEVD(), training_snapshots(T, q0, p0; step = 10), NamedTuple(),
        training_parameter_space(T); k = k)
    h = AnalyticHamiltonian(T(ANALYTIC_β))

    @test size(b.Φ) == (n, k)
    # `Φ` is orthonormal to round-off, which is what the round trip below costs. The cotangent EVD
    # is the fixture here and not the SVD: the analytic data is numerically rank 2, so the null
    # columns of `Φ` — and with them `ΦᵀΦ - I` — are the decomposition's own round-off, and the SVD's
    # is the larger. Measured in `Float64` at `n = 200` and `step = 10`, whose training snapshot
    # holds `606` columns, `‖ΦᵀΦ - I‖` is `44.4 eps` in Frobenius for the SVD against `11.9 eps` for
    # the EVD, and the round trip below is `12.8 eps` against `2.9 eps`. The SVD fixture cannot meet
    # the bound of §5V.P37 check 5, and the bound is not loosened for it.
    @test norm(b.Φ' * b.Φ - Matrix{T}(I, k, k)) ≤ sqrt(eps(T))

    # the two states of the check: §4's `x₀`, which lies in the training span, and a seeded state
    # that does not. The projection of a stacked state is the lift of its projection, which is what
    # the two `Φ Φᵀ` blocks of `lift` and `project` are.
    offspan = randn(Xoshiro(2), T, 2n)
    @test norm(offspan - lift(b, project(b, offspan))) ≥ T(0.5) * norm(offspan)

    for (what, X0) in (("§4's x₀", [q0; p0]), ("the off-span x₀", offspan))
        q, p = X0[1:n], X0[(n + 1):(2n)]
        prob = reduced_pode(b, h, (zero(T), chunk_span(T, nsteps, dt)), dt, project(b, X0))
        @test ntime(prob) == chunk_steps(T, nsteps, dt)

        sol = integrate(prob, LobattoIIIAIIIB(2); f_abstol = lobatto_abstol(T, k))
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

    # the round trip of the reduced coordinates, at 100 seeded `z`
    Z = randn(Xoshiro(3), T, 2k, 100)
    err = maximum(j -> norm(project(b, lift(b, Z[:, j])) - Z[:, j]) / norm(Z[:, j]),
        axes(Z, 2))
    @test err ≤ 10 * eps(T)
end

@testset "the gradient is inferred and allocates nothing, $T" for T in (Float64, Float32)
    h = AnalyticHamiltonian(T(ANALYTIC_β))
    x = randn(Xoshiro(5), T, 1000)
    g = zeros(T, 1000)

    @test @inferred(h(g, zero(T), x, NamedTuple())) === g

    if Base.JLOptions().code_coverage == 0
        @test _alloc_gradient(h, g, zero(T), x, NamedTuple()) == 0
    end
    @test _alloc_gradient((g, t, x, params) -> similar(g), g, zero(T), x, NamedTuple()) > 0

    # the `Dual` states of Newton's `JacobianAutodiff`, which the Lobatto run above passes
    xd = ForwardDiff.Dual.(x, one(T))
    gd = similar(xd)
    @test @inferred(h(gd, zero(T), xd, NamedTuple())) === gd

    if Base.JLOptions().code_coverage == 0
        @test _alloc_gradient(h, gd, zero(T), xd, NamedTuple()) == 0
    end
    @test _alloc_gradient((g, t, x, params) -> similar(g), gd, zero(T), xd, NamedTuple()) >
          0
end
