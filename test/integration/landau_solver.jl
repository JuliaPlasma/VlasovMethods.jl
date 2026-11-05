using ParticleMethods: ParticleList
using SimpleSplines: TensorProductQuadrature
using Test
using VlasovMethods

const VM = VlasovMethods

# No other test calls `Picard_iterate_Landau_nls!`. It is called once on each side of its
# `ti ≥ 4` branch, so an undefined global in either branch throws here.
function landau_solver_setup(::Type{T} = Float64) where {T}
    n1 = 16
    npart = n1^2
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pdist = ParticleDistribution(1, 2, ParticleList(zeros(T, 4, npart); variables = vars))
    # Cell centres inside the velocity domain, so the projected density stays positive.
    g = [T(-2) + (T(i) - T(1) / 2) * T(4) / n1 for i in 1:n1]
    pdist.particles.v[1, :] .= repeat(g; outer = n1)
    pdist.particles.v[2, :] .= repeat(g; inner = n1)
    pdist.particles.w[1, :] .= one(T) / npart

    axis = BSplineBasis(UniformMesh(8, T(-2) .. T(2)), 3, Free())
    B = TensorProductBasis((axis, axis))
    sdist = SplineDistribution(1, 2, B, TensorProductQuadrature(B), zeros(T, size(B)...))
    return Landau(pdist, CollisionEntropy(sdist)), copy(pdist.particles.v)
end

# One step at `ti`, with `stdout` redirected into a file; returns what the step returns and what
# it wrote to the file.
function landau_step(landau, v_prev::AbstractMatrix{T}, ti) where {T}
    rhs_prev = zeros(T, 2, size(v_prev, 2), 2)
    path = tempname()
    sol = open(path, "w") do io
        redirect_stdout(io) do
            VM.Picard_iterate_Landau_nls!(
                landau, T(1e-6), ti, zero(T), v_prev, copy(v_prev),
                rhs_prev; abstol = 1e-12)
        end
    end
    return sol, read(path, String)
end

# The kinetic energy `½ Σ_α w_α ‖v_α‖₂²` and the momentum `Σ_α w_α v_α`.
landau_energy(v, w) = sum(w .* vec(sum(abs2, v; dims = 1))) / 2
landau_momentum(v, w) = v * w

@testset "the Landau Picard solver runs" begin
    landau1, v_prev1 = landau_solver_setup()
    sol1, out1 = landau_step(landau1, v_prev1, 1)
    landau4, v_prev4 = landau_solver_setup()
    sol4, out4 = landau_step(landau4, v_prev4, 4)

    for sol in (sol1, sol4)
        @test sol isa NamedTuple{(:iterations, :residual)}
        @test sol.iterations ≥ 1
        @test sol.residual ≤ 1e-12
    end

    # the step moves the particles, from the midpoint guess of `ti < 4` and from the Hermite
    # guess of `ti ≥ 4`
    @test landau1.dist.particles.v != v_prev1
    @test landau4.dist.particles.v != v_prev4

    # nothing is written to `stdout`
    @test isempty(out1)
    @test isempty(out4)
end

@testset "the step converges and conserves" begin
    Δt, abstol = 1e-3, 1e-12
    landau, v_prev = landau_solver_setup()
    w = vec(landau.dist.particles.w)
    rhs_prev = zeros(2, size(v_prev, 2), 2)
    sol = VM.Picard_iterate_Landau_nls!(landau, Δt, 4, 0.0, v_prev, copy(v_prev), rhs_prev;
        abstol)
    v_new = landau.dist.particles.v
    f = rhs_prev[:, :, 1]

    @test sol.residual ≤ abstol
    @test sol.iterations ≤ 20
    @test v_new != v_prev
    # the accepted state is the update of the stored derivative
    @test v_new == v_prev .+ Δt .* f

    # momentum: every Picard update `v_prev + Δt f` keeps `Σ w v` to round-off
    vnorms = sum(w .* vec(sqrt.(sum(abs2, v_prev; dims = 1))))
    @test all(abs.(landau_momentum(v_new, w) .- landau_momentum(v_prev, w)) .≤
              64 * eps() * vnorms)

    # energy: with `r_k = v_k − v_new` and `v_new − v_prev = Δt f_k`, where `f_k` is the vector
    # field at the midpoint `m_k = (v_k + v_prev)/2`,
    #   E(v_new) − E(v_prev) = Σ w (v_new − v_prev)·(v_new + v_prev)/2
    #                        = Δt Σ w f_k·(m_k − r_k/2) = −(Δt/2) Σ w f_k·r_k,
    # because `Σ w f·m = 0` for the vector field at its own midpoint. So the drift is at most
    # `(Δt/2) · abstol · Σ w ‖f_k‖₁`, plus the round-off of the two energies.
    E = landau_energy(v_prev, w)
    bound = Δt / 2 * abstol * sum(w .* vec(sum(abs, f; dims = 1))) + 64 * eps() * E
    @test abs(landau_energy(v_new, w) - E) ≤ bound

    # the identity behind the bound: `Σ w f(v)·v = 0` for the Landau vector field at `v`
    params = (dist = landau.dist, ent = landau.entropy)
    fv = similar(v_prev)
    VM.collisional_vectorfield!(fv, v_prev, params, landau)
    @test abs(sum(w .* vec(sum(fv .* v_prev; dims = 1)))) ≤
          64 * eps() * sum(w .* vec(sum(abs.(fv .* v_prev); dims = 1)))

    # the control: five Picard iterations and the last iterate kept, as a fixed loop takes them.
    # With `v_prev_2 == v_prev` and `rhs_prev == 0`, the Hermite guess is `v_prev`.
    v, g = copy(v_prev), similar(v_prev)
    VM.collisional_vectorfield!(g, (v .+ v_prev) ./ 2, params, landau)
    for _ in 1:5
        v .= v_prev .+ Δt .* g
        VM.collisional_vectorfield!(g, (v .+ v_prev) ./ 2, params, landau)
    end
    @test abs(landau_energy(v, w) - E) > bound
end

@testset "the residuals are debug messages" begin
    landau, v_prev = landau_solver_setup()
    Δt, abstol = 1e-3, 1e-12
    logger = Test.TestLogger(; min_level = Base.CoreLogging.Debug)
    sol = Base.CoreLogging.with_logger(logger) do
        VM.Picard_iterate_Landau_nls!(landau, Δt, 4, 0.0, v_prev, copy(v_prev),
            zeros(2, size(v_prev, 2), 2); abstol)
    end
    records = filter(r -> r.message == "Landau Picard residual", logger.logs)
    # one message per vector-field evaluation, with its iterate's max-norm residual
    @test length(records) == sol.iterations
    @test all(r -> r.level == Base.CoreLogging.Debug, records)
    @test [r.kwargs[:iteration] for r in records] == 0:(sol.iterations - 1)
    @test all(r -> r.kwargs[:residual] > abstol, records[1:(end - 1)])
    @test records[end].kwargs[:residual] == sol.residual
    # With `v_prev_2 == v_prev` and `rhs_prev == 0`, the guess is `v_prev`, so its residual is
    # the max-norm of `Δt f(v_prev)`, up to the rounding of `v_prev + Δt f` (relative
    # `eps · ‖v‖∞ / (Δt ‖f‖∞)`, about 1e-13 here).
    landau0, _ = landau_solver_setup()
    f0 = similar(v_prev)
    VM.collisional_vectorfield!(f0, v_prev, (dist = landau0.dist, ent = landau0.entropy), landau0)
    @test isapprox(records[1].kwargs[:residual], Δt * maximum(abs, f0); rtol = 1e-10)
end

@testset "the step shifts the stored vector fields" begin
    landau, v_prev = landau_solver_setup()
    Δt = 1e-6
    rhs_prev = zeros(2, size(v_prev, 2), 2)
    rhs_prev[:, :, 1] .= 1
    rhs_prev[:, :, 2] .= 2
    VM.Picard_iterate_Landau_nls!(landau, Δt, 4, 0.0, v_prev, copy(v_prev), rhs_prev;
        abstol = 1e-12)
    # the accepted vector field enters column 1, and column 1 moves to column 2
    @test landau.dist.particles.v == v_prev .+ Δt .* rhs_prev[:, :, 1]
    @test all(==(1), rhs_prev[:, :, 2])
end

@testset "no convergence throws" begin
    # At Δt = 1e-2 the residuals of this fixture do not fall, and the 13th evaluation meets
    # particles that left [−2, 2]², so `max_iters` stays below 13.
    landau, v_prev = landau_solver_setup()
    rhs_prev = zeros(2, size(v_prev, 2), 2)
    logger = Test.TestLogger(; min_level = Base.CoreLogging.Debug)
    err = Base.CoreLogging.with_logger(logger) do
        try
            VM.Picard_iterate_Landau_nls!(landau, 1e-2, 4, 0.0, v_prev, copy(v_prev),
                rhs_prev; abstol = 1e-12, max_iters = 10)
            nothing
        catch e
            e
        end
    end
    # `max_iters` vector-field evaluations, each with its residual
    @test count(r -> r.message == "Landau Picard residual", logger.logs) == 10
    @test err isa ErrorException
    @test occursin("did not converge", err.msg)
    @test occursin("10", err.msg)
    # the message names the last residual
    @test occursin(string(logger.logs[end].kwargs[:residual]), err.msg)
    @test landau.dist.particles.v == v_prev
    @test all(iszero, rhs_prev)
end

@testset "refused arguments" begin
    for ti in (1, 4),
        kw in ((abstol = 0.0,), (abstol = -1e-12,), (abstol = NaN,), (abstol = Inf,),
            (abstol = 1e-12, max_iters = 0))

        landau, v_prev = landau_solver_setup()
        rhs_prev = zeros(2, size(v_prev, 2), 2)
        # every vector-field evaluation, the guess's of `ti < 4` too, writes `K1` of the cache
        K1 = copy(landau.cache[Float64].K1)
        logger = Test.TestLogger(; min_level = Base.CoreLogging.Debug)
        Base.CoreLogging.with_logger(logger) do
            @test_throws ArgumentError VM.Picard_iterate_Landau_nls!(
                landau, 1e-3, ti, 0.0, v_prev, copy(v_prev), rhs_prev; kw...)
        end
        @test landau.dist.particles.v == v_prev
        @test all(iszero, rhs_prev)
        # no vector field was evaluated
        @test !any(r -> r.message == "Landau Picard residual", logger.logs)
        @test isequal(landau.cache[Float64].K1, K1)
    end
end

@testset "a residual equal to abstol is accepted" begin
    landau, v_prev = landau_solver_setup()
    sol = VM.Picard_iterate_Landau_nls!(landau, 1e-3, 4, 0.0, v_prev, copy(v_prev),
        zeros(2, size(v_prev, 2), 2); abstol = 1e-12)
    # the same step with `abstol` at the residual it accepted stops at the same evaluation
    landau2, _ = landau_solver_setup()
    sol2 = VM.Picard_iterate_Landau_nls!(landau2, 1e-3, 4, 0.0, v_prev, copy(v_prev),
        zeros(2, size(v_prev, 2), 2); abstol = sol.residual)
    @test sol2 == sol
    @test landau2.dist.particles.v == landau.dist.particles.v
end

@testset "the residual is the max-norm" begin
    # `v - (v_prev + Δt v̇)` is `[1 -3]` here: the max-norm is 3, the largest entry 1
    @test VM._landau_picard_residual([0.0 0.0], [0.0 0.0], 1.0, [-1.0 3.0]) == 3.0
    # a `NaN` entry anywhere makes the residual `NaN`, which fails every `≤ abstol`
    @test isnan(VM._landau_picard_residual([0.0 0.0 0.0], [0.0 0.0 0.0], 1.0, [0.0 NaN 1.0
                                                                               ]))
    @test isnan(VM._landau_picard_residual([0.0 0.0], [0.0 0.0], 1.0, [1.0 NaN]))
end

# The bytes of one call after a warm call, through barriers of fixed arity whose arguments have
# concrete types: a step of the solver, the guess of its first three steps, and a least-squares
# solve of the vector field.
function step_bytes(f::F, l::L, Δt::DT, ti::Int, t::TT, a::A, b::B, r::R,
        abstol::Float64) where {F, L, DT, TT, A, B, R}
    f(l, Δt, ti, t, a, b, r; abstol)
    local sol
    bytes = @allocated sol = f(l, Δt, ti, t, a, b, r; abstol)
    return bytes, sol.iterations
end
function bytes6(f::F, a::A, b::B, c::C, d::D, e::E, g::G) where {F, A, B, C, D, E, G}
    f(a, b, c, d, e, g)
    return @allocated f(a, b, c, d, e, g)
end
function bytes4(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    f(a, b, c, d)
    return @allocated f(a, b, c, d)
end
solve(K, LJ, _, _) = K \ LJ
copy_first(a, b, c, d, e, g) = copy(a)

@testset "the Landau Picard solver allocates within its bounds, $T" for T in (Float32, Float64)
    landau, v_prev = landau_solver_setup(T)
    Δt, t = T(1e-6), zero(T)
    abstol = Float64(64 * eps(T))
    rhs_prev = zeros(T, 2, size(v_prev, 2), 2)
    cache = landau.cache[T]
    params = (dist = landau.dist, ent = landau.entropy)

    # the bytes of the measured call and its number of vector-field evaluations
    picard, iterations = step_bytes(
        VM.Picard_iterate_Landau_nls!, landau, Δt, 4, t, v_prev,
        copy(v_prev), rhs_prev, abstol)
    guess = bytes6(
        VM._landau_startup_guess!, similar(v_prev), v_prev, t, Δt, params, landau)
    # the ceiling of one vector-field evaluation: its two solves, on the matrices that the
    # last evaluation leaves in the cache
    vf = bytes4(solve, cache.K1, cache.LJ, nothing, nothing) +
         bytes4(solve, cache.K2, cache.LJ, nothing, nothing)

    # coverage instrumentation allocates on every counted line, so under it the assertions do
    # not run rather than being weakened
    if Base.JLOptions().code_coverage == 0
        @test picard ≤ iterations * vf + 1024
        @test guess ≤ 37 * vf + 200_000
    end
    # the controls: each barrier sees an allocation
    @test vf > 0
    @test bytes6(copy_first, v_prev, v_prev, t, Δt, params, landau) > 0
    @test picard > 0
end
