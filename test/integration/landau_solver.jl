using LinearAlgebra: norm
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

# One step at `ti`, with `stdout` redirected into a file; returns the histories and what the step
# wrote to the file.
function landau_step(landau, v_prev::AbstractMatrix{T}, ti) where {T}
    rhs_prev = zeros(T, 2, size(v_prev, 2), 2)
    path = tempname()
    sol = open(path, "w") do io
        redirect_stdout(io) do
            VM.Picard_iterate_Landau_nls!(landau, 1e-10, 1e-10, 0.5, T(1e-6), ti, zero(T),
                v_prev, copy(v_prev), rhs_prev, 3, 1)
        end
    end
    return sol, read(path, String)
end

@testset "the Landau Picard solver runs" begin
    landau1, v_prev1 = landau_solver_setup()
    sol1, out1 = landau_step(landau1, v_prev1, 1)
    @test sol1.v isa Vector
    @test length(sol1.v) == 6
    @test sol1.v̇ isa Vector

    landau4, v_prev4 = landau_solver_setup()
    sol4, out4 = landau_step(landau4, v_prev4, 4)
    @test sol4.v isa Vector
    @test length(sol4.v) == 6
    @test sol4.v̇ isa Vector

    # `sol.v[1]` is the extrapolated initial guess. The `ti ≥ 4` call starts from a Hermite guess
    # and the `ti < 4` call from a midpoint one, so the two first entries differ.
    @test sol4.v[1] != sol1.v[1]

    # the step moves the particles to its last iterate
    @test landau4.dist.particles.v != v_prev4
    @test landau4.dist.particles.v == sol4.v[end]

    # nothing is written to `stdout`
    @test isempty(out1)
    @test isempty(out4)
end

@testset "the residuals are debug messages" begin
    landau, v_prev = landau_solver_setup()
    Δt = 1e-6
    logger = Test.TestLogger(; min_level = Base.CoreLogging.Debug)
    sol = Base.CoreLogging.with_logger(logger) do
        VM.Picard_iterate_Landau_nls!(landau, 1e-10, 1e-10, 0.5, Δt, 4, 0.0, v_prev,
            copy(v_prev), zeros(2, size(v_prev, 2), 2), 3, 1)
    end
    records = filter(r -> r.message == "Landau Picard residual", logger.logs)
    @test length(records) == 6
    @test all(r -> r.level == Base.CoreLogging.Debug, records)
    # the guess and each of the five iterates, with the residual of the iterate and the vector
    # field of its midpoint, which the histories hold
    for (i, r) in enumerate(records)
        @test r.kwargs[:iteration] == i - 1
        @test r.kwargs[:residual] == norm(sol.v[i] .- (v_prev .+ Δt .* sol.v̇[i]))
    end
end

@testset "the step shifts the stored vector fields" begin
    landau, v_prev = landau_solver_setup()
    rhs_prev = zeros(2, size(v_prev, 2), 2)
    rhs_prev[:, :, 1] .= 1
    rhs_prev[:, :, 2] .= 2
    sol = VM.Picard_iterate_Landau_nls!(landau, 1e-10, 1e-10, 0.5, 1e-6, 4, 0.0, v_prev,
        copy(v_prev), rhs_prev, 3, 1)
    # the vector field of the last midpoint enters column 1, and column 1 moves to column 2
    @test rhs_prev[:, :, 1] == sol.v̇[end]
    @test all(==(1), rhs_prev[:, :, 2])
end

# The bytes of one call after a warm call, through barriers of fixed arity whose arguments have
# concrete types: a step of the solver, the guess of its first three steps, and a least-squares
# solve of the vector field.
function step_bytes(f::F, l::L, Δt::DT, ti::Int, t::TT, a::A, b::B,
        r::R) where {F, L, DT, TT, A, B, R}
    f(l, 1e-10, 1e-10, 0.5, Δt, ti, t, a, b, r, 3, 1)
    return @allocated f(l, 1e-10, 1e-10, 0.5, Δt, ti, t, a, b, r, 3, 1)
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
    rhs_prev = zeros(T, 2, size(v_prev, 2), 2)
    cache = landau.cache[T]
    params = (dist = landau.dist, ent = landau.entropy)

    picard = step_bytes(VM.Picard_iterate_Landau_nls!, landau, Δt, 4, t, v_prev,
        copy(v_prev), rhs_prev)
    sol = VM.Picard_iterate_Landau_nls!(landau, 1e-10, 1e-10, 0.5, Δt, 4, t, v_prev,
        copy(v_prev), rhs_prev, 3, 1)
    guess = bytes6(
        VM._landau_startup_guess!, similar(v_prev), v_prev, t, Δt, params, landau)
    # the ceiling of one vector-field evaluation: its two solves, on the matrices that the
    # last evaluation leaves in the cache
    vf = bytes4(solve, cache.K1, cache.LJ, nothing, nothing) +
         bytes4(solve, cache.K2, cache.LJ, nothing, nothing)

    # coverage instrumentation allocates on every counted line, so under it the assertions do
    # not run rather than being weakened
    if Base.JLOptions().code_coverage == 0
        @test picard ≤ 6 * vf + Base.summarysize(sol) + 1024
        @test guess ≤ 37 * vf + 200_000
    end
    # the controls: each barrier sees an allocation
    @test vf > 0
    @test bytes6(copy_first, v_prev, v_prev, t, Δt, params, landau) > 0
    @test picard > 0
end
