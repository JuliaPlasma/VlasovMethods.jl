using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, TensorProductBasis,
                     TensorProductQuadrature
using StableRNGs
using StaticArrays: SVector
using Test
using VlasovMethods

include("../helpers/allocation_barriers.jl")
include("../helpers/pairwise_landau.jl")

copy_first(a, b, c, d) = copy(a)

@testset "the model holds a Landau" begin
    v, w = fixture_particles(Float64)
    model = build(v, w; γ = 0, ν = 0.7)
    @test model.landau isa Landau
    @test model.landau.γ == 0
    @test model.landau.ν == 0.7
    @test model.landau.dist.particles.v == v
    @test PairwiseLandau(model.landau.dist, model.landau.entropy).landau.dist ===
          model.landau.dist
    @test_throws ArgumentError build(v, w; γ = -6)
end

@testset "the vector field is eq:eom_line3, γ = $γ" for γ in (-3, 0)
    model = fixture(; γ, ν = 0.7)
    v = copy(model.landau.dist.particles.v)
    w = weights(model)
    v̇ = similar(v)
    @test VM.collisional_vectorfield!(v̇, v, nothing, model) === nothing

    ref = eom_line3(v, model)
    @test maximum(abs, v̇ .- ref) ≤ 1e-13 * maximum(abs, v̇)

    # momentum and energy, against the scales of their sums
    Pscale = sum(w .* vec(sum(abs, v̇; dims = 1)))
    @test all(abs.(momentum(v̇, w)) .≤ 1e-14 * Pscale)
    vv̇ = vec(sum(v .* v̇; dims = 1))
    @test abs(sum(w .* vv̇)) ≤ 1e-14 * sum(w .* abs.(vv̇))

    # S_h falls, and the anisotropy T_x/T_y = 1.2146 relaxes
    @test sum(gradS(v, model) .* v̇) < 0
    u = momentum(v, w) ./ sum(w)
    @test sum(w .* (v[1, :] .- u[1]) .* v̇[1, :]) < 0 <
          sum(w .* (v[2, :] .- u[2]) .* v̇[2, :])

    # the control: the transposed weight ratio w_γ/w_β breaks the momentum
    v̇c = 0.7 .* G_written(v, gradS(v, model), w, kernelof(model); transposed = true)
    Pc = sum(w .* vec(sum(abs, v̇c; dims = 1)))
    @test any(abs.(momentum(v̇c, w)) .≥ 1e-2 * Pc)
end

@testset "the kernel of 𝔾, γ = $γ" for γ in (-3, 0)
    rng = StableRNG(11)
    v = randn(rng, 2, 9)
    w = 0.2 .+ rand(rng, 9)
    model = build(v, w; γ)

    G = zeros(18, 18)
    y = zeros(2, 9)
    for j in 1:18
        e = zeros(2, 9)
        e[j] = 1
        VM.apply_G!(y, v, e, model)
        G[:, j] .= vec(y)
    end
    @test maximum(abs, G .- G') ≤ 1e-14 * maximum(abs, G)
    λ = eigvals(Symmetric(G))
    @test maximum(λ) ≤ 1e-12 * maximum(abs, λ)
    s = svdvals(G)
    @test count(s .< 1e-10 * s[1]) == 3
    @test sort(s)[4] > 1e-3 * s[1]

    # the control: U = I/|u| keeps the momentum directions, not the energy direction
    Gc = zeros(18, 18)
    for j in 1:18
        e = zeros(2, 9)
        e[j] = 1
        Gc[:, j] .= vec(G_written(v, e, w, (a, b) -> I(2) ./ norm(a .- b)))
    end
    sc = svdvals(Gc)
    @test count(sc .< 1e-10 * sc[1]) == 2
end

@testset "equal velocities throw" begin
    model = fixture()
    v = model.landau.dist.particles.v
    v[:, 23] .= v[:, 7]
    before = copy(v)
    msg = "particles 7 and 23 have equal velocities"
    x = ones(2, 144)
    @test throws_with(() -> VM.apply_G!(similar(before), before, x, model), DomainError, msg)
    @test throws_with(
        () -> VM.collisional_vectorfield!(similar(before), before, nothing, model),
        DomainError, msg)
    @test model.landau.dist.particles.v == before

    # one equal component is no equal velocity
    v = copy(fixture().landau.dist.particles.v)
    v[1, 42] = v[1, 30]
    y = similar(v)
    VM.apply_G!(y, v, x, build(v, weights(model)))
    @test all(isfinite, y)
end

# Two different velocities at the distance d: the Coulomb kernel is -Inf and NaN where 1/|u|³
# overflows (d = 1e-13 in Float32, 1e-110 in Float64), and the zero matrix where |u|² underflows
# (1e-30, 1e-200).
@testset "near-equal velocities throw, $T, d = $d" for (T, d) in (
    (Float32, 1.0f-13), (Float32, 1.0f-30), (Float64, 1e-110), (Float64, 1e-200))
    v, w = fixture_particles(T)
    v[:, 7] .= 0
    v[:, 23] .= (d, 0)
    model = build(v, w)
    before = copy(v)
    msg = "particles 7 and 23 have velocities too close"
    x = ones(T, 2, 144)
    @test throws_with(() -> VM.apply_G!(similar(before), before, x, model), DomainError, msg)
    @test throws_with(
        () -> VM.collisional_vectorfield!(similar(before), before, nothing, model),
        DomainError, msg)
    @test model.landau.dist.particles.v == before
end

# The threshold τ of |v_γ - v_β|²: at the smallest d with d² ≥ τ the kernel is finite and not
# zero and apply_G! returns finite values; at the largest d with d² < τ apply_G! throws; at
# d² ≈ τ/16 the kernel itself is no longer finite, so τ is within a factor 16 of the failure.
@testset "the threshold of near-equal velocities, $T, γ = $γ" for T in (Float32, Float64),
    γ in (-3, -4.5)

    v, w = fixture_particles(T)
    model = build(v, w; γ)
    τ = VM._close_threshold(T, model.landau)
    @test floatmin(T) ≤ τ < 1
    U(d) = VM.kernel(SVector(zero(T), zero(T)), SVector(d, zero(T)), model.landau)
    x = ones(T, 2, 144)

    above = sqrt(τ)
    while above * above < τ
        above = nextfloat(above)
    end
    below = prevfloat(above)
    @test below * below < τ ≤ above * above
    @test all(isfinite, U(above)) && !iszero(U(above))
    @test !all(isfinite, U(sqrt(τ) / 4))

    v[:, 7] .= 0
    v[:, 23] .= (above, 0)
    y = similar(v)
    VM.apply_G!(y, v, x, build(v, w; γ))
    @test all(isfinite, y)
    v[:, 23] .= (below, 0)
    @test throws_with(() -> VM.apply_G!(y, v, x, build(v, w; γ)), DomainError,
        "particles 7 and 23 have velocities too close")
end

# For γ ≥ -2 the kernel tends to a finite value or to zero as |u| → 0, and is not refused.
@testset "γ = 0 accepts near-equal velocities" begin
    v, w = fixture_particles(Float64)
    v[:, 7] .= 0
    v[:, 23] .= (1e-200, 0)
    model = build(v, w; γ = 0)
    @test VM._close_threshold(Float64, model.landau) == 0
    y = similar(v)
    VM.apply_G!(y, v, ones(2, 144), model)
    @test all(isfinite, y)
    VM.collisional_vectorfield!(y, v, nothing, model)
    @test all(isfinite, y)
end

@testset "one particle, no particle, NaN" begin
    one_particle = build(reshape([0.1, -0.2], 2, 1), [1.0])
    v̇ = ones(2, 1)
    VM.collisional_vectorfield!(v̇, [0.1 -0.2]', nothing, one_particle)
    @test iszero(v̇)

    none = build(zeros(2, 0), zeros(0))
    @test_throws DomainError VM.collisional_vectorfield!(zeros(2, 0), zeros(2, 0), nothing, none)

    model = fixture()
    model.landau.dist.particles.v[1, 5] = NaN
    before = copy(model.landau.dist.particles.v)
    @test_throws DomainError VM.collisional_vectorfield!(similar(before), before, nothing, model)
    @test isequal(model.landau.dist.particles.v, before)
end

@testset "particle weights and sizes" begin
    # the scheme divides by w_γ: a weight that is not positive and finite is refused
    v, w = fixture_particles(Float64)
    for bad in (0.0, -1e-3, NaN, Inf)
        wb = copy(w)
        wb[17] = bad
        @test throws_with(() -> build(v, wb), ArgumentError, "weight of particle 17")
    end

    # the model holds the distribution by reference, so apply_G! and the vector field check
    # the weights as they are at the call, the vector field before it projects the particles
    for bad in (0.0, -1e-3, NaN, Inf)
        model = fixture()
        model.landau.dist.particles.w[1, 17] = bad
        vm = copy(model.landau.dist.particles.v)
        @test throws_with(() -> VM.apply_G!(similar(vm), vm, ones(2, 144), model),
            ArgumentError, "weight of particle 17")
        @test throws_with(
            () -> VM.collisional_vectorfield!(similar(vm), vm, nothing, model),
            ArgumentError, "weight of particle 17")
    end

    # an array of another number of columns than the model has particles
    model = fixture()
    v = copy(model.landau.dist.particles.v)
    v100 = v[:, 1:100]
    msg = "the model has 144 particles"
    @test throws_with(() -> VM.entropy_gradient!(similar(v100), v100, model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.entropy_gradient!(zeros(2, 100), v, model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.entropy_gradient!(similar(v), v100, model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.apply_G!(similar(v), v100, ones(2, 144), model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.apply_G!(similar(v100), v100, ones(2, 100), model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.apply_G!(similar(v), v, ones(2, 100), model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.apply_G!(zeros(2, 100), v, ones(2, 144), model),
        DimensionMismatch, msg)
    @test throws_with(
        () -> VM.collisional_vectorfield!(similar(v100), v100, nothing, model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.collisional_vectorfield!(zeros(2, 100), v, nothing, model),
        DimensionMismatch, msg)
    v145 = hcat(v, [0.05, 0.05])
    @test throws_with(() -> VM.entropy_gradient!(similar(v145), v145, model),
        DimensionMismatch, msg)
    @test throws_with(() -> VM.apply_G!(similar(v145), v145, ones(2, 145), model),
        DimensionMismatch, msg)

    # the vector field names its own argument
    @test throws_with(() -> VM.collisional_vectorfield!(zeros(2, 145), v, nothing, model),
        DimensionMismatch, "v̇ has size (2, 145)")

    # an array of another number of rows than the two velocity components
    v3 = vcat(v, zeros(1, 144))
    @test throws_with(() -> VM.entropy_gradient!(similar(v), v3, model),
        DimensionMismatch, "v has size (3, 144)")
    @test throws_with(() -> VM.entropy_gradient!(zeros(3, 144), v, model),
        DimensionMismatch, "∇S has size (3, 144)")
    @test throws_with(() -> VM.apply_G!(similar(v), v3, ones(2, 144), model),
        DimensionMismatch, "v has size (3, 144)")
    @test throws_with(() -> VM.apply_G!(similar(v), v, ones(3, 144), model),
        DimensionMismatch, "x has size (3, 144)")
    @test throws_with(() -> VM.apply_G!(zeros(3, 144), v, ones(2, 144), model),
        DimensionMismatch, "y has size (3, 144)")
    @test throws_with(() -> VM.collisional_vectorfield!(zeros(3, 144), v, nothing, model),
        DimensionMismatch, "v̇ has size (3, 144)")
    @test throws_with(() -> VM.collisional_vectorfield!(zeros(2, 144), v3, nothing, model),
        DimensionMismatch, "v has size (3, 144)")
    @test model.landau.dist.particles.v == v
end

@testset "Float32, allocations and inference" begin
    v32, w32 = fixture_particles(Float32)
    m32 = build(v32, w32)
    m64 = build(Float64.(v32), Float64.(w32))
    @test m32.landau.ν isa Float64

    f32 = similar(v32)
    VM.collisional_vectorfield!(f32, v32, nothing, m32)
    f64 = zeros(2, 144)
    VM.collisional_vectorfield!(f64, Float64.(v32), nothing, m64)
    @test eltype(f32) == Float32
    @test maximum(abs, f32 .- f64) ≤ 1e-4 * maximum(abs, f64)
    Pscale = sum(w32 .* vec(sum(abs, f32; dims = 1)))
    @test all(abs.(momentum(f32, w32)) .≤ 64 * eps(Float32) * Pscale)
    vf = vec(sum(v32 .* f32; dims = 1))
    @test abs(sum(w32 .* vf)) ≤ 64 * eps(Float32) * sum(w32 .* abs.(vf))

    for (model, v) in ((m32, v32), (m64, Float64.(v32)))
        x = ones(eltype(v), size(v))
        y = similar(v)
        @test _alloc4(VM.apply_G!, y, v, x, model) == 0
        @test _alloc4(VM.collisional_vectorfield!, y, v, nothing, model) == 0
        @test _alloc4(copy_first, y, v, x, model) > 0
        @inferred VM.apply_G!(y, v, x, model)
        @inferred VM.collisional_vectorfield!(y, v, nothing, model)
    end
end
