using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, TensorProductBasis,
                     TensorProductQuadrature
using StableRNGs
using StaticArrays: SVector
using Test
using VlasovMethods

include("../helpers/pairwise_landau.jl")

@testset "the Gonzalez gradient" begin
    model = fixture()
    v0 = copy(model.landau.dist.particles.v)
    g0 = similar(v0)
    VM.entropy_gradient!(g0, v0, model)
    @test g0 ≈ gradS(v0, model) rtol=1e-14

    # Δv = 0: the gradient at v0, bit for bit
    gd = similar(v0)
    VM.discrete_gradient!(gd, v0, copy(v0), model, GonzalezGradient())
    @test gd == g0
    @test all(isfinite, gd)

    # Δv ≠ 0: the discrete gradient identity, and the midpoint gradient as the control
    d = randn(StableRNG(5), 2, 144)
    v1 = v0 .+ 1e-3 .* d ./ norm(d)
    VM.discrete_gradient!(gd, v0, v1, model, GonzalezGradient())
    ΔS = Sh(v1, model) - Sh(v0, model)
    @test ΔS < 0
    @test abs(ΔS - dot(v1 .- v0, gd)) ≤ 1e-12 * abs(ΔS)
    gm = similar(v0)
    VM.entropy_gradient!(gm, (v0 .+ v1) ./ 2, model)
    @test abs(ΔS - dot(v1 .- v0, gm)) > 1e-10 * abs(ΔS)

    # the base point is the midpoint: the formula of eq:disc_grad_props, written out
    Δv = v1 .- v0
    g = gradS((v0 .+ v1) ./ 2, model)
    gref = g .+ Δv .* ((Sh(v1, model) - Sh(v0, model) - dot(Δv, g)) / sum(abs2, Δv))
    @test gd ≈ gref rtol=1e-12

    # only a computed |Δv|² of zero takes the midpoint gradient: at |Δv|² = 1e-10 the
    # correction applies
    v1 = v0 .+ 1e-5 .* d ./ norm(d)
    VM.discrete_gradient!(gd, v0, v1, model, GonzalezGradient())
    VM.entropy_gradient!(gm, (v0 .+ v1) ./ 2, model)
    @test gd != gm
    ΔS = Sh(v1, model) - Sh(v0, model)
    @test abs(ΔS - dot(v1 .- v0, gd)) ≤ 1e-12 * abs(ΔS)
end

@testset "one step conserves" begin
    model = fixture()
    vn = copy(model.landau.dist.particles.v)
    w = weights(model)
    Δt, abstol = 2e-3, 1e-12

    r = discrete_gradient_step!(model, Δt; abstol)
    v1 = copy(model.landau.dist.particles.v)
    @test r.residual ≤ abstol
    @test 1 ≤ r.iterations ≤ 10
    @test v1 != vn
    @test all(abs.(momentum(v1, w) .- momentum(vn, w)) .≤ momentum_bound(vn, w))
    @test abs(energy(v1, w) - energy(vn, w)) ≤ energy_bound(vn, v1, Δt, abstol, w)

    # the control: 𝔾(vⁿ) in place of 𝔾(m_k) breaks the energy, at order Δt²
    ΔE = map((2e-3, 1e-3)) do h
        vc = written_step(model, vn, h, abstol; G_at_old = true)
        e = abs(energy(vc, w) - energy(vn, w))
        @test e > energy_bound(vn, vc, h, abstol, w)
        e
    end
    @test 3.5 ≤ ΔE[1] / ΔE[2] ≤ 4.5

    # ν = 0.7 scales f: the step is the written-out step
    m7 = fixture(; ν = 0.7)
    discrete_gradient_step!(m7, Δt; abstol)
    @test maximum(abs, m7.landau.dist.particles.v .- written_step(m7, vn, Δt, abstol)) ≤
          abstol
end

@testset "the accepted state is the update" begin
    model = fixture()
    vn = copy(model.landau.dist.particles.v)
    w = weights(model)
    Δt, abstol = 8e-3, 1e-6

    r = discrete_gradient_step!(model, Δt; abstol)
    v1 = copy(model.landau.dist.particles.v)
    @test r.iterations == 3
    @test abs(energy(v1, w) - energy(vn, w)) ≤ energy_bound(vn, v1, Δt, abstol, w)

    # the control: the iterate v_k accepted in place of vⁿ + Δt f_k
    vc = written_step(model, vn, Δt, abstol; accept_iterate = true)
    @test abs(energy(vc, w) - energy(vn, w)) > energy_bound(vn, vc, Δt, abstol, w)
end

@testset "S_h decreases over 20 steps" begin
    model = fixture()
    w = weights(model)
    Δt, abstol = 2e-3, 1e-12
    ν = model.landau.ν
    v0 = copy(model.landau.dist.particles.v)
    gd = similar(v0)
    gm = similar(v0)
    Gg = similar(v0)
    Ebound = 0.0
    for n in 1:20
        vn = copy(model.landau.dist.particles.v)
        discrete_gradient_step!(model, Δt; abstol)
        v1 = copy(model.landau.dist.particles.v)
        Ebound += energy_bound(vn, v1, Δt, abstol, w)

        ΔS = Sh(v1, model) - Sh(vn, model)
        Δv = v1 .- vn
        VM.discrete_gradient!(gd, vn, v1, model, GonzalezGradient())
        VM.apply_G!(Gg, (vn .+ v1) ./ 2, gd, model)
        @test ΔS < 0
        @test abs(ΔS - dot(Δv, gd)) ≤ 1e-12 * abs(ΔS)
        @test abs(ΔS - Δt * ν * dot(gd, Gg)) ≤ 1e-9 * abs(ΔS)

        # the control: the midpoint gradient breaks the discrete gradient identity
        VM.entropy_gradient!(gm, (vn .+ v1) ./ 2, model)
        @test abs(ΔS - dot(Δv, gm)) > 1e-10 * abs(ΔS)
    end
    v20 = model.landau.dist.particles.v
    @test all(abs.(momentum(v20, w) .- momentum(v0, w)) .≤ 20 * momentum_bound(v0, w))
    @test abs(energy(v20, w) - energy(v0, w)) ≤ Ebound
    @test anisotropy(v20, w) < anisotropy(v0, w)
end

@testset "refused arguments" begin
    model = fixture()
    before = copy(model.landau.dist.particles.v)
    for Δt in (0.0, -2e-3, NaN, Inf)
        @test throws_with(() -> discrete_gradient_step!(model, Δt; abstol = 1e-12),
            ArgumentError, "Δt")
    end
    for abstol in (0.0, -1e-12, NaN, Inf)
        @test throws_with(
            () -> discrete_gradient_step!(model, 2e-3; abstol), ArgumentError,
            "abstol")
    end
    @test throws_with(
        () -> discrete_gradient_step!(model, 2e-3; abstol = 1e-12, max_iters = 0),
        ArgumentError, "max_iters")
    @test model.landau.dist.particles.v == before

    # Δt is checked in the element type of the velocities: a Float64 Δt that is positive and
    # finite, but zero or Inf in Float32, is refused, and the particles do not move
    m32 = fixture(Float32)
    before32 = copy(m32.landau.dist.particles.v)
    for Δt in (1e-50, 1e300)
        @test throws_with(() -> discrete_gradient_step!(m32, Δt; abstol = 1.0f-4),
            ArgumentError, "Δt must be finite and positive in Float32, got $Δt")
    end
    @test m32.landau.dist.particles.v == before32
    # the control: the same Δt is accepted by a Float64 model
    r = discrete_gradient_step!(fixture(), 1e-50; abstol = 1e-12)
    @test r.iterations ≥ 1
end

@testset "no convergence throws" begin
    model = fixture()
    before = copy(model.landau.dist.particles.v)
    @test throws_with(
        () -> discrete_gradient_step!(model, 2e-3; abstol = 1e-12, max_iters = 2),
        ErrorException, "did not converge in 2 evaluations")
    @test model.landau.dist.particles.v == before

    # a residual equal to abstol is accepted: the last residual of the thrown step, printed
    # in full, as abstol of the same step
    err = try
        discrete_gradient_step!(model, 2e-3; abstol = 1e-12, max_iters = 2)
    catch e
        e
    end
    r2 = parse(Float64, match(r"the last residual is (\S+),", err.msg).captures[1])
    @test discrete_gradient_step!(model, 2e-3; abstol = r2, max_iters = 2) ==
          (iterations = 2, residual = r2)

    # a residual that is not finite never meets abstol: ν = NaN makes f_0 and r_0 NaN (a
    # second evaluation would project the NaN iterate, which throws a DomainError)
    v, w = fixture_particles(Float64)
    nan_model = build(v, w; ν = NaN)
    @test throws_with(
        () -> discrete_gradient_step!(nan_model, 2e-3; abstol = 1e-12, max_iters = 1),
        ErrorException, "did not converge in 1 evaluations")
    @test nan_model.landau.dist.particles.v == v
end

@testset "the errors of the model pass through the step" begin
    # equal velocities
    model = fixture()
    v = model.landau.dist.particles.v
    v[:, 23] .= v[:, 7]
    before = copy(v)
    @test throws_with(
        () -> discrete_gradient_step!(model, 2e-3; abstol = 1e-12), DomainError,
        "particles 7 and 23 have equal velocities")
    @test model.landau.dist.particles.v == before

    # near-equal velocities, where the Coulomb kernel is not finite or is zero
    for (T, d) in ((Float32, 1.0f-13), (Float32, 1.0f-30), (Float64, 1e-110), (
        Float64, 1e-200))
        v, w = fixture_particles(T)
        v[:, 7] .= 0
        v[:, 23] .= (d, 0)
        model = build(v, w)
        @test throws_with(
            () -> discrete_gradient_step!(model, T(2e-3); abstol = T(1e-4)), DomainError,
            "particles 7 and 23 have velocities too close")
        @test model.landau.dist.particles.v == v
    end

    # a weight that is not positive and finite, set after construction
    for bad in (0.0, -1e-3, NaN, Inf)
        model = fixture()
        before = copy(model.landau.dist.particles.v)
        model.landau.dist.particles.w[1, 17] = bad
        @test throws_with(() -> discrete_gradient_step!(model, 2e-3; abstol = 1e-12),
            ArgumentError, "weight of particle 17")
        @test model.landau.dist.particles.v == before
    end

    # one particle does not move; no particle and a NaN velocity throw
    one_particle = build(reshape([0.1, -0.2], 2, 1), [1.0])
    @test discrete_gradient_step!(one_particle, 2e-3; abstol = 1e-12) ==
          (iterations = 1, residual = 0.0)
    @test one_particle.landau.dist.particles.v == reshape([0.1, -0.2], 2, 1)
    @test_throws DomainError discrete_gradient_step!(build(zeros(2, 0), zeros(0)), 2e-3;
        abstol = 1e-12)
    model = fixture()
    model.landau.dist.particles.v[1, 5] = NaN
    before = copy(model.landau.dist.particles.v)
    @test_throws DomainError discrete_gradient_step!(model, 2e-3; abstol = 1e-12)
    @test isequal(model.landau.dist.particles.v, before)
end

@testset "Float32" begin
    v32, w32 = fixture_particles(Float32)
    m32 = build(v32, w32)
    Δt, abstol = 2.0f-3, 1.0f-4
    r = discrete_gradient_step!(m32, Δt; abstol)
    v1 = copy(m32.landau.dist.particles.v)
    @test r.residual isa Float32
    @test r.residual ≤ abstol
    @test eltype(v1) == Float32
    @test v1 != v32
    @test all(abs.(momentum(v1, w32) .- momentum(v32, w32)) .≤ momentum_bound(v32, w32))
    @test abs(energy(v1, w32) - energy(v32, w32)) ≤ energy_bound(v32, v1, Δt, abstol, w32)
end
