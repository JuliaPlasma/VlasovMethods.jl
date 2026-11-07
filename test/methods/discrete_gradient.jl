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

# The bytes of one call of a five-argument method after a warm call, through a barrier whose
# arguments have concrete types; `copy_first_5` is the control that the barrier sees an allocation.
function bytes5(f::F, a::A, b::B, c::C, d::D, e::E) where {F, A, B, C, D, E}
    f(a, b, c, d, e)
    return @allocated f(a, b, c, d, e)
end
copy_first_5(a, b, c, d, e) = copy(a)

# The relative increments ‖v_{k+1} − v_k‖ / ‖v_{k+1} − vⁿ‖ of the fixed-point iteration of the
# step, written out with the rule `gradient` and no acceptance: from v_0 = vⁿ,
# f_k = ν 𝔾(m_k) ∇̄S(vⁿ, v_k), m_k = (vⁿ + v_k)/2 and v_{k+1} = vⁿ + Δt f_k.
function increments(model, vn, Δt, gradient; iters = 40)
    T = eltype(vn)
    h = T(Δt)
    ν = T(model.landau.ν)
    vk = copy(vn)
    g = similar(vk)
    m = similar(vk)
    f = similar(vk)
    out = zeros(T, iters)
    for k in 1:iters
        VM.discrete_gradient!(g, vn, vk, model, gradient)
        m .= (vn .+ vk) ./ 2
        VM.apply_G!(f, m, g, model)
        f .*= ν
        vnext = vn .+ h .* f
        out[k] = norm(vnext .- vk, Inf) / norm(vnext .- vn, Inf)
        vk .= vnext
    end
    return out
end

@testset "the quadrature rule" begin
    for nq in (1, 2, 3, 5, 8)
        rule = AVFGradient(nq)
        τ, ω = rule.nodes, rule.weights
        @test rule.nq == nq
        @test length(τ) == nq == length(ω)
        @test eltype(τ) === Float64 === eltype(ω)
        @test all(t -> 0 < t < 1, τ)
        @test issorted(τ)
        @test all(i -> abs(τ[i] + τ[nq + 1 - i] - 1) ≤ 1e-15, 1:nq)
        @test all(>(0), ω)
        @test abs(sum(ω) - 1) ≤ 1e-15
        # exact on τ^k for every k ≤ 2nq − 1
        @test maximum(k -> abs(sum(ω .* τ .^ k) - 1 / (k + 1)), 0:(2nq - 1)) ≤ 1e-14
    end
    @test AVFGradient(1).nodes == [0.5]
    @test AVFGradient(1).weights == [1.0]

    # the default is the five-point rule
    @test AVFGradient().nq == 5
    @test AVFGradient().nodes == AVFGradient(5).nodes

    # the rule sums its nodes with the weights, written out from its own fields
    model = fixture()
    v0 = copy(model.landau.dist.particles.v)
    d = randn(StableRNG(5), 2, 144)
    d ./= norm(d, Inf)
    v1 = v0 .+ 1e-3 .* d
    g = similar(v0)
    gt = similar(v0)
    ref = similar(v0)
    for nq in (1, 2, 5)
        rule = AVFGradient(nq)
        VM.discrete_gradient!(g, v0, v1, model, rule)
        fill!(ref, 0)
        for i in 1:nq
            VM.entropy_gradient!(gt, v0 .+ rule.nodes[i] .* (v1 .- v0), model)
            ref .+= rule.weights[i] .* gt
        end
        @test maximum(abs, g .- ref) ≤ 1e-14 * maximum(abs, ref)
    end
end

@testset "nq = 1 is the midpoint gradient" begin
    model = fixture()
    v0 = copy(model.landau.dist.particles.v)
    d = randn(StableRNG(5), 2, 144)
    d ./= norm(d, Inf)
    v1 = v0 .+ 1e-3 .* d
    g = similar(v0)
    gm = similar(v0)
    VM.discrete_gradient!(g, v0, v1, model, AVFGradient(1))
    VM.entropy_gradient!(gm, (v0 .+ v1) ./ 2, model)
    @test g == gm
    ΔS = Sh(v1, model) - Sh(v0, model)
    @test abs(ΔS - dot(v1 .- v0, g)) == abs(ΔS - dot(v1 .- v0, gm))

    # the single node is the midpoint of the segment, bit for bit
    @test v0 .+ AVFGradient(1).nodes[1] .* (v1 .- v0) == (v0 .+ v1) ./ 2
end

@testset "the identity at a step's displacement" begin
    model = fixture()
    vn = copy(model.landau.dist.particles.v)
    Δt, abstol = 2e-3, 1e-12
    discrete_gradient_step!(model, Δt; abstol)
    v1 = copy(model.landau.dist.particles.v)
    ΔS = Sh(v1, model) - Sh(vn, model)
    @test ΔS < 0
    @test maximum(abs, v1 .- vn) ≈ 5.863e-4 atol = 1e-7
    @test abs(ΔS) ≈ 1.973e-4 atol = 1e-7
    g = similar(vn)
    for nq in (2, 5, 8)
        VM.discrete_gradient!(g, vn, v1, model, AVFGradient(nq))
        @test abs(ΔS - dot(v1 .- vn, g)) ≤ 1e-10 * abs(ΔS)
    end

    # the control: the midpoint gradient breaks the bound
    VM.discrete_gradient!(g, vn, v1, model, AVFGradient(1))
    @test abs(ΔS - dot(v1 .- vn, g)) > 1e-10 * abs(ΔS)
end

@testset "the identity falls with nq" begin
    model = fixture()
    v0 = copy(model.landau.dist.particles.v)
    d = randn(StableRNG(7), 2, 144)
    d ./= norm(d, Inf)
    v1 = v0 .+ 1e-1 .* d
    ΔS = Sh(v1, model) - Sh(v0, model)
    g = similar(v0)
    rel = map((1, 2, 3, 5)) do nq
        VM.discrete_gradient!(g, v0, v1, model, AVFGradient(nq))
        abs(ΔS - dot(v1 .- v0, g)) / abs(ΔS)
    end
    @test rel[2] ≤ rel[1] / 5
    @test rel[3] ≤ rel[2] / 5
    @test rel[4] ≤ 1e-5
end

@testset "Δv = 0" begin
    model = fixture()
    v0 = copy(model.landau.dist.particles.v)
    g0 = similar(v0)
    VM.entropy_gradient!(g0, v0, model)
    scale = maximum(abs, g0)
    g = similar(v0)
    for nq in (5, 1)
        VM.discrete_gradient!(g, v0, copy(v0), model, AVFGradient(nq))
        @test maximum(abs, g .- g0) ≤ 8 * eps() * scale
    end
    VM.discrete_gradient!(g, v0, copy(v0), model, AVFGradient(1))
    @test g == g0
end

@testset "one AVF step" begin
    model = fixture()
    vn = copy(model.landau.dist.particles.v)
    w = weights(model)
    Δt, abstol = 2e-3, 1e-12
    r = discrete_gradient_step!(model, Δt; abstol, gradient = AVFGradient(5))
    v1 = copy(model.landau.dist.particles.v)
    @test r.iterations == 5
    @test r.residual ≤ abstol
    @test v1 != vn
    @test all(abs.(momentum(v1, w) .- momentum(vn, w)) .≤ momentum_bound(vn, w))
    @test abs(energy(v1, w) - energy(vn, w)) ≤ energy_bound(vn, v1, Δt, abstol, w)

    # the step is the one written out from the rule's own nodes and weights
    written = written_step(model, vn, Δt, abstol; gradient = AVFGradient(5))
    @test maximum(abs, written .- v1) ≤ abstol

    # the control: 𝔾(vⁿ) in place of 𝔾 at the midpoint breaks the energy bound at both steps, by
    # the square of the ratio of the step sizes
    dE = map((2e-3, 1e-3)) do Δt
        control = fixture()
        vc0 = copy(control.landau.dist.particles.v)
        vc1 = written_step(control, vc0, Δt, abstol;
            gradient = AVFGradient(5), G_at_old = true)
        @test abs(energy(vc1, w) - energy(vc0, w)) > energy_bound(vc0, vc1, Δt, abstol, w)
        abs(energy(vc1, w) - energy(vc0, w))
    end
    @test dE[1] / dE[2] ≈ 4 rtol = 0.25
end

@testset "the domain of the accepted state" begin
    # The nodes of the five-point rule lie below 0.954 of the segment, so none of them leaves the
    # domain of the spline basis while the accepted state does: the projection at the accepted
    # state throws, as the entropy of `v1` does for the Gonzalez rule, and the particles stay as
    # they were.
    v0, w = fixture_particles(Float64)
    L = maximum(abs, v0) + 1e-3
    model = build(v0, w; L)
    pv = model.landau.dist.particles.v
    f = similar(pv)
    gap, frac = 1e-5, 1.03
    found = false
    for a in axes(pv, 2), c in 1:2

        s = sign(v0[c, a])
        pv .= v0
        pv[c, a] = s * (L - gap)
        VM.collisional_vectorfield!(f, pv, nothing, model)
        s * f[c, a] ≥ 1e-3 || continue
        Δt = frac * gap / (s * f[c, a])
        for rule in (AVFGradient(5), GonzalezGradient())
            pv .= v0
            pv[c, a] = s * (L - gap)
            before = copy(pv)
            @test throws_with(
                () -> discrete_gradient_step!(model, Δt; abstol = 1e-12,
                    gradient = rule),
                DomainError, "domain")
            @test pv == before
        end
        found = true
        break
    end
    @test found
end

@testset "the relaxed state, both rules" begin
    model = fixture()
    w = weights(model)
    for _ in 1:250
        discrete_gradient_step!(model, 8e-3; abstol = 1e-12)
    end
    vn = copy(model.landau.dist.particles.v)
    @test 1.14 < anisotropy(vn, w) < 1.143

    # 16 starts: the relaxed state, and 15 copies with each entry moved by −2 to 2 ulps
    starts = map(0:15) do s
        v = copy(vn)
        s == 0 && return v
        rng = StableRNG(s)
        for i in eachindex(v)
            k = floor(Int, 5 * rand(rng)) - 2
            for _ in 1:abs(k)
                v[i] = k > 0 ? nextfloat(v[i]) : prevfloat(v[i])
            end
        end
        return v
    end

    # the largest relative increment of the evaluations 11 to 40 of each rule, from each start:
    # at most 6.68e-13 for the averaged-vector-field rule, and 0 on 10 of the 16 starts
    tail(x) = maximum(x[11:40])
    a_max = maximum(v -> tail(increments(model, v, 8e-3, AVFGradient(5))), starts)
    g_max = map(v -> tail(increments(model, v, 8e-3, GonzalezGradient())), starts)
    @test a_max ≤ 1.5e-12
    # the median of the 16 values: the mean of the eighth and the ninth
    # measured 2.77e-11 for the Gonzalez rule, with the 16 starts from 1.50e-11 to 4.71e-11
    gs = sort(g_max)
    @test (gs[8] + gs[9]) / 2 ≥ 1e-11

    # one step from the relaxed state with each rule, and the residual of the state it returns
    # against the fixed point of the averaged-vector-field map there
    ∇̄S = similar(vn)
    f = similar(vn)
    function map_residual(gradient)
        m = fixture()
        m.landau.dist.particles.v .= vn
        discrete_gradient_step!(m, 8e-3; abstol = 1e-12, gradient)
        v1 = copy(m.landau.dist.particles.v)
        VM.discrete_gradient!(∇̄S, vn, v1, m, AVFGradient(5))
        VM.apply_G!(f, (vn .+ v1) ./ 2, ∇̄S, m)
        return maximum(abs, vn .+ 8e-3 * m.landau.ν .* f .- v1)
    end
    @test map_residual(AVFGradient(5)) ≤ 1e-12
    @test map_residual(GonzalezGradient()) > 1e-10
end

@testset "twenty AVF steps" begin
    model = fixture()
    w = weights(model)
    Δt, abstol = 2e-3, 1e-12
    v0 = copy(model.landau.dist.particles.v)
    Ebound = 0.0
    smallest, worst, worst_control = Inf, 0.0, 0.0
    g = similar(v0)
    for n in 1:20
        vn = copy(model.landau.dist.particles.v)
        discrete_gradient_step!(model, Δt; abstol, gradient = AVFGradient(5))
        v1 = copy(model.landau.dist.particles.v)
        Ebound += energy_bound(vn, v1, Δt, abstol, w)

        ΔS = Sh(v1, model) - Sh(vn, model)
        @test ΔS < 0
        smallest = min(smallest, -ΔS)
        VM.discrete_gradient!(g, vn, v1, model, AVFGradient(5))
        worst = max(worst, abs(ΔS - dot(v1 .- vn, g)) / abs(ΔS))
        @test abs(ΔS - dot(v1 .- vn, g)) ≤ 1e-7 * abs(ΔS)

        # the control: the midpoint gradient is one node of the rule, and breaks the bound
        VM.discrete_gradient!(g, vn, v1, model, AVFGradient(1))
        worst_control = max(worst_control, abs(ΔS - dot(v1 .- vn, g)) / abs(ΔS))
    end
    v20 = model.landau.dist.particles.v
    @test smallest > 0
    @test worst > 0
    @test worst_control > 1e-7
    @test all(abs.(momentum(v20, w) .- momentum(v0, w)) .≤ 20 * momentum_bound(v0, w))
    @test abs(energy(v20, w) - energy(v0, w)) ≤ Ebound
    @test anisotropy(v20, w) < anisotropy(v0, w)
    @test 1.211 < anisotropy(v20, w) < 1.214
end

@testset "Float32, nq = 0 and the hot path" begin
    v32, w32 = fixture_particles(Float32)
    m32 = build(v32, w32)
    v0 = copy(m32.landau.dist.particles.v)
    d = randn(StableRNG(5), Float32, 2, 144)
    d ./= norm(d, Inf)
    v1 = v0 .+ 1.0f-3 .* d
    ΔS = Sh(v1, m32) - Sh(v0, m32)
    g = similar(v0)
    gt = similar(v0)
    ref = similar(v0)
    for nq in (2, 5)
        rule = AVFGradient(nq)
        VM.discrete_gradient!(g, v0, v1, m32, rule)
        @test eltype(g) == Float32
        @test abs(ΔS - dot(v1 .- v0, g)) ≤ 1e-5

        # the rule stays in Float32: its sum is the one written out in Float32, bit for bit
        fill!(ref, 0)
        for i in 1:nq
            VM.entropy_gradient!(gt, v0 .+ Float32(rule.nodes[i]) .* (v1 .- v0), m32)
            ref .+= Float32(rule.weights[i]) .* gt
        end
        @test g == ref
    end

    Δt, abstol = 2.0f-3, 1.0f-4
    r = discrete_gradient_step!(m32, Δt; abstol, gradient = AVFGradient(5))
    y1 = copy(m32.landau.dist.particles.v)
    @test eltype(y1) == Float32
    @test r.iterations == 2
    @test r.residual ≤ abstol
    @test abs(energy(y1, w32) - energy(v0, w32)) ≤ energy_bound(v0, y1, Δt, abstol, w32)

    # an nq below 1 is refused, and an nq that is not an integer has no method
    for bad in (0, -1)
        @test throws_with(() -> AVFGradient(bad), ArgumentError, "nq")
    end
    @test_throws MethodError AVFGradient(2.5)

    # the hot path: no allocation after a warm call and a concrete return, in both precisions
    m64 = build(Float64.(v0), Float64.(w32))
    for (model, v) in ((m32, v0), (m64, Float64.(v0)))
        x = v .+ eltype(v)(1e-3) .* eltype(v).(d)
        y = similar(v)
        for nq in (1, 5)
            rule = AVFGradient(nq)
            @test bytes5(VM.discrete_gradient!, y, v, x, model, rule) == 0
            @inferred VM.discrete_gradient!(y, v, x, model, rule)
        end
        @test bytes5(copy_first_5, y, v, x, model, AVFGradient(5)) > 0
    end
end

@testset "a large nq" begin
    model = fixture()
    vn = copy(model.landau.dist.particles.v)
    r = discrete_gradient_step!(model, 2e-3; abstol = 1e-12, gradient = AVFGradient(40))
    @test r.iterations == 5
    @test r.residual ≤ 1e-12
    @test copy(model.landau.dist.particles.v) != vn

    d = randn(StableRNG(5), 2, 144)
    d ./= norm(d, Inf)
    v1 = vn .+ 1e-3 .* d
    g = similar(vn)
    VM.discrete_gradient!(g, vn, v1, model, AVFGradient(40))
    ΔS = Sh(v1, model) - Sh(vn, model)
    @test abs(ΔS - dot(v1 .- vn, g)) ≤ 1e-10 * abs(ΔS)
end

@testset "the labels, the cost and the pages" begin
    d = string(@doc AVFGradient)
    for word in ("eq:disc_grad_props", "quadrature error", "O(N²)", "nq", "midpoint")
        @test occursin(word, d)
    end
    # the measured cost of one call, against the Gonzalez rule
    @test occursin("271", d) && occursin("2.24", d)
    ds = string(@doc discrete_gradient_step!)
    for word in ("AVFGradient", "GonzalezGradient", "quadrature error")
        @test occursin(word, ds)
    end
    @test :AVFGradient in names(VlasovMethods)
end
