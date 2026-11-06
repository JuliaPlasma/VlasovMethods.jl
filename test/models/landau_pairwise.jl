using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, TensorProductBasis,
                     TensorProductQuadrature
using StableRNGs
using StaticArrays: SVector
using Test
using VlasovMethods

const VM = VlasovMethods

# A model of the particles `v` (2 × N) with the weights `w`, on quadratic splines of 8 × 8 cells
# of [-4, 4]².
function build(v::AbstractMatrix{T}, w::AbstractVector{T}; γ = -3, ν = 1.0) where {T}
    N = size(v, 2)
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pd = ParticleDistribution(1, 2, ParticleList(zeros(T, 4, N); variables = vars))
    pd.particles.v .= v
    pd.particles.w[1, :] .= w
    axis = BSplineBasis(UniformMesh(8, T(-4) .. T(4)), 2, Free())
    B = TensorProductBasis((axis, axis))
    sd = SplineDistribution(1, 2, B, TensorProductQuadrature(B), zeros(T, size(B)...))
    return PairwiseLandau(pd, CollisionEntropy(sd); γ, ν)
end

# The fixture: a 12 × 12 lattice on [-4 + 1/3, 4 - 1/3]², each velocity moved by
# 0.02 randn, with the unequal weights exp(-x²/(2·1.3²) - y²/(2·0.8²)) + 0.08 of the lattice
# point, normalised to Σ w = 1, so that T_x/T_y = 1.2146 and f_s > 0 at every node.
function fixture_particles(::Type{T}) where {T}
    rng = StableRNG(1234)
    xs = range(-4 + 1 / 3, 4 - 1 / 3; length = 12)
    v = zeros(T, 2, 144)
    w = zeros(144)
    m = 0
    for y in xs, x in xs

        m += 1
        v[1, m] = T(x + 0.02 * randn(rng))
        v[2, m] = T(y + 0.02 * randn(rng))
        w[m] = exp(-x^2 / (2 * 1.3^2) - y^2 / (2 * 0.8^2)) + 0.08
    end
    return v, T.(w ./ sum(w))
end
fixture(::Type{T} = Float64; kw...) where {T} = build(fixture_particles(T)...; kw...)

weights(model) = collect(vec(model.landau.dist.particles.w))
momentum(v, w) = v * w
energy(v, w) = sum(w .* vec(sum(abs2, v; dims = 1))) / 2
function anisotropy(v, w)
    u = momentum(v, w) ./ sum(w)
    return sum(w .* (v[1, :] .- u[1]) .^ 2) / sum(w .* (v[2, :] .- u[2]) .^ 2)
end

# ∂S_h/∂v_{c,α} = Σ_k K_c[k, α] J_k (eq:entropy_derivative), and S_h, from the pieces of the
# appendix model that the pairwise model holds.
function gradS(v::AbstractMatrix{T}, model) where {T}
    landau = model.landau
    cache = landau.cache[T]
    sd = cache.sdist
    VM._projection!(cache.rhs, cache.vals, cache.mwork, v, landau.dist, sd)
    VM.compute_K!(cache.K1, cache.K2, v, sd, landau)
    VM.compute_J!(cache.J, sd, landau)
    return permutedims(hcat(transpose(cache.K1) * cache.J, transpose(cache.K2) * cache.J))
end
function Sh(v::AbstractMatrix{T}, model) where {T}
    landau = model.landau
    cache = landau.cache[T]
    VM._projection!(cache.rhs, cache.vals, cache.mwork, v, landau.dist, cache.sdist)
    return CollisionEntropy(cache.sdist)()
end

function kernelof(model)
    (a, b) -> VM.kernel(SVector(a[1], a[2]), SVector(b[1], b[2]), model.landau)
end

# eq:G_semi-disc, written out: (𝔾 x)_γ = −Σ_{β≠γ} U(v_γ, v_β) (w_β/w_γ x_γ − x_β); with
# `transposed`, the weight ratio w_γ/w_β of the control.
function G_written(v, x, w, U; transposed = false)
    y = zero(x)
    for γ in axes(v, 2), β in axes(v, 2)

        β == γ && continue
        r = transposed ? w[γ] / w[β] : w[β] / w[γ]
        y[:, γ] .-= U(v[:, γ], v[:, β]) * (r .* x[:, γ] .- x[:, β])
    end
    return y
end

# eq:eom_line3, written out: v̇_γ = −ν Σ_{α≠γ} w_α U(v_γ, v_α) (ψ_γ − ψ_α), ψ_α = ∂S_h/∂v_α / w_α.
function eom_line3(v, model)
    w = weights(model)
    ψ = gradS(v, model) ./ w'
    U = kernelof(model)
    out = zero(v)
    for γ in axes(v, 2), α in axes(v, 2)

        α == γ && continue
        out[:, γ] .-= w[α] .* (U(v[:, γ], v[:, α]) * (ψ[:, γ] .- ψ[:, α]))
    end
    return model.landau.ν .* out
end

# One step of eq:disc_grad_update_rule, written out with the Gonzalez gradient: from v_0 = vⁿ,
# f_k = ν 𝔾(m_k) ∇̄S(vⁿ, v_k), m_k = (vⁿ + v_k)/2, accepted at ‖v_k − vⁿ − Δt f_k‖∞ ≤ abstol.
# The controls: `G_at_old` takes 𝔾(vⁿ), `accept_iterate` returns v_k in place of vⁿ + Δt f_k.
function written_step(model, vn, Δt, abstol; G_at_old = false, accept_iterate = false)
    w = weights(model)
    S0 = Sh(vn, model)
    vk = copy(vn)
    for k in 0:99
        vm = (vn .+ vk) ./ 2
        g = gradS(vm, model)
        Δ = vk .- vn
        n2 = sum(abs2, Δ)
        iszero(n2) || (g = g .+ Δ .* ((Sh(vk, model) - S0 - dot(Δ, g)) / n2))
        f = model.landau.ν .* G_written(G_at_old ? vn : vm, g, w, kernelof(model))
        vnew = vn .+ Δt .* f
        maximum(abs, vk .- vnew) ≤ abstol && return accept_iterate ? vk : vnew
        vk = vnew
    end
    error("the written-out step did not converge")
end

# The bound on the energy change of one step. With 𝔾 and the iterate at the same midpoint m_k,
# Σ_γ w_γ m_γ·(𝔾(m) x)_γ = 0 for every x, so Σ w f_k·m_k = 0. The accepted state is
# vⁿ⁺¹ = vⁿ + Δt f_k, and (vⁿ⁺¹ + vⁿ)/2 = m_k − r_k/2 with r_k = v_k − vⁿ − Δt f_k, so
# E(vⁿ⁺¹) − E(vⁿ) = Δt Σ w f_k·(vⁿ⁺¹ + vⁿ)/2 = −(Δt/2) Σ w f_k·r_k, at most
# (Δt/2) ‖r_k‖∞ Σ_α w_α ‖f_α‖₁ ≤ (Δt/2) abstol Σ_α w_α ‖f_α‖₁; 64 eps E covers the rounding.
function energy_bound(vn, v1, Δt, abstol, w)
    f = (v1 .- vn) ./ Δt
    return Δt / 2 * abstol * sum(w .* vec(sum(abs, f; dims = 1))) +
           64 * eps(eltype(vn)) * energy(vn, w)
end
# The bound on each component of the momentum change of one step: rounding only.
momentum_bound(v, w) = 64 * eps(eltype(v)) * sum(w .* vec(sqrt.(sum(abs2, v; dims = 1))))

function throws_with(f, ::Type{E}, msg) where {E}
    err = try
        f()
        nothing
    catch e
        e
    end
    return err isa E && occursin(msg, sprint(showerror, err))
end

# The bytes of one call after a warm call, through a barrier whose arguments have concrete types.
function bytes4(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    f(a, b, c, d)
    return @allocated f(a, b, c, d)
end
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
    @test throws_with(
        () -> discrete_gradient_step!(model, 2e-3; abstol = 1e-12), DomainError,
        msg)
    @test model.landau.dist.particles.v == before

    # one equal component is no equal velocity
    v = copy(fixture().landau.dist.particles.v)
    v[1, 42] = v[1, 30]
    y = similar(v)
    VM.apply_G!(y, v, x, build(v, weights(model)))
    @test all(isfinite, y)
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

@testset "one particle, no particle, NaN" begin
    one_particle = build(reshape([0.1, -0.2], 2, 1), [1.0])
    v̇ = ones(2, 1)
    VM.collisional_vectorfield!(v̇, [0.1 -0.2]', nothing, one_particle)
    @test iszero(v̇)
    @test discrete_gradient_step!(one_particle, 2e-3; abstol = 1e-12) ==
          (iterations = 1, residual = 0.0)
    @test one_particle.landau.dist.particles.v == reshape([0.1, -0.2], 2, 1)

    none = build(zeros(2, 0), zeros(0))
    @test_throws DomainError VM.collisional_vectorfield!(zeros(2, 0), zeros(2, 0), nothing, none)
    @test_throws DomainError discrete_gradient_step!(none, 2e-3; abstol = 1e-12)

    model = fixture()
    model.landau.dist.particles.v[1, 5] = NaN
    before = copy(model.landau.dist.particles.v)
    @test_throws DomainError VM.collisional_vectorfield!(similar(before), before, nothing, model)
    @test_throws DomainError discrete_gradient_step!(model, 2e-3; abstol = 1e-12)
    @test isequal(model.landau.dist.particles.v, before)
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

    Δt, abstol = 2.0f-3, 1.0f-4
    r = discrete_gradient_step!(m32, Δt; abstol)
    v1 = copy(m32.landau.dist.particles.v)
    @test r.residual isa Float32
    @test r.residual ≤ abstol
    @test eltype(v1) == Float32
    @test v1 != v32
    @test all(abs.(momentum(v1, w32) .- momentum(v32, w32)) .≤ momentum_bound(v32, w32))
    @test abs(energy(v1, w32) - energy(v32, w32)) ≤ energy_bound(v32, v1, Δt, abstol, w32)

    for (model, v) in ((m32, v32), (m64, Float64.(v32)))
        x = ones(eltype(v), size(v))
        y = similar(v)
        @test bytes4(VM.apply_G!, y, v, x, model) == 0
        @test bytes4(VM.collisional_vectorfield!, y, v, nothing, model) == 0
        @test bytes4(copy_first, y, v, x, model) > 0
        @inferred VM.apply_G!(y, v, x, model)
        @inferred VM.collisional_vectorfield!(y, v, nothing, model)
    end
end
