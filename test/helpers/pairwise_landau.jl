# The fixture and the written-out forms of the particle-pair Landau scheme, shared by
# `test/models/landau_pairwise.jl` and `test/methods/discrete_gradient.jl`. The including file
# loads LinearAlgebra, ParticleMethods, SimpleSplines, StableRNGs, StaticArrays and
# VlasovMethods.

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
