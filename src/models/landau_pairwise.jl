@noinline function _throw_weight(i, wi)
    throw(ArgumentError("the weight of particle $i is $wi; the pairwise scheme divides by " *
                        "the weights, which must be positive and finite"))
end
function _check_weights(w)
    i = findfirst(x -> !(isfinite(x) && x > 0), w)
    isnothing(i) || _throw_weight(i, w[i])
    return nothing
end

@doc raw"""
    PairwiseLandau(dist, ent; ν = 1.0, γ = -3)

The particle-pair Landau scheme of the Landau manuscript's main text: the gradient system

```math
\dot{v} = \nu \, \mathbb{G}(v) \, \nabla S_h(v)
```

of `eq:landau_sd_gradient_form`, with the matrix ``\mathbb{G}`` of `eq:G_semi-disc`,

```math
\big( \mathbb{G}(v) \, x \big)_\gamma
  = - \sum_{\beta \ne \gamma} U(v_\gamma, v_\beta)
      \left( \frac{w_\beta}{w_\gamma} \, x_\gamma - x_\beta \right) ,
```

and the gradient ``\partial S_h / \partial v_\alpha = w_\alpha \sum_k \mathbb{L}_k \nabla
\varphi_k(v_\alpha)`` of the entropy ``S_h = \int f_s \log |f_s| \, dv`` of the particles'
projection. Written out, this is `eq:eom_line3`,

```math
\dot{v}_\gamma = - \nu \sum_{\alpha \ne \gamma} w_\alpha \, U(v_\gamma, v_\alpha)
  \left( \psi_\gamma - \psi_\alpha \right) , \qquad
\psi_\alpha = \sum_k \mathbb{L}_k \nabla \varphi_k(v_\alpha) .
```

``\mathbb{G}`` is symmetric and negative semi-definite, so ``S_h`` decreases; ``\sum_\gamma
w_\gamma \mathbb{G}_{\gamma\sigma} = 0`` conserves the momentum ``\sum_\alpha w_\alpha v_\alpha``
and ``\sum_\gamma w_\gamma v_\gamma \cdot \mathbb{G}_{\gamma\sigma} = 0`` the energy
``\frac{1}{2} \sum_\alpha w_\alpha |v_\alpha|^2``. For distinct velocities that do not all lie
on one line the kernel of ``\mathbb{G}`` has the dimension ``d + 1 = 3``. Only two velocity
dimensions are supported.

The model holds a `Landau` of the same arguments in its field `landau`, and takes from it
the particle distribution, the projection onto the spline distribution of `ent`, the matrices
``\mathbb{K}`` of [`compute_K!`](@ref) and ``\mathbb{L}_k`` of [`compute_J!`](@ref), the
[`kernel`](@ref) of the exponent `γ`, its check of `γ` and the buffers of its cache. The weights
are row 1 of `dist.particles.w`. ``\mathbb{G}`` divides by them, so a weight that is not positive
and finite throws an `ArgumentError`, in the constructor and in each [`apply_G!`](@ref), which
reads the weights of the distribution as they are at the call.

The particle velocities must be distinct: ``U(v, v)`` is not defined, and two particles of equal
velocities throw a `DomainError`. For ``\gamma < -2``, where ``U`` diverges as ``|u| \to 0``, so
do two different velocities whose squared distance is below a threshold of the element type `T`
and of ``\gamma``, at which ``U`` in `T` stops being finite and not zero: for ``\gamma = -3`` a
distance below ``(2 / \mathtt{floatmax}(T)^{2/3})^{1/2}``, `2.0e-13` in `Float32` and `2.5e-103`
in `Float64`. For ``\gamma \ge -2`` ``U`` stays finite, and only equal velocities throw.

# Cost

One product ``\mathbb{G}(v) \, x`` of [`apply_G!`](@ref) costs `O(N²)` kernel evaluations,
one per ordered pair of particles. On the fixture of `test/helpers/pairwise_landau.jl` widened to
a lattice of `N` particles, quadratic splines on 8 × 8 cells (`M = 100`, `Q = 576`), `Float64`,
the minimum of 5 calls after a warm call in a cold process: `apply_G!` takes 0.36 ms at
`N = 400` and 5.75 ms at `N = 1600`, 2.25 ns per pair, at a load of about 10; the gradient
``\nabla S_h`` takes 0.08 ms and 0.27 ms, and ``S_h`` 0.05 ms and 0.10 ms, at a load of about 11.

The time integrator is [`discrete_gradient_step!`](@ref).
"""
struct PairwiseLandau{LT <: Landau} <: VlasovModel
    landau::LT

    function PairwiseLandau(dist::DistributionFunction, ent::Entropy; ν = 1.0, γ = -3)
        _check_weights(view(dist.particles.w, 1, :))
        landau = Landau(dist, ent; γ, ν)
        new{typeof(landau)}(landau)
    end
end

@doc raw"""
    entropy_gradient!(∇S, v, model::PairwiseLandau)

The gradient ``\partial S_h / \partial v_{c,\alpha} = \sum_k \mathbb{K}_c[k, \alpha] \,
\mathbb{L}_k`` of the entropy ``S_h = \int f_s \log |f_s| \, dv`` of the projection of the
particles at the velocities `v` (`2 × N`), written into `∇S` (`2 × N`). The projection,
``\mathbb{K}`` and ``\mathbb{L}_k`` are those of `model.landau`, with the buffers of its cache of
the element type of `v`; the projection overwrites the coefficients of that cache's spline
distribution. `v` or `∇S` of another size than `2 × N`, with `N` the number of particles of the
model, throws a `DimensionMismatch`.
"""
function entropy_gradient!(∇S, v::AbstractMatrix{T}, model::PairwiseLandau) where {T}
    _check_size(:v, v, model)
    _check_size(:∇S, ∇S, model)
    landau = model.landau
    cache = landau.cache[T]
    sdist = cache.sdist
    _projection!(cache.rhs, cache.vals, cache.mwork, v, landau.dist, sdist)
    compute_K!(cache.K1, cache.K2, v, sdist, landau)
    compute_J!(cache.J, sdist, landau)
    mul!(view(∇S, 1, :), transpose(cache.K1), cache.J)
    mul!(view(∇S, 2, :), transpose(cache.K2), cache.J)
    return ∇S
end

# The arrays of `entropy_gradient!`, `apply_G!` and the vector field are `2 × N`: one row per
# velocity component and one column per particle of the model.
@noinline function _throw_size(name, sz, N)
    throw(DimensionMismatch("$name has size $sz, and the model has $N particles of 2 " *
                            "velocity components"))
end
function _check_size(name, a, model::PairwiseLandau)
    N = size(model.landau.dist.particles.w, 2)
    size(a) == (2, N) || _throw_size(name, size(a), N)
    return nothing
end

# S_h of the projection of the particles at the velocities `v`, on the cache's spline
# distribution, evaluated right after its own projection.
function _entropy(v::AbstractMatrix{T}, model::PairwiseLandau) where {T}
    landau = model.landau
    cache = landau.cache[T]
    _projection!(cache.rhs, cache.vals, cache.mwork, v, landau.dist, cache.sdist)
    return CollisionEntropy(cache.sdist)()
end

# `apply_G!` meets each repeated pair first with `γ < β`: the pair of the smallest index with a
# repeated velocity, and the smallest of its partners.
@noinline function _throw_equal_velocities(γ, β)
    throw(DomainError((γ, β),
        "particles $γ and $β have equal velocities, where the kernel U(v, v) of the " *
        "pairwise Landau scheme is not defined"))
end

@noinline function _throw_close_velocities(γ, β, n2, τ)
    throw(DomainError((γ, β),
        "particles $γ and $β have velocities too close for the kernel U of the pairwise " *
        "Landau scheme in $(typeof(τ)): |v_$γ - v_$β|² = $n2 is below $τ"))
end

# The threshold τ of the squared distance n2 = |u|², u = v_α - v_β, at and above which the
# `kernel` of `landau` in `T` is finite and not zero. With F = floatmax(T) and e = γ + 2:
# - e ≥ 0: τ = 0. U tends to a finite value or to zero as |u| → 0, and stays finite.
# - e = -1, the Coulomb branch: 1/|u|³ = n2^(-3/2) is finite for n2 > F^(-2/3), and then so are
#   1/|u| and the entries, each at most 1/|u| in size; U₁₁ + U₂₂ = 1/|u| > 0. τ = 2 F^(-2/3)
#   leaves n2^(-3/2) ≤ F / 2^(3/2).
# - the other e < 0, the hypot branch: |u|^e is finite for |u| > F^(1/e), that is for
#   n2 > F^(2/e), and the larger diagonal entry is at least |u|^e / 2. τ = 2 F^(2/e) leaves
#   |u|^e ≤ 2^(e/2) F; where τ exceeds floatmin(T), |e| ≥ 2 and so |u|^e ≤ F / 2. Where it does
#   not, τ = floatmin(T), above the bound F^(2/e).
# τ > 0, so an n2 that underflows to zero for two different velocities is refused too.
function _close_threshold(::Type{T}, landau::Landau) where {T}
    e = T(landau.γ) + 2
    e < 0 || return zero(T)
    e == -1 && return 2 * floatmax(T)^(T(-2) / 3)
    return max(floatmin(T), 2 * floatmax(T)^(2 / e))
end

@doc raw"""
    apply_G!(y, v, x, model::PairwiseLandau)

The product ``y = \mathbb{G}(v) \, x`` of the matrix of `eq:G_semi-disc` at the velocities `v`
with `x`, all three `2 × N`:

```math
y_\gamma = - \sum_{\beta \ne \gamma} U(v_\gamma, v_\beta)
  \left( \frac{w_\beta}{w_\gamma} \, x_\gamma - x_\beta \right) ,
```

with ``U`` the [`kernel`](@ref) of `model.landau` and the weights ``w`` of row 1 of
`model.landau.dist.particles.w`. The pair ``\beta = \gamma`` is skipped by its index. Two
particles ``\gamma \ne \beta`` of equal velocities throw a `DomainError` that names both, the
smaller index first. For ``\gamma < -2`` so do two different velocities whose squared distance is
below the threshold of the element type of `v` at which ``U`` stops being finite and not zero.
`y` must not alias `v` or `x`. An array of another size than `2 × N`, with `N` the number of
particles of the model, throws a `DimensionMismatch`, and a weight that is not positive and
finite an `ArgumentError`.

The cost is ``O(N^2)`` kernel evaluations, one per ordered pair; the [`PairwiseLandau`](@ref)
docstring gives measured times. A call allocates nothing.
"""
function apply_G!(y, v::AbstractMatrix{T}, x, model::PairwiseLandau) where {T}
    _check_size(:y, y, model)
    _check_size(:v, v, model)
    _check_size(:x, x, model)
    landau = model.landau
    w = landau.dist.particles.w
    _check_weights(view(w, 1, :))
    τ = _close_threshold(T, landau)
    for γ in axes(v, 2)
        vγ = SVector(v[1, γ], v[2, γ])
        xγ = SVector(x[1, γ], x[2, γ])
        wγ = w[1, γ]
        acc = zero(SVector{2, T})
        for β in axes(v, 2)
            β == γ && continue
            vβ = SVector(v[1, β], v[2, β])
            vγ == vβ && _throw_equal_velocities(γ, β)
            u = vγ - vβ
            n2 = u[1] * u[1] + u[2] * u[2]
            n2 < τ && _throw_close_velocities(γ, β, n2, τ)
            U = kernel(vγ, vβ, landau)
            acc += U * (w[1, β] / wγ * xγ - SVector(x[1, β], x[2, β]))
        end
        y[1, γ] = -acc[1]
        y[2, γ] = -acc[2]
    end
    return y
end

@doc raw"""
    collisional_vectorfield!(v̇, v, params, model::PairwiseLandau)

The vector field ``\dot{v} = \nu \, \mathbb{G}(v) \, \nabla S_h(v)`` of
`eq:landau_sd_gradient_form`, which is `eq:eom_line3`, at the velocities `v`, written into `v̇`
(both `2 × N`). `params` is not read. The gradient goes into the buffer `v̇` of
`model.landau.cache[T]`, with `T` the element type of `v`, so a call allocates nothing; `ν` is
converted to `T`.

The projection throws a `DomainError` for a velocity outside the domain of the spline basis,
`NaN` and `Inf` included, and [`compute_J!`](@ref) for a projection that is zero at a quadrature
node, as for no particle; [`apply_G!`](@ref) throws for two equal velocities and for a weight
that is not positive and finite. `v` or `v̇` of another size than `2 × N`, with `N` the number
of particles of the model, throws a `DimensionMismatch` that names it.
"""
function collisional_vectorfield!(v̇::AbstractArray{ST}, v::AbstractArray{ST}, params,
        model::PairwiseLandau) where {ST}
    _check_size(:v̇, v̇, model)
    ∇S = model.landau.cache[ST].v̇
    entropy_gradient!(∇S, v, model)
    apply_G!(v̇, v, ∇S, model)
    v̇ .*= ST(model.landau.ν)
    return nothing
end
