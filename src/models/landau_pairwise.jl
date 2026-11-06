@noinline function _throw_weight(i, wi)
    throw(ArgumentError("the weight of particle $i is $wi; the pairwise scheme divides by " *
                        "the weights, which must be positive and finite"))
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
``\frac{1}{2} \sum_\alpha w_\alpha |v_\alpha|^2``. For distinct velocities the kernel of
``\mathbb{G}`` has the dimension ``d + 1 = 3``. Only two velocity dimensions are supported.

The model holds a `Landau` of the same arguments in its field `landau`, and takes from it
the particle distribution, the projection onto the spline distribution of `ent`, the matrices
``\mathbb{K}`` of [`compute_K!`](@ref) and ``\mathbb{L}_k`` of [`compute_J!`](@ref), the
[`kernel`](@ref) of the exponent `γ`, its check of `γ` and the buffers of its cache. The weights
are row 1 of `dist.particles.w`. ``\mathbb{G}`` divides by them, so a weight that is not positive
and finite throws an `ArgumentError`.

The particle velocities must be distinct: ``U(v, v)`` is not defined, and two particles of equal
velocities throw a `DomainError`.

# Cost

One product ``\mathbb{G}(v) \, x`` of [`apply_G!`](@ref) costs `O(N²)` kernel evaluations,
one per ordered pair of particles. On the fixture of `test/models/landau_pairwise.jl` widened to
a lattice of `N` particles, quadratic splines on 8 × 8 cells (`M = 100`, `Q = 576`), `Float64`,
the minimum of 5 calls after a warm call in a cold process at a load of about 11: `apply_G!`
takes 0.35 ms at `N = 400` and 5.7 ms at `N = 1600`, 2.2 ns per pair; the gradient
``\nabla S_h`` takes 0.08 ms and 0.27 ms, and ``S_h`` 0.05 ms and 0.10 ms.

The time integrator is [`discrete_gradient_step!`](@ref).
"""
struct PairwiseLandau{LT <: Landau} <: VlasovModel
    landau::LT

    function PairwiseLandau(dist::DistributionFunction, ent::Entropy; ν = 1.0, γ = -3)
        w = view(dist.particles.w, 1, :)
        i = findfirst(x -> !(isfinite(x) && x > 0), w)
        isnothing(i) || _throw_weight(i, w[i])
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
distribution. `v` or `∇S` of another number of columns than the model has particles throws a
`DimensionMismatch`.
"""
function entropy_gradient!(∇S, v::AbstractMatrix{T}, model::PairwiseLandau) where {T}
    _check_columns(:v, v, model)
    _check_columns(:∇S, ∇S, model)
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

# The arrays of `entropy_gradient!` and `apply_G!` have one column per particle of the model.
@noinline function _throw_columns(name, n, N)
    throw(DimensionMismatch("$name has $n columns, and the model has $N particles"))
end
function _check_columns(name, a, model::PairwiseLandau)
    N = size(model.landau.dist.particles.w, 2)
    size(a, 2) == N || _throw_columns(name, size(a, 2), N)
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
smaller index first. `y` must not alias `v` or `x`. An array of another number of columns than
the model has particles throws a `DimensionMismatch`.

The cost is ``O(N^2)`` kernel evaluations, one per ordered pair; the [`PairwiseLandau`](@ref)
docstring gives measured times. A call allocates nothing.
"""
function apply_G!(y, v::AbstractMatrix{T}, x, model::PairwiseLandau) where {T}
    _check_columns(:y, y, model)
    _check_columns(:v, v, model)
    _check_columns(:x, x, model)
    landau = model.landau
    w = landau.dist.particles.w
    for γ in axes(v, 2)
        vγ = SVector(v[1, γ], v[2, γ])
        xγ = SVector(x[1, γ], x[2, γ])
        wγ = w[1, γ]
        acc = zero(SVector{2, T})
        for β in axes(v, 2)
            β == γ && continue
            vβ = SVector(v[1, β], v[2, β])
            vγ == vβ && _throw_equal_velocities(γ, β)
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
node, as for no particle; [`apply_G!`](@ref) throws for two equal velocities. `v` or `v̇` of
another number of columns than the model has particles throws a `DimensionMismatch`.
"""
function collisional_vectorfield!(v̇::AbstractArray{ST}, v::AbstractArray{ST}, params,
        model::PairwiseLandau) where {ST}
    ∇S = model.landau.cache[ST].v̇
    entropy_gradient!(∇S, v, model)
    apply_G!(v̇, v, ∇S, model)
    v̇ .*= ST(model.landau.ν)
    return nothing
end
