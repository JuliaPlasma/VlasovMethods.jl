"""
    DiscreteGradient

The discrete gradients ``\\bar{\\nabla} S(v^0, v^1)`` of the entropy that
[`discrete_gradient_step!`](@ref) takes, each a method of `discrete_gradient!`. A discrete
gradient is defined by the two conditions of the Landau manuscript's `eq:disc_grad_props`:
``S(v^1) - S(v^0) = (v^1 - v^0) \\cdot \\bar{\\nabla} S(v^0, v^1)`` and
``\\bar{\\nabla} S(v, v) = \\nabla S(v)``. The [`GonzalezGradient`](@ref) is the one this package
implements.
"""
abstract type DiscreteGradient end

@doc raw"""
    GonzalezGradient()

The Gonzalez discrete gradient, the Landau manuscript's choice in `eq:disc_grad_props`:

```math
\bar{\nabla} S(v^0, v^1) = \nabla S(\bar{v})
  + \frac{S(v^1) - S(v^0) - \Delta v \cdot \nabla S(\bar{v})}{|\Delta v|^2} \, \Delta v ,
\qquad \bar{v} = \frac{v^0 + v^1}{2} , \quad \Delta v = v^1 - v^0 ,
```

with the dot product and the norm over all `2N` entries. It satisfies
``S(v^1) - S(v^0) = \Delta v \cdot \bar{\nabla} S`` to rounding. Where the computed
``|\Delta v|^2`` is zero it is ``\nabla S(\bar{v})``, without a division.
"""
struct GonzalezGradient <: DiscreteGradient end

"""
    discrete_gradient!(∇̄S, v0, v1, model::PairwiseLandau, ::GonzalezGradient)

The [`GonzalezGradient`](@ref) of the entropy `S_h` of `model` between the velocities `v0` and
`v1` (`2 × N`), written into `∇̄S`. At `v1 == v0` it is `entropy_gradient!` at `v0`, bit for
bit. A call projects the particles at the midpoint for `∇S_h`, and, where `v1 ≠ v0`, at `v1`
and at `v0` for `S_h`. The midpoint goes into the buffer `v` of `model.landau.cache[T]`. At
`v1 == v0` a call allocates nothing; otherwise the two evaluations of `S_h` by
[`CollisionEntropy`](@ref) allocate.
"""
function discrete_gradient!(∇̄S, v0::AbstractMatrix{T}, v1::AbstractMatrix{T},
        model::PairwiseLandau, ::GonzalezGradient) where {T}
    v̄ = model.landau.cache[T].v
    v̄ .= (v0 .+ v1) ./ 2
    entropy_gradient!(∇̄S, v̄, model)

    n2 = zero(T)
    Δv∇S = zero(T)
    for i in eachindex(v0, v1, ∇̄S)
        Δ = v1[i] - v0[i]
        n2 += Δ * Δ
        Δv∇S += Δ * ∇̄S[i]
    end
    iszero(n2) && return ∇̄S

    c = (_entropy(v1, model) - _entropy(v0, model) - Δv∇S) / n2
    ∇̄S .+= c .* (v1 .- v0)
    return ∇̄S
end

@doc raw"""
    discrete_gradient_step!(model::PairwiseLandau, Δt; abstol, max_iters = 100,
        gradient = GonzalezGradient())

Take one step of the discrete gradient scheme `eq:disc_grad_update_rule` of the Landau
manuscript,

```math
\frac{v^{n+1} - v^n}{\Delta t} = \nu \, \mathbb{G}\!\left( \frac{v^n + v^{n+1}}{2} \right)
  \bar{\nabla} S_h(v^n, v^{n+1}) ,
```

with ``\mathbb{G}`` at the midpoint (`eq:midpoint_G`) and the discrete gradient `gradient`, from
the velocities ``v^n`` of `model.landau.dist.particles.v`; write the accepted velocities into the
particles and return `(iterations = …, residual = …)`.

From ``v_0 = v^n``, iteration ``k = 0, 1, …`` evaluates
``f_k = \nu \, \mathbb{G}(m_k) \, \bar{\nabla} S_h(v^n, v_k)`` with ``m_k = (v^n + v_k)/2``, and
the residual ``r_k = v_k - v^n - \Delta t \, f_k``. The first iterate with
``\|r_k\|_\infty \le`` `abstol` is accepted: the particles take ``v^n + \Delta t \, f_k``, and
the step returns `iterations = k + 1`, the number of evaluations of ``f``, and
`residual` ``= \|r_k\|_\infty``. Otherwise ``v_{k+1} = v^n + \Delta t \, f_k``. The first
evaluation, at ``v_0 = v^n``, takes ``\bar{\nabla} S_h(v^n, v^n) = \nabla S_h(v^n)``, the
explicit Euler step.

# Structure

The accepted state conserves the momentum ``\sum_\alpha w_\alpha v_\alpha`` to rounding, because
``\sum_\gamma w_\gamma \mathbb{G}_{\gamma\sigma} = 0``. With ``\mathbb{G}`` and the iterate at the
same midpoint ``m_k``, ``\sum_\alpha w_\alpha f_\alpha \cdot m_\alpha = 0``, so the energy
changes by ``-(\Delta t / 2) \sum_\alpha w_\alpha f_\alpha \cdot r_\alpha``, at most
``(\Delta t / 2) \, \text{abstol} \sum_\alpha w_\alpha \|f_\alpha\|_1``. The discrete gradient
gives ``S_h(v^{n+1}) - S_h(v^n) = (v^{n+1} - v^n) \cdot \bar{\nabla} S_h`` (`eq:disc_grad_props`),
and so the discrete H-theorem: ``S_h`` does not increase, because ``\mathbb{G}`` is negative
semi-definite and ``\Delta t > 0``. The H-theorem holds for the fixed point of the iteration;
the accepted state differs from it by a residual of at most `abstol`.

# Tolerance

The residual stalls at the rounding of one evaluation, and `abstol` must lie above that floor,
or the step runs to `max_iters` and throws. On the fixture of `test/helpers/pairwise_landau.jl`
(144 particles on ``[-4, 4]^2``, `γ = -3`), the residuals of the 20th to the 30th evaluation lie
between `1.1e-15` and `1.6e-14` in `Float64` and between `3.6e-7` and `1.5e-5` in `Float32`, for
`Δt` from `5e-5` to `8e-3`; for `γ = 0` at `Δt = 5e-5`, between `1.1e-14` and `5.3e-14` in
`Float64` and up to `8.1e-5` in `Float32`. The iteration contracts by a factor of about 100 per
evaluation at `Δt = 2e-3` and `γ = -3`, and more slowly for `γ = 0`. The manuscript's stall "at
about ``10^{-11}``" is the same floor in the relative increment
``\|v_{k+1} - v_k\| / \|v_{k+1} - v^n\|``.

# Errors

`Δt`, converted to the element type of the velocities, and `abstol` must be finite and positive
and `max_iters` at least 1, or an `ArgumentError` that names the argument is thrown before any
evaluation: a `Float64` `Δt` that is zero or `Inf` in `Float32` is refused. So is a particle
weight that is not positive and finite, with an `ArgumentError` that names the particle. When
no iterate meets `abstol` within `max_iters` evaluations, an `ErrorException` is thrown that
names the number of evaluations and the last residual. A residual that is not finite never
meets `abstol`. The exceptions of the vector field pass through: a `DomainError` for two equal
velocities and, for ``\gamma < -2``, for two velocities too close, and the projection's
`DomainError` once an iterate leaves the velocity domain. In each case the particles are not
written.

# Cost

One evaluation of ``f`` is one ``\nabla S_h``, two ``S_h`` and one [`apply_G!`](@ref) of
``O(N^2)``; the first evaluation takes no ``S_h``. A step allocates four `2 × N` arrays, and
each ``S_h`` of [`CollisionEntropy`](@ref) allocates: 92 kB on 8 × 8 quadratic cells.
"""
function discrete_gradient_step!(model::PairwiseLandau, Δt; abstol,
        max_iters::Integer = 100, gradient::DiscreteGradient = GonzalezGradient())
    vn = model.landau.dist.particles.v
    T = eltype(vn)
    h = T(Δt)
    isfinite(h) && h > 0 ||
        throw(ArgumentError("Δt must be finite and positive in $T, got $Δt"))
    isfinite(abstol) && abstol > 0 ||
        throw(ArgumentError("abstol must be finite and positive, got $abstol"))
    max_iters ≥ 1 || throw(ArgumentError("max_iters must be at least 1, got $max_iters"))
    _check_weights(view(model.landau.dist.particles.w, 1, :))

    ν = T(model.landau.ν)

    vk = copy(vn)
    m = similar(vk)
    ∇̄S = similar(vk)
    f = similar(vk)

    residual = zero(T)
    for k in 0:(max_iters - 1)
        discrete_gradient!(∇̄S, vn, vk, model, gradient)
        m .= (vn .+ vk) ./ 2
        apply_G!(f, m, ∇̄S, model)
        f .*= ν

        residual = zero(T)
        for i in eachindex(vk, vn, f)
            residual = max(residual, abs(vk[i] - (vn[i] + h * f[i])))
        end

        if residual ≤ abstol
            vn .+= h .* f
            return (iterations = k + 1, residual = residual)
        end

        vk .= vn .+ h .* f
    end

    error("the discrete gradient iteration did not converge in $max_iters evaluations: ",
        "the last residual is $residual, above abstol = $abstol")
end
