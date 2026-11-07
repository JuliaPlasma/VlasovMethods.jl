"""
    DiscreteGradient

The discrete gradients ``\\bar{\\nabla} S(v^0, v^1)`` of the entropy that
[`discrete_gradient_step!`](@ref) takes, each a method of `discrete_gradient!`. A discrete
gradient is defined by the two conditions of the Landau manuscript's `eq:disc_grad_props`:
``S(v^1) - S(v^0) = (v^1 - v^0) \\cdot \\bar{\\nabla} S(v^0, v^1)`` and
``\\bar{\\nabla} S(v, v) = \\nabla S(v)``. The package implements two of them: the
[`GonzalezGradient`](@ref), the manuscript's, and the [`AVFGradient`](@ref).
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

# The `n`-point Gauß–Legendre rule on `[0, 1]`: the nodes in increasing order and the weights,
# in `Float64`, by Newton's method on the Legendre polynomial, which converges in a few
# iterations from the asymptotic position of each root.
function _gauss_legendre_rule(n::Int)
    ξ = zeros(Float64, n)
    ω = zeros(Float64, n)
    for i in 1:((n + 1) ÷ 2)
        z = cos(Float64(π) * (i - 0.25) / (n + 0.5))
        dP = one(Float64)
        for _ in 1:100
            P1, P2 = one(Float64), zero(Float64)
            for j in 1:n
                P3 = P2
                P2 = P1
                P1 = ((2j - 1) * z * P2 - (j - 1) * P3) / j
            end
            dP = n * (z * P1 - P2) / (z * z - 1)
            z1 = z
            z = z1 - P1 / dP
            # The last Newton step is the rounding of the root, so the stop is a few ulp of `z`;
            # a fixed 1e-16 lies below it in [0.5, 1) and fires only where `z` repeats exactly.
            abs(z - z1) ≤ 4 * eps(z) && break
        end
        ξ[i] = -z
        ξ[n + 1 - i] = z
        ω[i] = 2 / ((1 - z * z) * dP * dP)
        ω[n + 1 - i] = ω[i]
    end
    return (ξ .+ 1) ./ 2, ω ./ 2
end

@doc raw"""
    AVFGradient(nq::Integer = 5)

The averaged-vector-field discrete gradient: the mean of the entropy gradient along the segment
between the two velocities,

```math
\bar{\nabla} S_h(v^0, v^1) = \int_0^1 \nabla S_h\big( (1 - \tau) \, v^0 + \tau \, v^1 \big) \,
  d\tau ,
```

taken by the `nq`-point Gauß–Legendre rule on ``[0, 1]``,

```math
\bar{\nabla} S_h(v^0, v^1) = \sum_{i = 1}^{nq} \omega_i \,
  \nabla S_h\big( v^0 + \tau_i \, \Delta v \big) , \qquad \Delta v = v^1 - v^0 ,
```

with the nodes ``\tau_i`` in increasing order and symmetric about ``1/2`` and the weights
``\omega_i > 0`` summing to one. The rule integrates ``\tau^k`` exactly for ``k \le 2nq - 1``,
because Gauß–Legendre of `nq` nodes is exact for polynomials of degree ``2nq - 1``; the
constructor computes it in `Float64` by Newton's method on the Legendre polynomial, once per
rule, and converts a node and a weight to the element type of the velocities at each use.

The two conditions of the Landau manuscript's `eq:disc_grad_props` hold to the quadrature error
and not to rounding, because ``S_h`` is not a polynomial in the velocities:
``S_h(v^1) - S_h(v^0) = \Delta v \cdot \bar{\nabla} S_h`` holds to the error of the rule, and
``\bar{\nabla} S_h(v, v) = \nabla S_h(v)`` to the rounding of ``\sum_i \omega_i``. A test of the
identity states its tolerance above the quadrature error, which the rounding does not bound. At
`nq = 1` the single node is ``1/2``, so the rule is the midpoint gradient
``\nabla S_h \big( (v^0 + v^1)/2 \big)``.

# Cost

A call takes `nq` evaluations of ``\nabla S_h`` and no ``S_h``, where the
[`GonzalezGradient`](@ref) takes one ``\nabla S_h`` and two ``S_h``; one
[`apply_G!`](@ref) of `O(N²)` follows either in the step. On the fixture of
`test/helpers/pairwise_landau.jl` (144 particles, quadratic splines on 8 × 8 cells,
`Float64`, at a load of about 7), ``\nabla S_h`` takes 54 µs and ``S_h`` 34 µs, so one call
of the [`GonzalezGradient`](@ref) takes 121 µs and one of `AVFGradient(5)` 271 µs, 2.24 times
as long. A call allocates nothing after a warm call: the node velocities go into the buffer `v`
of `model.landau.cache[T]`, the accumulated gradient into `∇̄S`, and the gradient of each later
node into the buffer `v̇`.

# Errors

An `nq` below 1 throws an `ArgumentError` that names it, before any node is computed; `nq` that is
not an integer has no method.
"""
struct AVFGradient <: DiscreteGradient
    nq::Int
    nodes::Vector{Float64}
    weights::Vector{Float64}

    function AVFGradient(nq::Integer = 5)
        nq ≥ 1 || throw(ArgumentError("nq must be at least 1, got $nq"))
        nodes, weights = _gauss_legendre_rule(Int(nq))
        return new(Int(nq), nodes, weights)
    end
end

@doc raw"""
    discrete_gradient!(∇̄S, v0::AbstractMatrix{T}, v1::AbstractMatrix{T}, model::PairwiseLandau, ::AVFGradient) where {T}

The [`AVFGradient`](@ref) of the entropy `S_h` of `model` between the velocities `v0` and `v1`
(`2 × N`), written into `∇̄S`: the `nq`-point Gauß–Legendre rule of `∇S_h` along the segment, with
the node and the weight converted to the element type of `v0`. A call evaluates no `S_h`. The
node velocities go into the buffer `v` of `model.landau.cache[T]` and the gradient of each node
after the first into the buffer `v̇`, so a call allocates nothing after a warm call.
"""
function discrete_gradient!(∇̄S, v0::AbstractMatrix{T}, v1::AbstractMatrix{T},
        model::PairwiseLandau, rule::AVFGradient) where {T}
    τ = rule.nodes
    ω = rule.weights
    cache = model.landau.cache[T]
    vnode = cache.v
    tmp = cache.v̇
    vnode .= v0 .+ T(τ[1]) .* (v1 .- v0)
    entropy_gradient!(∇̄S, vnode, model)
    ∇̄S .*= T(ω[1])
    for i in 2:length(τ)
        vnode .= v0 .+ T(τ[i]) .* (v1 .- v0)
        entropy_gradient!(tmp, vnode, model)
        ∇̄S .+= T(ω[i]) .* tmp
    end
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

# The discrete gradients

The keyword `gradient` takes either rule. [`GonzalezGradient`](@ref), the default, satisfies the
identity ``S_h(v^{n+1}) - S_h(v^n) = (v^{n+1} - v^n) \cdot \bar{\nabla} S_h`` of
`eq:disc_grad_props` to rounding, because it corrects the midpoint gradient by the two entropies
and the displacement. [`AVFGradient`](@ref) averages the entropy gradient along the segment with
an `nq`-point Gauß–Legendre rule, and takes that identity to the quadrature error and not to
rounding, because ``S_h`` is not a polynomial in the velocities: the rounding of its sum does not
bound that error, so a test of the identity with this rule states its tolerance above the
quadrature error. The conservation above holds for both rules, because it comes from
``\mathbb{G}`` at the midpoint and the weights, not from the gradient. The H-theorem holds for
[`GonzalezGradient`](@ref) to rounding, and for [`AVFGradient`](@ref) to the quadrature error of
its identity, which is not zero: the entropy falls where
``\Delta v \cdot \bar{\nabla} S_h`` exceeds that error. The two rules differ in cost, in the
section below.

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
`DomainError` for a velocity outside the domain of the spline basis, which the step evaluates at
the accepted state before the particles take it. A rule's evaluations lie at its own points, and
the state they accept can be outside the domain where none of those points is. In each case the
particles are not written.

# Cost

One evaluation of ``f`` is one [`apply_G!`](@ref) of ``O(N^2)`` and the gradient it multiplies:
one ``\nabla S_h`` and, from the second evaluation, two ``S_h`` for the
[`GonzalezGradient`](@ref); `nq` ``\nabla S_h`` and no ``S_h`` for the [`AVFGradient`](@ref). A
step allocates four `2 × N` arrays and takes one projection of the accepted state, for either
rule, and each ``S_h`` of [`CollisionEntropy`](@ref) allocates: 92 kB on 8 × 8 quadratic cells.
"""
function discrete_gradient_step!(model::PairwiseLandau, Δt; abstol,
        max_iters::Integer = 100, gradient::DiscreteGradient = GonzalezGradient())
    vn = model.landau.dist.particles.v
    T = eltype(vn)
    landau = model.landau
    cache = landau.cache[T]
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
            # The accepted state is projected before the particles take it: a rule evaluates at
            # its own points, and the state it accepts can lie outside the domain where none of
            # them does.
            m .= vn .+ h .* f
            _projection!(cache.rhs, cache.vals, cache.mwork, m, landau.dist, cache.sdist)
            vn .= m
            return (iterations = k + 1, residual = residual)
        end

        vk .= vn .+ h .* f
    end

    error("the discrete gradient iteration did not converge in $max_iters evaluations: ",
        "the last residual is $residual, above abstol = $abstol")
end
