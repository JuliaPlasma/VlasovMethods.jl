@doc raw"""
    LowRankParityDistribution{T, ST}

A 1D1V distribution function on a velocity grid symmetric about 0 in the parity form

```math
f = 1 \otimes f_0 + f^+ + f^- , \qquad
f^+ = X^+ S^+ (V^+)^\mathsf{T} , \qquad f^- = X^- S^- (V^-)^\mathsf{T} ,
```

a `DistributionFunction{T, 1, 1}`, the state of [`LowRankParitySplitting`](@ref). Its fields:

  - `f0`, the even background, a `Vector{T}` over the `nv` velocity nodes;
  - `plus`, the [`LowRankDistribution`](@ref) of the even part ``f^+``: every column of `V` is
    even, and its first two columns, `nfixed = 2`, span ``\{1, v^2\}``;
  - `minus`, the [`LowRankDistribution`](@ref) of the odd part ``f^-``: every column of `V` is
    odd, and its first column, `nfixed = 1`, is ``\pm v / \|v\|``.

The parity of a column holds bit for bit: `V[j, :] == V[nv + 1 - j, :]` in `plus` and
`V[j, :] == -V[nv + 1 - j, :]` in `minus`. `size(Y)` is `(nx, nv)` and `eltype(Y)` is `T`.

    LowRankParityDistribution(dist::GridDistribution, rank::Tuple{Integer, Integer};
                              background = nothing)

Factorise the node values `f = dist.values` on [`GridFactorSpace`](@ref)`(dist)` at the ranks
`rank = (r⁺, r⁻)`. `background = nothing` takes the x-mean of the even part of `f`,
``f_{0,j} = \frac{1}{n_x} \sum_i (f_{ij} + f_{i, n_v + 1 - j}) / 2``; a given vector is
replaced by ``(f_0 + \mathrm{reverse}(f_0)) / 2``. With
`(f⁺, f⁻) = `[`parity_split`](@ref)`(f, f0)` and ``h_x``, ``h_v`` the steps of the grid:

 1. the first columns of `plus.V` are the weighted orthonormal basis ``V_m`` of
    ``\{1, v^2\}``, then follow the leading right singular vectors of
    ``\sqrt{h_x h_v} (f^+ - h_v f^+ V_m V_m^\mathsf{T})``, `r⁺ - 2` of them, so that the
    vectors of vanishing singular values complete the even basis to `r⁺` columns;
 2. `plus.X, plus.S = orthonormalize_x(space, h_v f⁺ V)`, so `plus.S` is upper triangular;
 3. `minus` is built the same way from ``f^-``, with the fixed mode ``v`` and `r⁻` columns.

Each basis comes from the weighted Householder QR on the half grid of its parity: the rows
``j \le \lceil n_v / 2 \rceil`` with the weight ``2 h_v``, ``h_v`` on the middle row of an odd
`nv`, for an even basis, and the rows ``j \le \lfloor n_v / 2 \rfloor`` with the weight
``2 h_v`` for an odd basis, mirrored onto the other half. So the completion columns have their
parity too. With full ranks, `(min(nx, cld(nv, 2)), min(nx, fld(nv, 2)))` for `nx ≥ cld(nv, 2)`,
the factors reproduce `f` to round-off; at lower ranks they keep the density, the momentum and
the kinetic energy of `f` at every x node.

Throws an `ArgumentError` for a v-grid that is not symmetric,
``|v_j + v_{n_v + 1 - j}| > 10 \, \mathrm{eps}(T) \max(|v_1|, |v_{n_v}|)`` for some ``j``; for a
background that is a matrix or neither `nothing` nor a vector, that is not even
(``10 \, \mathrm{eps}(T) \max |f_0|``), or that holds a `NaN` or an `Inf`; for values that hold
a `NaN` or an `Inf`; and for `r⁺` outside `2:min(nx, cld(nv, 2))` or `r⁻` outside
`1:min(nx, fld(nv, 2))`. Throws a `DimensionMismatch` for a background of another length than `nv`.
"""
struct LowRankParityDistribution{T, ST <: AbstractFactorSpace{T}} <:
       DistributionFunction{T, 1, 1}
    f0::Vector{T}
    plus::LowRankDistribution{T, ST}
    minus::LowRankDistribution{T, ST}
end

function LowRankParityDistribution(
        dist::GridDistribution{T}, rank::Tuple{Integer, Integer};
        background = nothing) where {T}
    space = GridFactorSpace(dist)
    nx, nv = size(space)
    values = dist.values
    f0 = _parity_checks(values, space.v, background, "the low-rank parity distribution",
        "the values")
    _check_parity_ranks(rank, nx, nv)

    fp, fm = parity_split(values, f0)
    plus = _parity_factors(space, fp, +1, rank[1])
    minus = _parity_factors(space, fm, -1, rank[2])
    return LowRankParityDistribution(f0, plus, minus)
end

# the ranks (r⁺, r⁻) of the two parts, at most the dimensions of the x-space and of the even and
# the odd subspace of the v-grid
function _check_parity_ranks(rank, nx, nv)
    r⁺max, r⁻max = min(nx, cld(nv, 2)), min(nx, fld(nv, 2))
    2 ≤ rank[1] ≤ r⁺max || throw(ArgumentError(
        "the rank $(rank[1]) of f⁺ is not in 2:$r⁺max on a grid of $nx × $nv nodes"))
    1 ≤ rank[2] ≤ r⁻max || throw(ArgumentError(
        "the rank $(rank[2]) of f⁻ is not in 1:$r⁻max on a grid of $nx × $nv nodes"))
    return nothing
end

Base.eltype(::LowRankParityDistribution{T}) where {T} = T
Base.size(Y::LowRankParityDistribution) = size(Y.plus)

"""
    GridDistribution(Y::LowRankParityDistribution)

The [`GridDistribution`](@ref) of the values
[`parity_merge`](@ref)`(Y.f0, X⁺ S⁺ V⁺ᵀ, X⁻ S⁻ V⁻ᵀ)` on the grid of the [`GridFactorSpace`](@ref)
of `Y`. It forms the `nx × nv` array, so it serves tests and diagnostics.
"""
function GridDistribution(Y::LowRankParityDistribution{T, <:GridFactorSpace}) where {T}
    p, m = Y.plus, Y.minus
    space = p.space
    return GridDistribution(space.x, space.v,
        parity_merge(Y.f0, p.X * p.S * p.V', m.X * m.S * m.V'))
end

# An orthonormal basis [Vfix Vnew] of span [Vfix, A] in the subspace of the given parity, +1 even
# and -1 odd, of at most the dimension of that subspace: the parity part of A, weighted on the
# half grid so that the sum over the half grid is the rectangle rule of the whole grid, by the
# Householder QR, mirrored onto the other half. `Vfix` is orthonormal, of that parity, and kept
# as it is in the leading columns.
function _parity_basis(space::GridFactorSpace, Vfix::AbstractMatrix{T},
        A::AbstractMatrix{T}, parity) where {T}
    nv, hv = length(space.v), step(space.v)
    nh = parity > 0 ? cld(nv, 2) : fld(nv, 2)
    m, k = size(Vfix, 2), size(Vfix, 2) + size(A, 2)
    s = [sqrt(hv * (2j == nv + 1 ? one(T) : T(2))) for j in 1:nh]
    C = Matrix{T}(undef, nh, k)
    for j in 1:nh
        for c in 1:m
            C[j, c] = s[j] * Vfix[j, c]
        end
        for c in axes(A, 2)
            C[j, m + c] = s[j] * (A[j, c] + parity * A[nv + 1 - j, c]) / 2
        end
    end
    n = min(nh, k)
    V = zeros(T, nv, n)
    V[:, 1:m] .= Vfix
    Q = Matrix(qr(C).Q)
    for c in (m + 1):n, j in 1:nh

        V[j, c] = Q[j, c] / s[j]
    end
    return _mirror!(V, parity)
end

# Rows nv + 1 - j of V set to the parity times rows j ≤ nv/2. A product `V Q` of a basis of one
# parity has that parity in exact arithmetic, and its rows j and nv + 1 - j round alike only where
# the BLAS sums them alike; the copy makes the parity exact. The middle row of an odd basis on an
# odd nv is zero in every input, and so in every product.
function _mirror!(V, parity)
    nv = size(V, 1)
    for j in 1:fld(nv, 2)
        V[nv + 1 - j, :] .= parity .* view(V, j, :)
    end
    return V
end

# the weighted orthonormal fixed modes: {1, v²} for the even part, v for the odd part
function _parity_fixed_modes(space::GridFactorSpace{T}, parity) where {T}
    v = collect(space.v)
    A = parity > 0 ? hcat(ones(T, length(v)), v .^ 2) : reshape(v, :, 1)
    return _parity_basis(space, zeros(T, length(v), 0), A, parity)
end

# the factors of the part f of the given parity at the rank r: the fixed modes, then the basis that
# the half-grid QR makes of the r - m leading right singular vectors of the part of f outside
# them, which completes the parity subspace where the singular values vanish
function _parity_factors(space::GridFactorSpace, f, parity, r)
    hx, hv = step(space.x), step(space.v)
    Vm = _parity_fixed_modes(space, parity)
    m = size(Vm, 2)
    g = f - (hv .* (f * Vm)) * Vm'
    V = _parity_basis(space, Vm, svd(sqrt(hx * hv) .* g).V[:, 1:(r - m)] ./ sqrt(hv), parity)
    X, S = orthonormalize_x(space, hv .* (f * V))
    return LowRankDistribution(space, X, S, V; nfixed = m)
end

@doc raw"""
    _truncate_keeping_moments(space, K, V, nfixed; tol, minrank, rank)

The truncation of ``f = K V^\mathsf{T}``, with `V` orthonormal in the weighted product of
`space` and its first `m = nfixed` columns the fixed modes, that keeps the projection
``K_{:, 1:m} V_{:, 1:m}^\mathsf{T}`` of ``f`` onto ``L^2_x \otimes \mathrm{span}\, V_{:, 1:m}``
exactly. The rest ``K_\perp V_\perp^\mathsf{T}`` is truncated by the weighted SVD
``\sqrt{h_x} K_\perp = U \Sigma Q^\mathsf{T}`` to ``r_1 - m`` columns, with
`r₁ = clamp(r*, minrank, rank)`, capped at `m` plus the number of singular values, and `r*` the
least `r ≥ m` with ``(\sum_{j > r - m} \sigma_j^2)^{1/2} \le`` `tol` ``\|f\|``,
``\|f\| = \sqrt{h_x} \|K\|_F``: `tol` is relative to the norm of `f`. The whole tail is at most
``\|f\|``, so a `tol ≥ 1`, `Inf` included, gives `r₁ = minrank` with the same cap.

Returns `(X, S, V′)`: `X, S = orthonormalize_x(space, [K[:, 1:m]  K⊥ Q₁])`, so `S` is upper
triangular, and `V′ = [V[:, 1:m]  V⊥ Q₁]`. So ``X S V'^\mathsf{T} h_v V_{:, 1:m} = K_{:, 1:m}``
to round-off. Returns `nothing` when ``\|f\|`` is not finite: when `K` holds a `NaN` or an `Inf`,
or when its entries are finite and its norm overflows.
"""
function _truncate_keeping_moments(space::GridFactorSpace, K::AbstractMatrix{T},
        V::AbstractMatrix{T}, nfixed; tol, minrank, rank) where {T}
    hx = step(space.x)
    ν = sqrt(hx) * norm(K)
    isfinite(ν) || return nothing
    m = nfixed
    K⊥ = K[:, (m + 1):end]
    F = svd(sqrt(hx) .* K⊥)
    n = length(F.S)
    lo, hi = min(minrank - m, n), min(rank - m, n)
    # the whole tail is at most ‖f‖, so a tol ≥ 1 keeps the least rank; tol ν can be Inf ⋅ 0
    r = tol ≥ 1 ? lo : clamp(_least_rank(F.S, tol * ν), lo, hi)
    Q₁ = F.V[:, 1:r]
    X, S = orthonormalize_x(space, hcat(K[:, 1:m], K⊥ * Q₁))
    return X, S, hcat(V[:, 1:m], V[:, (m + 1):end] * Q₁)
end

"""
    _closed_space(space, X̂, G)

The x-basis, orthonormal in the weighted product of `space`, of the span of `[X̂ 1 G X̂]`, with
`1` the constant: the spatial trial space of the S-step of [`LowRankParitySplitting`](@ref). It
has `size(X̂, 2) + 1 + size(X̂, 2)` columns, or `nx`, the whole x-space, when that number is `nx`
or more.
"""
function _closed_space(space::GridFactorSpace{T}, X̂, G) where {T}
    return _basis(hcat(X̂, ones(T, length(space.x)), G * X̂), step(space.x))
end

# the implicit midpoint step y₁ = y₀ + τ R((y₀ + y₁)/2) by two fixed-point iterations from
# y₁ = y₀. R depends on y only through the density of y, and the density rate does not depend on
# the potential, so the derivative J of the map has J² = 0 and two iterations solve it exactly.
function _midpoint(R, y₀, τ)
    y₁ = y₀ .+ τ .* R(y₀)
    return y₀ .+ τ .* R((y₀ .+ y₁) ./ 2)
end

# the periodic shifts S_t V and (S_q k) ∘ (S_t V), q, t ∈ {-1, 0, 1}: the columns whose span holds
# the rows of the bracket term P(X Vᵀ, ϕ ⊗ 1 + 1 ⊗ k) for any X and ϕ
function _velocity_range(kinetic, V)
    return reduce(hcat,
        [q == 2 ? circshift(V, (t, 0)) : circshift(kinetic, q) .* circshift(V, (t, 0))
         for t in -1:1 for q in -1:2])
end

# the potential of 1 ⊗ f0 + K Vᵀ
function _parity_potential(field::LowRankFields1D1V{T}, f0, K, V) where {T}
    nx = size(K, 1)
    return _potential!(Vector{T}(undef, nx), hcat(ones(T, nx), K), hcat(f0, V), field)
end

# The f⁻ substep of length τ: f⁻ + τ F(g), g = 1 ⊗ f0 + f⁺, on V, the odd orthonormal basis of
# [V⁻ R] with R the velocity range of [f0 V⁺], the fixed mode first, then the truncation that
# keeps the fixed mode. The rows of F(g) are odd and lie in span R, so F(g) = K_inc Vᵀ with
# K_inc = F(g) hv V, and the sum is exact on V. When the norm of the truncated K is not finite,
# the truncation returns nothing and the substep leaves a NaN core.
function _parity_fminus_substep!(Y::LowRankParityDistribution{T}, method, field, τ) where {T}
    plus, minus = Y.plus, Y.minus
    space = minus.space
    nx = length(space.x)
    hv = step(space.v)
    r⁺ = size(plus.S, 1)
    X⁻, S⁻, V⁻ = minus.X, minus.S, minus.V
    Vg = hcat(Y.f0, plus.V)
    V = _parity_basis(space, V⁻[:, 1:1],
        hcat(V⁻[:, 2:end], _velocity_range(field.kinetic, Vg)), -1)
    # the last columns of the K-field on [Kg 0] and [Vg V] are K_inc = F(g) hv V
    Kg = hcat(ones(T, nx), plus.X * plus.S, zeros(T, nx, size(V, 2)))
    K̇ = _kfield!(similar(Kg), Kg, hcat(Vg, V), field)
    Kinc = K̇[:, (r⁺ + 2):end]

    K = X⁻ * (S⁻ * (hv .* (V⁻' * V))) .+ τ .* Kinc
    t = _truncate_keeping_moments(space, K, V, 1;
        tol = T(method.tol[2]), minrank = method.minrank[2], rank = method.rank[2])
    if t === nothing
        minus.S = fill(T(NaN), size(S⁻))
        return Y
    end
    X, S, V′ = t
    minus.X, minus.S, minus.V, minus.nfixed = X, S, _mirror!(V′, -1), 1
    return Y
end

# The f⁺ step of length τ with f⁻ frozen: the augmented BUG step of ḟ⁺ = F_ϕ(f⁻), ϕ the potential
# of 1 ⊗ f0 + f⁺, with the S-step on the closed space and every substep an implicit midpoint step,
# then the truncation that keeps the fixed modes. The cross terms of f⁻ and the bases of f⁺ are
# the blocks of the fields on stacked factors. When the norm of the truncated Xa S₁ is not finite,
# the truncation returns nothing and the step leaves a NaN core.
function _parity_fplus_step!(Y::LowRankParityDistribution{T}, method, field, τ) where {T}
    plus, minus = Y.plus, Y.minus
    space = plus.space
    nx, nv = size(space)
    hx, hv = step(space.x), step(space.v)
    f0 = Y.f0
    X₀, S₀, V₀ = plus.X, plus.S, plus.V
    X⁻, S⁻, V⁻ = minus.X, minus.S, minus.V
    r₀, r⁻ = size(S₀, 1), size(S⁻, 1)
    tail(r) = (r⁻ + 1):(r⁻ + r)

    # the K-step from X₀ S₀: the last r₀ columns of the K-field on [X⁻S⁻ 0] and [V⁻ V₀]
    Kc, Vc = hcat(X⁻ * S⁻, zeros(T, nx, r₀)), hcat(V⁻, V₀)
    K̇ = similar(Kc)
    kstep(K) = _kfield!(K̇, Kc, Vc, field, _parity_potential(field, f0, K, V₀))[:, tail(r₀)]
    K₁ = _midpoint(kstep, X₀ * S₀, τ)

    # the L-step from V₀ S₀ᵀ: the last r₀ columns of the L-field on [X⁻ X₀] and [V⁻S⁻ᵀ 0]
    Xc, Lc = hcat(X⁻, X₀), hcat(V⁻ * S⁻', zeros(T, nv, r₀))
    L̇ = similar(Lc)
    lstep(L) = _lfield!(L̇, Xc, Lc, field, _parity_potential(field, f0, X₀, L))[:, tail(r₀)]
    L₁ = _midpoint(lstep, V₀ * S₀', τ)

    X̂ = _basis(hcat(K₁, X₀), hx)
    V̂ = _parity_basis(space, V₀[:, 1:2], hcat(V₀[:, 3:end], L₁), +1)
    Xa = _closed_space(space, X̂, field.G)

    # the S-step on Xa and V̂: a block of the S-field on [X⁻ Xa 0], blockdiag(S⁻, 0) and
    # [V⁻ V̂ 0], padded with zero columns to one rank
    na, nb = size(Xa, 2), size(V̂, 2)
    n = max(na, nb)
    Xs, Vs = hcat(X⁻, _pad(Xa, n)), hcat(V⁻, _pad(V̂, n))
    Ss = zeros(T, r⁻ + n, r⁻ + n)
    Ss[1:r⁻, 1:r⁻] .= S⁻
    Ṡ = similar(Ss)
    function sstep(S)
        ϕ = _parity_potential(field, f0, Xa * S, V̂)
        return _sfield!(Ṡ, Xs, Ss, Vs, field, ϕ)[tail(na), tail(nb)]
    end
    Ŝ₀ = (hx .* (Xa' * X₀)) * S₀ * (hv .* (V₀' * V̂))
    S₁ = _midpoint(sstep, Ŝ₀, τ)
    t = _truncate_keeping_moments(space, Xa * S₁, V̂, 2;
        tol = T(method.tol[1]), minrank = method.minrank[1], rank = method.rank[1])
    if t === nothing
        plus.S = fill(T(NaN), size(S₀))
        return Y
    end
    X, S, V′ = t
    plus.X, plus.S, plus.V, plus.nfixed = X, S, _mirror!(V′, +1), 2
    return Y
end

@doc raw"""
    LowRankParitySplitting(rank::Tuple{Integer, Integer}; tol = (0, 0), minrank = rank)

The low-rank parity splitting of 1D1V Vlasov–Poisson, a method for [`lowrank_step!`](@ref) on a
[`LowRankParityDistribution`](@ref) ``f = 1 \otimes f_0 + f^+ + f^-``. `rank`, `minrank` and
`tol` hold the largest ranks, the least ranks and the relative truncation tolerances of
``f^+`` (first) and ``f^-`` (second). With the default `minrank = rank` the ranks stay fixed.

Throws an `ArgumentError` unless `rank[1] ≥ 2`, `rank[2] ≥ 1`, `2 ≤ minrank[1] ≤ rank[1]`,
`1 ≤ minrank[2] ≤ rank[2]`, and each `tol` is finite and `≥ 0`.
"""
struct LowRankParitySplitting{TT <: Real}
    rank::NTuple{2, Int}
    minrank::NTuple{2, Int}
    tol::NTuple{2, TT}

    function LowRankParitySplitting(rank::Tuple{Integer, Integer};
            tol::Tuple{Real, Real} = (0, 0), minrank::Tuple{Integer, Integer} = rank)
        rank[1] ≥ 2 && rank[2] ≥ 1 || throw(ArgumentError(
            "LowRankParitySplitting takes a rank ≥ 2 of f⁺ and ≥ 1 of f⁻; got rank = $rank"))
        2 ≤ minrank[1] ≤ rank[1] && 1 ≤ minrank[2] ≤ rank[2] || throw(ArgumentError(
            "LowRankParitySplitting takes 2 ≤ minrank[1] ≤ rank[1] and 1 ≤ minrank[2] ≤ " *
            "rank[2]; got minrank = $minrank for rank = $rank"))
        all(t -> isfinite(t) && t ≥ 0, tol) || throw(ArgumentError(
            "LowRankParitySplitting takes finite tolerances tol ≥ 0; got tol = $tol"))
        tt = promote(tol...)
        new{typeof(tt[1])}(Int.(rank), Int.(minrank), tt)
    end
end

@doc raw"""
    lowrank_step!(Y::LowRankParityDistribution, method::LowRankParitySplitting,
                  field::LowRankFields1D1V, t, τ)

Advance `Y` by one step of length `τ` of the low-rank parity splitting, in place, and return
`Y`. `field` is the [`LowRankFields1D1V`](@ref) of the grid of `Y` with the
[`PeriodicVelocityClosure`](@ref). With ``P(f, h)`` the bracket term of the grid model,
``F_\phi(a) = P(a, \phi \otimes 1 + 1 \otimes v^2/2)``, ``\phi(g) = G h_v \sum_j g_{\cdot j}``,
``F(g) = F_{\phi(g)}(g)`` and ``F_0 = 1 \otimes f_0``, the step is

 1. the ``f^-`` substep of ``\tau/2``: with ``g = F_0 + f^+`` and `V` the odd orthonormal
    basis of ``V^-`` and of the periodic shifts ``S_t V_g`` and ``(S_q k) \circ (S_t V_g)``,
    ``q, t \in \{-1, 0, 1\}``, of ``V_g = [f_0 \; V^+]``, the fixed mode first, whose span holds
    the rows of ``F(g)``: ``f^- \leftarrow f^- + (\tau/2) K_\mathrm{inc} V^\mathsf{T}`` with
    ``K_\mathrm{inc} = h_v F(g) V``, exact, then the truncation that keeps the fixed mode, with
    `tol[2]`, `minrank[2]` and `rank[2]`;
 2. the ``f^+`` step of ``\tau`` with ``f^-`` frozen, an augmented BUG step of
    ``\dot{f}^+ = F_{\phi(F_0 + f^+)}(f^-)``: the K-step
    ``\dot{K} = h_v F_{\phi(F_0 + K V_0^\mathsf{T})}(f^-) V_0`` and the L-step
    ``\dot{L} = h_x F_{\phi(F_0 + X_0 L^\mathsf{T})}(f^-)^\mathsf{T} X_0``; `X̂` the
    orthonormal basis of ``[K(\tau) \; X_0]`` and `V̂` the even one of ``[V_0 \; L(\tau)]``, the
    fixed modes first; the closed space `Xa` of `[X̂ 1 G X̂]` ([`_closed_space`](@ref)); the
    Galerkin S-step ``\dot{S} = h_x h_v X_a^\mathsf{T} F_{\phi(F_0 + X_a S \hat{V}^\mathsf{T})}(f^-)
    \hat{V}`` from ``(h_x X_a^\mathsf{T} X_0) S_0 (h_v V_0^\mathsf{T} \hat{V})``; and the
    truncation that keeps the fixed modes, with `tol[1]`, `minrank[1]` and `rank[1]`;
 3. the ``f^-`` substep of ``\tau/2`` again.

Each substep of 2 is the implicit midpoint step, solved exactly by two fixed-point iterations:
its right-hand side depends on the state only through the density, whose rate does not depend on
the potential. The truncation is [`_truncate_keeping_moments`](@ref): it keeps the projection of
``f^+`` onto ``L^2_x \otimes \mathrm{span}\{1, v^2\}`` and of ``f^-`` onto
``L^2_x \otimes \mathrm{span}\{v\}``, so it keeps the density, the kinetic energy and the
momentum at every x node, and its `S` is upper triangular. The step keeps `nfixed`, 2 and 1, and
the parity of every column of `V` bit for bit.

The energy ``H = h_x h_v \sum f v^2/2 + (h_x/2) \rho^\mathsf{T} G \rho`` depends on ``f^+``
alone. Step 1 does not change ``f^+``; the trial space of the S-step holds
``G \rho \otimes 1 + 1 \otimes v^2/2`` at every point of it, so its midpoint step keeps ``H``;
the truncation keeps ``\rho`` and ``\sum v^2 f^+``. So the step keeps ``H`` to round-off, with no
projection or correction, and it keeps the mass. It evaluates the fields on the factors and
does not reconstruct ``f``; the sizes of its temporaries depend on the augmented ranks, and at
the ranks `(cld(nv, 2), fld(nv, 2))` some of them are ``n_x \times n_v``.

`τ` is converted to the element type `T` of `Y`, and a `τ` that is not finite in `T` throws the
`ArgumentError` of the generic [`lowrank_step!`](@ref) before any write; `τ = 0` returns `Y`
unchanged, with no truncation, and `τ < 0` steps back. A `field` with another closure than the
periodic one, a rank of `method` above `min(nx, cld(nv, 2))` for ``f^+`` or
`min(nx, fld(nv, 2))` for ``f^-``, and a state whose rank of ``f^+`` or ``f^-`` is below
`minrank[1]` or `minrank[2]` throw an `ArgumentError` before any write. A step of `τ ≠ 0`
truncates a state above `rank` down. A substep whose core is not finite does not throw: it
skips its truncation and writes a core of `NaN` at the input rank of its part, with the input
`X` and `V`, so `all(isfinite, Y.plus.S) && all(isfinite, Y.minus.S)` is the status. The core
is the matrix that the truncation reads, and a core whose norm is not finite takes the same
path: one with a `NaN` or `Inf` entry, and one whose entries are finite and whose norm
overflows.
"""
function lowrank_step!(Y::LowRankParityDistribution{T}, method::LowRankParitySplitting,
        field::LowRankFields1D1V{T}, t, τ) where {T}
    τ = _finite_step(τ, T)
    field.closure isa PeriodicVelocityClosure || throw(ArgumentError(
        "LowRankParitySplitting runs on the periodic velocity closure; got a field with " *
        "$(typeof(field.closure))"))
    nx, nv = size(Y)
    _check_parity_ranks(method.rank, nx, nv)
    r⁺, r⁻ = size(Y.plus.S, 1), size(Y.minus.S, 1)
    r⁺ ≥ method.minrank[1] || throw(ArgumentError(
        "the rank $r⁺ of f⁺ is below minrank[1] = $(method.minrank[1]) of the method"))
    r⁻ ≥ method.minrank[2] || throw(ArgumentError(
        "the rank $r⁻ of f⁻ is below minrank[2] = $(method.minrank[2]) of the method"))
    iszero(τ) && return Y
    _parity_fminus_substep!(Y, method, field, τ / 2)
    _parity_fplus_step!(Y, method, field, τ)
    _parity_fminus_substep!(Y, method, field, τ / 2)
    return Y
end
