@doc raw"""
    SplineLowRankFields1D1V(space::SplineFactorSpace)
    SplineLowRankFields1D1V(d::LowRankDistribution)

The workspace of the K-, L- and S-fields of 1D1V Vlasov–Poisson on the
[`SplineFactorSpace`](@ref) `space`, in its element type `T`: the space, the
`PoissonSolverSpline` of its x basis, the values and the first derivatives of the two bases at
the nodes of their quadratures, the weights of the quadratures, the kinetic energy ``v^2/2`` and
its derivative ``v`` at the velocity nodes, and the integrals ``\iota_x = (\int \varphi_k dx)_k``
and ``\iota_v = (\int \psi_j dv)_j`` of the basis functions. The second form is the first on the
space of `d`.

The fields are the projections of the field of the spline Galerkin model
[`VlasovPoisson`](@ref) on a [`TensorSplineDistribution`](@ref),
``\dot F = M_x^{-1} B(F) M_v^{-1}`` with ``B(F)_{ij} = \bar T(\varphi_i \psi_j, h, f)``, the
cyclic average [`_averaged_form`](@ref) and ``h = \phi \otimes 1 + 1 \otimes v^2/2``. The
methods of [`_kfield!`](@ref), [`_lfield!`](@ref), [`_sfield!`](@ref) and
[`_potential!`](@ref) take the workspace, and the projector splitting steps a state of the space
with it, through `lowrank_step!(Y, ProjectorSplitting(c), w, t, τ)`.

The constructor lets the `ArgumentError` of `PoissonSolverSpline` propagate, where PoissonSolvers
refuses the basis: a fine periodic grid in `Float32`, such as 768 cells of ``[0, 4\pi)``.
"""
struct SplineLowRankFields1D1V{T, ST <: SplineFactorSpace{T}, PT <: PoissonSolverSpline{T},
    BT <: AbstractMatrix{T}}
    space::ST
    solver::PT
    Φ0::BT
    Φ1::BT
    Ψ0::BT
    Ψ1::BT
    wx::Vector{T}
    wv::Vector{T}
    k0::Vector{T}
    k1::Vector{T}
    ιx::Vector{T}
    ιv::Vector{T}

    function SplineLowRankFields1D1V(space::SplineFactorSpace{T}) where {T}
        solver = PoissonSolverSpline(space.xbasis)
        xq, vq = space.xquadrature, space.vquadrature
        v = Vector{T}(quadrature_nodes(vq))
        tables = map(((q, d),) -> permutedims(basis_values(q, d)),
            ((xq, 0), (xq, 1), (vq, 0), (vq, 1)))
        new{T, typeof(space), typeof(solver), eltype(tables)}(space, solver, tables...,
            Vector{T}(quadrature_weights(xq)), Vector{T}(quadrature_weights(vq)),
            v .^ 2 ./ 2, v, Vector{T}(basis_integrals(xq)), Vector{T}(basis_integrals(vq)))
    end
end

function SplineLowRankFields1D1V(d::LowRankDistribution{T, <:SplineFactorSpace}) where {T}
    SplineLowRankFields1D1V(d.space)
end

# The factored form of the cyclic average. For a test family a = a_x ⊗ a_v, a factor
# f = Σ_α f_x,α ⊗ f_v,α and h = Σ_t h_x,t ⊗ h_v,t,
#   T(a₁⊗a₂, b₁⊗b₂, c₁⊗c₂) = (∫a₁ b₁' c₁)(∫a₂ b₂ c₂') − (∫a₁ b₁ c₁')(∫a₂ b₂' c₂),
# so each cyclic term is a sum of products of one-dimensional integrals of three factors. With
# I(d_a, d_h, d_f) the integral of the derivatives of orders d_a, d_h, d_f, and P = I(0, 1, 0),
# Q = I(0, 0, 1), R = I(1, 0, 0) in x and p, q, r in v, the average of one term of h is
#   3 B = P (q − r)ᵀ + Q (r − p)ᵀ + R (p − q)ᵀ.
# The term ϕ ⊗ 1 of h has p = 0, the term 1 ⊗ v²/2 has P = 0, which leaves five products.

# the x integrals P, Q, R of the term ϕ ⊗ 1 and Q, R of the term 1 ⊗ v²/2, for the test family
# of values a0 and derivatives a1 at the x nodes and the x factors of values f0 and derivatives
# f1, each a matrix of one column per test function and per factor
function _x_integrals(w, a0, a1, f0, f1, ϕ)
    wx = w.wx
    wϕ0, wϕ1 = wx .* (w.Φ0 * ϕ), wx .* (w.Φ1 * ϕ)
    P1 = a0' * (wϕ1 .* f0)
    Q1 = a0' * (wϕ0 .* f1)
    R1 = a1' * (wϕ0 .* f0)
    Q2 = a0' * (wx .* f1)
    R2 = a1' * (wx .* f0)
    return (P1, Q1, R1, Q2, R2)
end

# the v integrals q, r of the term ϕ ⊗ 1 and p, q, r of the term 1 ⊗ v²/2
function _v_integrals(w, a0, a1, f0, f1)
    wv = w.wv
    wk0, wk1 = wv .* w.k0, wv .* w.k1
    q1 = a0' * (wv .* f1)
    r1 = a1' * (wv .* f0)
    p2 = a0' * (wk1 .* f0)
    q2 = a0' * (wk0 .* f1)
    r2 = a1' * (wk0 .* f0)
    return (q1, r1, p2, q2, r2)
end

# the five products a bᵀ of 3 B, from the x and the v integrals
function _averaged_terms((P1, Q1, R1, Q2, R2), (q1, r1, p2, q2, r2))
    return ((P1, q1 .- r1), (Q1, r1), (R1, .-q1), (Q2, r2 .- p2), (R2, p2 .- q2))
end

# B = Σ a bᵀ / 3 over the products
function _sum_products!(B, products)
    fill!(B, zero(eltype(B)))
    for (a, b) in products
        mul!(B, a, transpose(b), true, true)
    end
    B ./= 3
    return B
end

# B ← M⁻¹ B for M = UᵀU
_mass_solve!(U, B) = ldiv!(UpperTriangular(U), ldiv!(transpose(UpperTriangular(U)), B))

function _check_spline_kfield(K̇, K, V, w, ϕ = nothing)
    nx, nv = size(w.space)
    r = size(K, 2)
    size(K̇) == size(K) == (nx, r) && size(V) == (nv, r) &&
    (ϕ === nothing || length(ϕ) == nx) ||
        throw(DimensionMismatch(
            "_kfield! on a spline space of $nx × $nv functions takes K̇ and K of $nx rows and " *
            "V of $nv rows, with $r columns each$(_ϕ_wanted(ϕ, nx)); got K̇ of size " *
            "$(size(K̇)), K of size $(size(K)) and V of size $(size(V))$(_ϕ_got(ϕ))"))
    return nothing
end

function _check_spline_lfield(L̇, X, L, w, ϕ = nothing)
    nx, nv = size(w.space)
    r = size(L, 2)
    size(L̇) == size(L) == (nv, r) && size(X) == (nx, r) &&
    (ϕ === nothing || length(ϕ) == nx) ||
        throw(DimensionMismatch(
            "_lfield! on a spline space of $nx × $nv functions takes L̇ and L of $nv rows and " *
            "X of $nx rows, with $r columns each$(_ϕ_wanted(ϕ, nx)); got L̇ of size " *
            "$(size(L̇)), X of size $(size(X)) and L of size $(size(L))$(_ϕ_got(ϕ))"))
    return nothing
end

function _check_spline_sfield(Ṡ, X, S, V, w, ϕ = nothing)
    nx, nv = size(w.space)
    r = size(X, 2)
    size(X) == (nx, r) && size(V) == (nv, r) && size(Ṡ) == size(S) == (r, r) &&
    (ϕ === nothing || length(ϕ) == nx) ||
        throw(DimensionMismatch(
            "_sfield! on a spline space of $nx × $nv functions takes X of $nx rows, V of $nv " *
            "rows, and Ṡ and S of size r × r, with r the columns of X and V" *
            "$(_ϕ_wanted(ϕ, nx)); got Ṡ of size $(size(Ṡ)), X of size $(size(X)), S of size " *
            "$(size(S)) and V of size $(size(V))$(_ϕ_got(ϕ))"))
    return nothing
end

@doc raw"""
    _kfield!(K̇, K, V, w::SplineLowRankFields1D1V)
    _kfield!(K̇, K, V, w::SplineLowRankFields1D1V, ϕ)

Write the K-field of the spline Galerkin model into `K̇` and return it:

```math
\dot{K} = \dot F(K V^\mathsf{T}) \, M_v V
= M_x^{-1} \big[ \bar T(\varphi_i \otimes v_\gamma, h, K V^\mathsf{T}) \big]_{i \gamma} ,
\qquad v_\gamma = \sum_j V_{j\gamma} \psi_j ,
```

with ``\dot F`` the field of [`SplineLowRankFields1D1V`](@ref), `K` and `K̇` of size `nx × r` and
`V` of size `nv × r`. Without `ϕ` the potential is that of ``f = K V^\mathsf{T}``, from
[`_potential!`](@ref); with `ϕ`, the coefficient `Vector` of length `nx` of a potential in the x
basis, ``h = \phi \otimes 1 + 1 \otimes v^2/2`` takes the given `ϕ`, frozen. The factors are
`Matrix`es of the element type `T` of `w`; `V` need not be orthonormal.

``h`` has rank 2 and ``f`` rank `r`, so the form is a sum of products of one-dimensional
quadratures of the factors, at a cost of ``O((n_x + n_v) r^2)`` beyond the Poisson solve, without
the ``n_x \times n_v`` array.

Throws a `DimensionMismatch` before any write when the sizes do not fit the space of `w` and the
rank `r = size(K, 2)`, and an `ArgumentError` when `K̇` shares memory with `K`, `V` or `ϕ`. At
`r = 0` it returns the empty `K̇`. A non-finite entry propagates to the result, unchecked.
"""
function _kfield!(K̇::Matrix{T}, K::Matrix{T}, V::Matrix{T},
        w::SplineLowRankFields1D1V{T}) where {T}
    _check_spline_kfield(K̇, K, V, w)
    ϕ = _potential!(Vector{T}(undef, size(w.space)[1]), K, V, w)
    return _kfield!(K̇, K, V, w, ϕ)
end

function _kfield!(K̇::Matrix{T}, K::Matrix{T}, V::Matrix{T}, w::SplineLowRankFields1D1V{T},
        ϕ::Vector{T}) where {T}
    _check_spline_kfield(K̇, K, V, w, ϕ)
    _check_alias("_kfield!", "K̇", K̇, (K = K, V = V, ϕ = ϕ))
    V0, V1 = w.Ψ0 * V, w.Ψ1 * V
    x = _x_integrals(w, w.Φ0, w.Φ1, w.Φ0 * K, w.Φ1 * K, ϕ)
    v = _v_integrals(w, V0, V1, V0, V1)
    _sum_products!(K̇, _averaged_terms(x, v))
    return _mass_solve!(w.space.Ux, K̇)
end

@doc raw"""
    _lfield!(L̇, X, L, w::SplineLowRankFields1D1V)
    _lfield!(L̇, X, L, w::SplineLowRankFields1D1V, ϕ)

Write the L-field of the spline Galerkin model into `L̇` and return it:

```math
\dot{L} = \dot F(X L^\mathsf{T})^\mathsf{T} M_x X
= M_v^{-1} \big[ \bar T(x_\gamma \otimes \psi_j, h, X L^\mathsf{T}) \big]_{j \gamma} ,
\qquad x_\gamma = \sum_i X_{i\gamma} \varphi_i ,
```

with `X` of size `nx × r`, and `L` and `L̇` of size `nv × r`. Without `ϕ` the potential is that
of ``f = X L^\mathsf{T}``, from [`_potential!`](@ref)`(ϕ, X, L, w)`; with `ϕ`, the potential is
the given `ϕ`, frozen. `X` need not be orthonormal. The evaluation, the edges and the errors are
those of the method of [`_kfield!`](@ref) for the workspace, with `L̇` in place of `K̇` and the
rank `r = size(L, 2)`.
"""
function _lfield!(L̇::Matrix{T}, X::Matrix{T}, L::Matrix{T},
        w::SplineLowRankFields1D1V{T}) where {T}
    _check_spline_lfield(L̇, X, L, w)
    ϕ = _potential!(Vector{T}(undef, size(w.space)[1]), X, L, w)
    return _lfield!(L̇, X, L, w, ϕ)
end

function _lfield!(L̇::Matrix{T}, X::Matrix{T}, L::Matrix{T}, w::SplineLowRankFields1D1V{T},
        ϕ::Vector{T}) where {T}
    _check_spline_lfield(L̇, X, L, w, ϕ)
    _check_alias("_lfield!", "L̇", L̇, (X = X, L = L, ϕ = ϕ))
    X0, X1 = w.Φ0 * X, w.Φ1 * X
    x = _x_integrals(w, X0, X1, X0, X1, ϕ)
    v = _v_integrals(w, w.Ψ0, w.Ψ1, w.Ψ0 * L, w.Ψ1 * L)
    _sum_products!(L̇, map(reverse, _averaged_terms(x, v)))
    return _mass_solve!(w.space.Uv, L̇)
end

@doc raw"""
    _sfield!(Ṡ, X, S, V, w::SplineLowRankFields1D1V)
    _sfield!(Ṡ, X, S, V, w::SplineLowRankFields1D1V, ϕ)

Write the S-field of the spline Galerkin model into `Ṡ` and return it:

```math
\dot{S} = X^\mathsf{T} M_x \dot F(X S V^\mathsf{T}) \, M_v V
= \big[ \bar T(x_\alpha \otimes v_\beta, h, X S V^\mathsf{T}) \big]_{\alpha \beta} ,
```

with `X` of size `nx × r`, `V` of size `nv × r`, and `S` and `Ṡ` of size `r × r`. Without `ϕ`
the potential is that of ``f = X S V^\mathsf{T}``, from
[`_potential!`](@ref)`(ϕ, X, S, V, w)`; with `ϕ`, the potential is the given `ϕ`, frozen. `X`
and `V` need not be orthonormal. At `r = 1` the field vanishes, up to round-off, by the total
antisymmetry of ``\bar T``. The evaluation, the edges and the errors are those of the method of
[`_kfield!`](@ref) for the workspace, with `Ṡ` in place of `K̇` and the rank `r = size(X, 2)`.
"""
function _sfield!(Ṡ::Matrix{T}, X::Matrix{T}, S::Matrix{T}, V::Matrix{T},
        w::SplineLowRankFields1D1V{T}) where {T}
    _check_spline_sfield(Ṡ, X, S, V, w)
    ϕ = _potential!(Vector{T}(undef, size(w.space)[1]), X, S, V, w)
    return _sfield!(Ṡ, X, S, V, w, ϕ)
end

function _sfield!(Ṡ::Matrix{T}, X::Matrix{T}, S::Matrix{T}, V::Matrix{T},
        w::SplineLowRankFields1D1V{T}, ϕ::Vector{T}) where {T}
    _check_spline_sfield(Ṡ, X, S, V, w, ϕ)
    _check_alias("_sfield!", "Ṡ", Ṡ, (X = X, S = S, V = V, ϕ = ϕ))
    X0, X1 = w.Φ0 * X, w.Φ1 * X
    V0, V1 = w.Ψ0 * V, w.Ψ1 * V
    x = _x_integrals(w, X0, X1, X0 * S, X1 * S, ϕ)
    v = _v_integrals(w, V0, V1, V0, V1)
    return _sum_products!(Ṡ, _averaged_terms(x, v))
end

@doc raw"""
    _potential!(ϕ, K, V, w::SplineLowRankFields1D1V)
    _potential!(ϕ, X, S, V, w::SplineLowRankFields1D1V)

Write the coefficients of the potential of factored data into the `Vector` `ϕ` of length `nx`
and return it, in the convention of [`VlasovPoisson`](@ref) on a
[`TensorSplineDistribution`](@ref): the `PoissonSolverSpline` solve of the load vector
``M_x \rho - \bar\rho \, \iota_x``, shifted to ``\iota_x^\mathsf{T} \phi = 0``, with
``\bar\rho = 1^\mathsf{T} M_x \rho / L`` for the length ``L = 1^\mathsf{T} \iota_x`` of the
x domain and the density coefficients

```math
\rho = K (V^\mathsf{T} \iota_v) \quad \text{of } f = K V^\mathsf{T} , \qquad
\rho = X S (V^\mathsf{T} \iota_v) \quad \text{of } f = X S V^\mathsf{T} ,
```

``\iota_v = (\int \psi_j dv)_j``, without the ``n_x \times n_v`` array. The potential of
``f = X L^\mathsf{T}`` is `_potential!(ϕ, X, L, w)`. The factors are `Matrix`es of the element
type `T` of `w`.

Throws a `DimensionMismatch` before any write when the sizes do not fit the space of `w` and the
rank, and an `ArgumentError` when `ϕ` shares memory with a factor. At the rank 0 it writes zeros.
A non-finite entry propagates to the result, unchecked.
"""
function _potential!(ϕ::Vector{T}, K::Matrix{T}, V::Matrix{T},
        w::SplineLowRankFields1D1V{T}) where {T}
    nx, nv = size(w.space)
    length(ϕ) == size(K, 1) == nx && size(V, 1) == nv && size(K, 2) == size(V, 2) ||
        throw(DimensionMismatch(
            "_potential! on a spline space of $nx × $nv functions takes ϕ of length $nx, K " *
            "of $nx rows and V of $nv rows, with as many columns each; got ϕ of length " *
            "$(length(ϕ)), K of size $(size(K)) and V of size $(size(V))"))
    _check_alias("_potential!", "ϕ", ϕ, (K = K, V = V))
    return _spline_density_potential!(ϕ, K * (V' * w.ιv), w)
end

function _potential!(ϕ::Vector{T}, X::Matrix{T}, S::Matrix{T}, V::Matrix{T},
        w::SplineLowRankFields1D1V{T}) where {T}
    nx, nv = size(w.space)
    r = size(X, 2)
    length(ϕ) == nx && size(X) == (nx, r) && size(V) == (nv, r) && size(S) == (r, r) ||
        throw(DimensionMismatch(
            "_potential! on a spline space of $nx × $nv functions takes ϕ of length $nx, X of " *
            "$nx rows, V of $nv rows and S of size r × r, with r the columns of X and V; got " *
            "ϕ of length $(length(ϕ)), X of size $(size(X)), S of size $(size(S)) and V of " *
            "size $(size(V))"))
    _check_alias("_potential!", "ϕ", ϕ, (X = X, S = S, V = V))
    return _spline_density_potential!(ϕ, X * (S * (V' * w.ιv)), w)
end

# ϕ of the density coefficients ρ: the mean-free solve of the load Mx ρ − ρ̄ ι_x, then the shift
# to ι_xᵀ ϕ = 0, as `_load_potential` of the model does
function _spline_density_potential!(ϕ, ρ, w)
    ι = w.ιx
    L = sum(ι)
    b = w.space.Mx * ρ
    b .-= (sum(b) / L) .* ι
    PoissonSolvers.solve!(ϕ, w.solver, b)
    ϕ .-= dot(ι, ϕ) / L
    return ϕ
end

# The projector splitting reads the fields of the workspace through these methods. The fields
# are autonomous, so the methods ignore t.
_lowrank_kfield!(K̇, K, V, w::SplineLowRankFields1D1V, t) = _kfield!(K̇, K, V, w)
_lowrank_sfield!(Ṡ, X, S, V, w::SplineLowRankFields1D1V, t) = _sfield!(Ṡ, X, S, V, w)
_lowrank_lfield!(L̇, X, L, w::SplineLowRankFields1D1V, t) = _lfield!(L̇, X, L, w)
_lowrank_lfield!(L̇, X, L, w::SplineLowRankFields1D1V, t, ϕ) = _lfield!(L̇, X, L, w, ϕ)

# the potential of f = X S Vᵀ, and of f = X Lᵀ
function _lowrank_potential(X, S, V, w::SplineLowRankFields1D1V{T}, t) where {T}
    _potential!(Vector{T}(undef, size(w.space)[1]), X, S, V, w)
end
function _lowrank_potential(X, L, w::SplineLowRankFields1D1V{T}, t) where {T}
    _potential!(Vector{T}(undef, size(w.space)[1]), X, L, w)
end
