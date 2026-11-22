# The interface through which a step reads its field: the K-, S- and L-right-hand sides of a
# factored state at the time t, and for the Strang composition the potential of a factored state
# and the L-right-hand side at a given, frozen potential. The S-right-hand side carries the plus
# sign of `_sfield!`, `Ṡ = hx hv Xᵀ F(X S Vᵀ, t) V`; the step integrates it backward. A field of
# another kind joins the step by methods of these functions. The fields of `LowRankFields1D1V`
# are autonomous, so their methods ignore t.

_lowrank_kfield!(K̇, K, V, w::LowRankFields1D1V, t) = _kfield!(K̇, K, V, w)
_lowrank_sfield!(Ṡ, X, S, V, w::LowRankFields1D1V, t) = _sfield!(Ṡ, X, S, V, w)
_lowrank_lfield!(L̇, X, L, w::LowRankFields1D1V, t) = _lfield!(L̇, X, L, w)
_lowrank_lfield!(L̇, X, L, w::LowRankFields1D1V, t, ϕ) = _lfield!(L̇, X, L, w, ϕ)

# the potential of f = X S Vᵀ, and of f = X Lᵀ
function _lowrank_potential(X, S, V, w::LowRankFields1D1V{T}, t) where {T}
    _potential!(Vector{T}(undef, w.arakawa.nx), X, S, V, w)
end
function _lowrank_potential(X, L, w::LowRankFields1D1V{T}, t) where {T}
    _potential!(Vector{T}(undef, w.arakawa.nx), X, L, w)
end

"""
    ProjectorSplitting(composition)

The projector-splitting integrator of dynamical low rank (Lubich and Oseledets 2014,
doi:10.1007/s10543-013-0454-0), a method for [`lowrank_step!`](@ref). `composition` is `:lie`, the
first-order K–S–L splitting, or `:strang`, the second-order splitting with the prediction of the
midpoint potential of Einkemmer and Lubich 2018 (doi:10.1137/18m116383x, Algorithm 1). The two
compositions are methods of two different types. Each substep is one step of the classical RK4
tableau over the substep's interval.

Throws an `ArgumentError` for any other `composition`, a `String` included.
"""
struct ProjectorSplitting{C}
    function ProjectorSplitting(composition)
        composition === :lie || composition === :strang ||
            throw(ArgumentError(
                "the composition of a projector splitting is :lie or :strang; got " *
                repr(composition)))
        new{composition}()
    end
end

# One step of the classical RK4 tableau of ẏ = rhs(y, s) over [t, t + τ], with
# `rhs!(ẏ, y, s)` writing ẏ and returning it. The substep solver of the projector splitting.
function _substep(::ProjectorSplitting, rhs!, y, t, τ)
    k1 = rhs!(similar(y), y, t)
    k2 = rhs!(similar(y), y .+ (τ / 2) .* k1, t + τ / 2)
    k3 = rhs!(similar(y), y .+ (τ / 2) .* k2, t + τ / 2)
    k4 = rhs!(similar(y), y .+ τ .* k3, t + τ)
    return y .+ (τ / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

# the backward S-substep: Ṡ = −hx hv Xᵀ F(X S Vᵀ, s) V
function _backward_sstep(method, X, S, V, field, t, τ)
    _substep(
        method, (Ṡ, S, s) -> (_lowrank_sfield!(Ṡ, X, S, V, field, s); Ṡ .= .-Ṡ), S, t, τ)
end

function _kstep(method, K, V, field, t, τ)
    _substep(method, (K̇, K, s) -> _lowrank_kfield!(K̇, K, V, field, s), K, t, τ)
end

function _lstep(method, X, L, field, t, τ)
    _substep(method, (L̇, L, s) -> _lowrank_lfield!(L̇, X, L, field, s), L, t, τ)
end

function _lstep(method, X, L, field, t, τ, ϕ)
    _substep(method, (L̇, L, s) -> _lowrank_lfield!(L̇, X, L, field, s, ϕ), L, t, τ)
end

# τ in the element type T of the state, refused where it is not finite in T
function _step_length(τ, ::Type{T}) where {T}
    τT = T(τ)
    isfinite(τT) ||
        throw(ArgumentError("the step τ of a projector splitting is finite; got $τT"))
    return τT
end

@doc raw"""
    lowrank_step!(Y::LowRankDistribution, method::ProjectorSplitting, field, t, τ)

Advance the low-rank state `Y` by one step of length `τ` from the time `t`, in place, and return
`Y`. `field` is the right-hand side, such as the [`LowRankFields1D1V`](@ref) of the state's grid.
With ``F`` the full-order right-hand side, the substeps are those of Lubich and Oseledets 2014:

```math
\dot{K} = h_v F(K V^\mathsf{T}) V , \qquad
\dot{S} = -h_x h_v X^\mathsf{T} F(X S V^\mathsf{T}) V , \qquad
\dot{L} = h_x F(X L^\mathsf{T})^\mathsf{T} X ,
```

each one RK4 step. After a K-substep, `X, S = orthonormalize_x(Y.space, K)`; after an L-substep,
`V, R = orthonormalize_v(Y.space, L)` and `S = Rᵀ`, so `X` and `V` are orthonormal in the
weighted products of the space, ``h_x X^\mathsf{T} X = I`` and ``h_v V^\mathsf{T} V = I``.

  - `:lie` runs the K-substep over ``[t, t + τ]`` from ``K = X S``, the backward S-substep over
    ``[t, t + τ]``, and the L-substep over ``[t, t + τ]`` from ``L = V S^\mathsf{T}``.
  - `:strang` runs Algorithm 1 of Einkemmer and Lubich 2018: the K- and the backward S-substep
    over ``[t, t + τ/2]``; the potential ``\phi^0`` of ``X S V^\mathsf{T}``; an L-substep over
    ``[t, t + τ/2]`` at the frozen ``\phi^0``, which predicts ``L^{1/2}``; the potential
    ``\phi^{1/2}`` of ``X (L^{1/2})^\mathsf{T}``; the L-substep over ``[t, t + τ]`` at the frozen
    ``\phi^{1/2}``; and the backward S- and the K-substep over ``[t + τ/2, t + τ]``. The K- and
    S-substeps take the self-consistent field.

The step needs no inverse of `S`, so a rank-deficient `S`, the zero state included, is taken. It
forms no ``n_x \times n_v`` array, writes only `Y`, assigns new factors to it, and sets
`Y.nfixed = 0`: the splitting does not keep fixed velocity modes in the leading columns of `V`.

`τ` is a signed step, converted to the element type `T` of `Y` before its first use, as is `t`;
`τ = 0` returns the state unchanged up to the round-off of the orthonormalisations. Throws an
`ArgumentError` before any write when `τ` is not finite in `T`. The entries of `Y` are not
checked: a non-finite entry propagates.

A field of another kind joins the step through the methods `_lowrank_kfield!(K̇, K, V, field, t)`,
`_lowrank_sfield!(Ṡ, X, S, V, field, t)` (with the plus sign, ``h_x h_v X^\mathsf{T} F V``) and
`_lowrank_lfield!(L̇, X, L, field, t)`, and for `:strang` also
`_lowrank_potential(X, S, V, field, t)`, `_lowrank_potential(X, L, field, t)` and
`_lowrank_lfield!(L̇, X, L, field, t, ϕ)`.
"""
function lowrank_step!(
        Y::LowRankDistribution, method::ProjectorSplitting{:lie}, field, t, τ)
    T = eltype(Y)
    τ = _step_length(τ, T)
    t = T(t)
    K = _kstep(method, Y.X * Y.S, Y.V, field, t, τ)
    X, S₁ = orthonormalize_x(Y.space, K)
    S₂ = _backward_sstep(method, X, S₁, Y.V, field, t, τ)
    L = _lstep(method, X, Y.V * S₂', field, t, τ)
    V, R = orthonormalize_v(Y.space, L)
    Y.X, Y.S, Y.V, Y.nfixed = X, permutedims(R), V, 0
    return Y
end

function lowrank_step!(
        Y::LowRankDistribution, method::ProjectorSplitting{:strang}, field, t, τ)
    T = eltype(Y)
    τ = _step_length(τ, T)
    t = T(t)
    τ½ = τ / 2
    t½ = t + τ½
    V₀ = Y.V
    K₁ = _kstep(method, Y.X * Y.S, V₀, field, t, τ½)
    X₁, S₁ = orthonormalize_x(Y.space, K₁)
    S₂ = _backward_sstep(method, X₁, S₁, V₀, field, t, τ½)
    L₀ = V₀ * S₂'
    ϕ₀ = _lowrank_potential(X₁, S₂, V₀, field, t)
    L½ = _lstep(method, X₁, L₀, field, t, τ½, ϕ₀)
    ϕ½ = _lowrank_potential(X₁, L½, field, t½)
    L₁ = _lstep(method, X₁, L₀, field, t, τ, ϕ½)
    V₁, R = orthonormalize_v(Y.space, L₁)
    S₃ = _backward_sstep(method, X₁, permutedims(R), V₁, field, t½, τ½)
    K₂ = _kstep(method, X₁ * S₃, V₁, field, t½, τ½)
    X₂, S₄ = orthonormalize_x(Y.space, K₂)
    Y.X, Y.S, Y.V, Y.nfixed = X₂, S₄, V₁, 0
    return Y
end
