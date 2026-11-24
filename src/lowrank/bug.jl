"""
    BUG()

The fixed-rank basis-update & Galerkin integrator of Ceruti & Lubich 2022
(doi 10.1007/s10543-021-00873-0), a reference method of dynamical low rank. One step of length
`τ` from `X₀ S₀ V₀ᵀ`, with ``F`` the right-hand side of [`_grid_1d1v_rhs!`](@ref) and the
weights ``h_x``, ``h_v`` of the grid:

 1. the K-step ``\\dot{K} = h_v F(K V_0^\\mathsf{T}) V_0`` from `X₀ S₀` and the L-step
    ``\\dot{L} = h_x F(X_0 L^\\mathsf{T})^\\mathsf{T} X_0`` from `V₀ S₀ᵀ`, both over `τ`;
 2. `X₁` and `V₁`, the weighted orthonormal bases of `K(τ)` and `L(τ)`;
 3. the S-step ``\\dot{S} = h_x h_v X_1^\\mathsf{T} F(X_1 S V_1^\\mathsf{T}) V_1`` over `τ` from
    `M S₀ Nᵀ`, with `M = h_x X₁ᵀX₀` and `N = h_v V₁ᵀV₀`.

The rank stays that of the input. [`lowrank_step!`](@ref) takes the step.
"""
struct BUG end

"""
    AugmentedBUG()
    AugmentedBUG(tol; rmin = 1, rmax)

The augmented basis-update & Galerkin integrator of Ceruti, Kusch & Lubich 2022
(doi 10.1007/s10543-021-00907-7), a reference method of dynamical low rank. One step of length
`τ` from `X₀ S₀ V₀ᵀ` takes the K- and L-steps of [`BUG`](@ref), then

 1. `X̂` and `V̂`, the weighted orthonormal bases of `[K(τ) X₀]` and `[L(τ) V₀]`, of
    `min(2r, nx)` and `min(2r, nv)` columns;
 2. the S-step on `X̂` and `V̂` over `τ` from `(h_x X̂ᵀX₀) S₀ (h_v V₀ᵀV̂)`;
 3. the SVD `Ŝ(τ) = P Σ Qᵀ` and the truncation `X₁ = X̂ P₁`, `S₁ = Diagonal(σ₁, …, σ_{r₁})`,
    `V₁ = V̂ Q₁` to the first `r₁` singular vectors, so `S₁` is diagonal and non-increasing.

`AugmentedBUG()` truncates to the rank `r` of the step's input. `AugmentedBUG(tol; rmin, rmax)`
truncates to `r₁ = clamp(r*, rmin, min(rmax, k))`, with `k` the number of singular values of
`Ŝ(τ)` and `r*` the least `r` with ``(\\sum_{j > r} \\sigma_j^2)^{1/2} \\le`` `tol`: `tol` is an
absolute tolerance in the weighted Frobenius norm, which is the Frobenius norm of `S` for
orthonormal bases. A step of this method throws an `ArgumentError` before any write when the
rank of its input lies outside `rmin:rmax`.

Throws an `ArgumentError` unless `tol` is a finite real `≥ 0` and `rmin` and `rmax` are integers
with `1 ≤ rmin ≤ rmax`.
[`lowrank_step!`](@ref) takes the step.
"""
struct AugmentedBUG{TT <: Union{Nothing, Real}}
    tol::TT
    rmin::Int
    rmax::Int

    AugmentedBUG() = new{Nothing}(nothing, 1, typemax(Int))

    function AugmentedBUG(tol; rmin = 1, rmax)
        tol isa Real && isfinite(tol) && tol ≥ 0 ||
            throw(ArgumentError("the tolerance of AugmentedBUG is a finite real ≥ 0; got $tol"))
        rmin isa Integer && rmax isa Integer && 1 ≤ rmin ≤ rmax || throw(ArgumentError(
            "AugmentedBUG takes integers 1 ≤ rmin ≤ rmax; got rmin = $rmin and rmax = $rmax"))
        new{typeof(tol)}(tol, rmin, rmax)
    end
end

"""
    MidpointBUG()

The midpoint basis-update & Galerkin integrator, "Midpoint BUG (4r)" of Ceruti, Einkemmer,
Kusch & Lubich 2024 (doi 10.1007/s10543-024-01032-x; arXiv 2402.08607v1, §3), a reference
method of dynamical low rank, truncated to the rank `r` of its input. One step of length `τ`
from `Y₀ = X₀ S₀ V₀ᵀ`:

 1. the step of [`AugmentedBUG`](@ref) over `τ/2`, without its truncation, gives `Û`, `Ŝ½`,
    `V̂` and `Ŷ½ = Û Ŝ½ V̂ᵀ`;
 2. `Ū` and `V̄`, the weighted orthonormal bases of `[Û, τ h_v F(Ŷ½) V̂]` and
    `[V̂, τ h_x F(Ŷ½)ᵀ Û]`, of up to `4r` columns, and the Galerkin step
    ``\\dot{S} = h_x h_v \\bar{U}^\\mathsf{T} F(\\bar{U} S \\bar{V}^\\mathsf{T}) \\bar{V}`` over `τ`
    from `(h_x ŪᵀX₀) S₀ (h_v V₀ᵀV̄)`;
 3. the SVD of `S̄(τ)` and the truncation to `r`, as in [`AugmentedBUG`](@ref).

[`lowrank_step!`](@ref) takes the step.
"""
struct MidpointBUG end

const BUGMethod = Union{BUG, AugmentedBUG, MidpointBUG}

_weights(space::GridFactorSpace) = (step(space.x), step(space.v))

# Q / √h of the Householder QR of √h A: h QᵀQ = I with min(size(A)...) columns. A matrix of more
# columns than rows gets the full Q of its grid direction, which `orthonormalize_x` refuses.
_basis(A, h) = Matrix(qr(sqrt(h) .* A).Q) ./ sqrt(h)

_pad(A, n) = hcat(A, zeros(eltype(A), size(A, 1), n - size(A, 2)))

# The S-step over [t, t + τ] from the p × q core S₀ on X of p and V of q columns. The S-field of
# the workspace takes a square core, so the narrower factor and S₀ get zero columns and rows up
# to n = max(p, q): for finite entries they add exact zeros to the field and stay zero, and the
# p × q block is returned. The fields are autonomous, so the substeps do not read the time.
function _sstep(X, S₀, V, w, t, τ)
    p, q = size(S₀)
    n = max(p, q)
    Xn, Vn = _pad(X, n), _pad(V, n)
    Sn = zeros(eltype(S₀), n, n)
    Sn[1:p, 1:q] .= S₀
    return _rk4((Ṡ, S, s) -> _sfield!(Ṡ, Xn, S, Vn, w), Sn, t, τ)[1:p, 1:q]
end

# the K-step from X S and the L-step from V Sᵀ, over [t, t + τ]
function _klsteps(X, S, V, w, t, τ)
    K = _rk4((K̇, K, s) -> _kfield!(K̇, K, V, w), X * S, t, τ)
    L = _rk4((L̇, L, s) -> _lfield!(L̇, X, L, w), V * S', t, τ)
    return K, L
end

# the augmented BUG step, untruncated
function _augmented(X, S, V, w, hx, hv, t, τ)
    K, L = _klsteps(X, S, V, w, t, τ)
    X̂, V̂ = _basis(hcat(K, X), hx), _basis(hcat(L, V), hv)
    Ŝ₀ = (hx .* (X̂' * X)) * S * (hv .* (V' * V̂))
    return X̂, _sstep(X̂, Ŝ₀, V̂, w, t, τ), V̂
end

# the least r with ‖σ[r+1:end]‖ ≤ tol; r = k gives the empty tail, so one exists
function _least_rank(σ, tol)
    first(r for r in 0:length(σ) if norm(view(σ, (r + 1):length(σ))) ≤ tol)
end

_truncation_rank(::Union{AugmentedBUG{Nothing}, MidpointBUG}, σ, r) = min(r, length(σ))
function _truncation_rank(m::AugmentedBUG{<:Real}, σ, r)
    clamp(_least_rank(σ, m.tol), m.rmin, min(m.rmax, length(σ)))
end

# The SVD of Ŝ and the truncation to the rank of `method`. `svd` throws for a non-finite matrix,
# so a non-finite Ŝ skips the SVD: the step installs the first r columns of the bases and an
# r × r core of NaN, a state of the input's rank whose status `all(isfinite, S)` is false.
function _truncate(method, X̂, Ŝ, V̂, r)
    all(isfinite, Ŝ) || return X̂[:, 1:r], fill(eltype(Ŝ)(NaN), r, r), V̂[:, 1:r]
    F = svd(Ŝ)
    r₁ = _truncation_rank(method, F.S, r)
    return X̂ * F.U[:, 1:r₁], Matrix(Diagonal(F.S[1:r₁])), V̂ * F.V[:, 1:r₁]
end

function _bug_step(::BUG, X, S, V, w, hx, hv, t, τ)
    K, L = _klsteps(X, S, V, w, t, τ)
    X₁, V₁ = _basis(K, hx), _basis(L, hv)
    S₀ = (hx .* (X₁' * X)) * S * (hv .* (V' * V₁))
    return X₁, _sstep(X₁, S₀, V₁, w, t, τ), V₁
end

function _bug_step(method::AugmentedBUG, X, S, V, w, hx, hv, t, τ)
    return _truncate(method, _augmented(X, S, V, w, hx, hv, t, τ)..., size(S, 1))
end

function _bug_step(method::MidpointBUG, X, S, V, w, hx, hv, t, τ)
    Û, Ŝ½, V̂ = _augmented(X, S, V, w, hx, hv, t, τ / 2)
    K½, L½ = Û * Ŝ½, V̂ * Ŝ½'
    FV = _kfield!(similar(K½), K½, V̂, w)       # h_v F(Ŷ½) V̂
    FU = _lfield!(similar(L½), Û, L½, w)       # h_x F(Ŷ½)ᵀ Û
    Ū, V̄ = _basis(hcat(Û, τ .* FV), hx), _basis(hcat(V̂, τ .* FU), hv)
    S̄₀ = (hx .* (Ū' * X)) * S * (hv .* (V' * V̄))
    return _truncate(method, Ū, _sstep(Ū, S̄₀, V̄, w, t, τ), V̄, size(S, 1))
end

function _check_rank_bounds(m::AugmentedBUG{<:Real}, r)
    m.rmin ≤ r ≤ m.rmax || throw(ArgumentError(
        "the rank $r of the state is not in rmin:rmax = $(m.rmin):$(m.rmax) of the method"))
    return nothing
end

"""
    lowrank_step!(Y, method::Union{BUG, AugmentedBUG, MidpointBUG}, field, t, τ)

The step of [`BUG`](@ref), [`AugmentedBUG`](@ref) or [`MidpointBUG`](@ref) on the
[`LowRankDistribution`](@ref) `Y` on a [`GridFactorSpace`](@ref). `field` is the workspace
[`VlasovMethods.LowRankFields1D1V`](@ref) of the same element type `T`, whose K-, L- and S-fields
evaluate the full-order field on the factors, without the `nx × nv` array. The fields are
autonomous, so `t` is not read. Each K-, L-, S- and Galerkin substep is one step of the classical
RK4 tableau over its interval. The contract of every step, the conversion of `t` and `τ` and the
refusal of a non-finite `τ`, is that of the generic [`lowrank_step!`](@ref).

`τ = 0` keeps `X S Vᵀ` to round-off, except for `AugmentedBUG(tol; rmin, rmax)`, whose truncation
still applies and changes `X S Vᵀ` by up to `tol` in the weighted Frobenius norm. The step counts
every column of `V` in the rank, the fixed modes of `Y` included, and sets `nfixed = 0`, because
no reference method keeps the fixed modes in the span of `V`.

An S-step or a Galerkin step that leaves a non-finite value does not throw. `BUG()` keeps its
factors of rank `r`, non-finite as they come; the augmented and the midpoint BUG skip the
truncation and install the first `r` columns of the augmented bases and an `r × r` core of `NaN`,
with `r` the rank of the input. The status that a caller tests is `all(isfinite, Y.S)`.
"""
function lowrank_step!(Y::LowRankDistribution{T, <:GridFactorSpace}, method::BUGMethod,
        field::LowRankFields1D1V{T}, t, τ::Real) where {T}
    h = _finite_step(τ, T)
    _check_rank_bounds(method, size(Y.S, 1))
    hx, hv = _weights(Y.space)
    Y.X, Y.S, Y.V = _bug_step(method, Y.X, Y.S, Y.V, field, hx, hv, T(t), h)
    Y.nfixed = 0
    return Y
end

# the columns of the stored factors of a run: the largest rank that the method can reach
_stored_rank(m::AugmentedBUG{<:Real}, Y) = min(m.rmax, size(Y)...)
