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

# one step of the classical RK4 tableau, b = (1, 2, 2, 1)/6, over τ for ẏ = g!(ẏ, y)
function _rk4(g!, y₀::Matrix{T}, τ::T) where {T}
    k₁ = g!(similar(y₀), y₀)
    k₂ = g!(similar(y₀), y₀ .+ (τ / 2) .* k₁)
    k₃ = g!(similar(y₀), y₀ .+ (τ / 2) .* k₂)
    k₄ = g!(similar(y₀), y₀ .+ τ .* k₃)
    return y₀ .+ (τ / 6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
end

_weights(space::GridFactorSpace) = (step(space.x), step(space.v))

# Q / √h of the Householder QR of √h A: h QᵀQ = I with min(size(A)...) columns. A matrix of more
# columns than rows gets the full Q of its grid direction, which `orthonormalize_x` refuses.
_basis(A, h) = Matrix(qr(sqrt(h) .* A).Q) ./ sqrt(h)

_pad(A, n) = hcat(A, zeros(eltype(A), size(A, 1), n - size(A, 2)))

# The S-step over τ from the p × q core S₀ on X of p and V of q columns. The S-field of the
# workspace takes a square core, so the narrower factor and S₀ get zero columns and rows up to
# n = max(p, q): for finite entries they add exact zeros to the field and stay zero, and the
# p × q block is returned.
function _sstep(X, S₀, V, w, τ)
    p, q = size(S₀)
    n = max(p, q)
    Xn, Vn = _pad(X, n), _pad(V, n)
    Sn = zeros(eltype(S₀), n, n)
    Sn[1:p, 1:q] .= S₀
    return _rk4((Ṡ, S) -> _sfield!(Ṡ, Xn, S, Vn, w), Sn, τ)[1:p, 1:q]
end

# the K-step from X S and the L-step from V Sᵀ, over τ
function _klsteps(X, S, V, w, τ)
    K = _rk4((K̇, K) -> _kfield!(K̇, K, V, w), X * S, τ)
    L = _rk4((L̇, L) -> _lfield!(L̇, X, L, w), V * S', τ)
    return K, L
end

# the augmented BUG step, untruncated
function _augmented(X, S, V, w, hx, hv, τ)
    K, L = _klsteps(X, S, V, w, τ)
    X̂, V̂ = _basis(hcat(K, X), hx), _basis(hcat(L, V), hv)
    Ŝ₀ = (hx .* (X̂' * X)) * S * (hv .* (V' * V̂))
    return X̂, _sstep(X̂, Ŝ₀, V̂, w, τ), V̂
end

# the least r with ‖σ[r+1:end]‖ ≤ tol; r = k gives the empty tail, so one exists
function _least_rank(σ, tol)
    first(r for r in 0:length(σ) if norm(view(σ, (r + 1):length(σ))) ≤ tol)
end

_truncation_rank(::Union{AugmentedBUG{Nothing}, MidpointBUG}, σ, r) = min(r, length(σ))
function _truncation_rank(m::AugmentedBUG{<:Real}, σ, r)
    clamp(_least_rank(σ, m.tol), m.rmin, min(m.rmax, length(σ)))
end

# The SVD of Ŝ and the truncation to the rank of `method`. A non-finite Ŝ skips the truncation,
# and the step leaves the untruncated factors: `svd` throws for a non-finite matrix.
function _truncate(method, X̂, Ŝ, V̂, r)
    all(isfinite, Ŝ) || return X̂, Ŝ, V̂
    F = svd(Ŝ)
    r₁ = _truncation_rank(method, F.S, r)
    return X̂ * F.U[:, 1:r₁], Matrix(Diagonal(F.S[1:r₁])), V̂ * F.V[:, 1:r₁]
end

function _bug_step(::BUG, X, S, V, w, hx, hv, τ)
    K, L = _klsteps(X, S, V, w, τ)
    X₁, V₁ = _basis(K, hx), _basis(L, hv)
    S₀ = (hx .* (X₁' * X)) * S * (hv .* (V' * V₁))
    return X₁, _sstep(X₁, S₀, V₁, w, τ), V₁
end

function _bug_step(method::AugmentedBUG, X, S, V, w, hx, hv, τ)
    return _truncate(method, _augmented(X, S, V, w, hx, hv, τ)..., size(S, 1))
end

function _bug_step(method::MidpointBUG, X, S, V, w, hx, hv, τ)
    Û, Ŝ½, V̂ = _augmented(X, S, V, w, hx, hv, τ / 2)
    K½, L½ = Û * Ŝ½, V̂ * Ŝ½'
    FV = _kfield!(similar(K½), K½, V̂, w)       # h_v F(Ŷ½) V̂
    FU = _lfield!(similar(L½), Û, L½, w)       # h_x F(Ŷ½)ᵀ Û
    Ū, V̄ = _basis(hcat(Û, τ .* FV), hx), _basis(hcat(V̂, τ .* FU), hv)
    S̄₀ = (hx .* (Ū' * X)) * S * (hv .* (V' * V̄))
    return _truncate(method, Ū, _sstep(Ū, S̄₀, V̄, w, τ), V̄, size(S, 1))
end

_check_rank_bounds(method, r) = nothing
function _check_rank_bounds(m::AugmentedBUG{<:Real}, r)
    m.rmin ≤ r ≤ m.rmax || throw(ArgumentError(
        "the rank $r of the state is not in rmin:rmax = $(m.rmin):$(m.rmax) of the method"))
    return nothing
end

# the step of a run, which `lowrank_step!` of a BUG method converts to the element type of Y
_check_run_step(method, Y, tstep) = nothing
function _check_run_step(::BUGMethod, Y, tstep)
    h = eltype(Y)(tstep)
    isfinite(h) && h > 0 || throw(ArgumentError(
        "the step $tstep is zero or not finite in the element type $(eltype(Y)) of the state"))
    return nothing
end

"""
    lowrank_step!(Y, method, field, t, τ)

Take one step of length `τ` from the time `t` of the dynamical low-rank `method` — [`BUG`](@ref),
[`AugmentedBUG`](@ref) or [`MidpointBUG`](@ref) — in place on the
[`LowRankDistribution`](@ref) `Y` on a [`GridFactorSpace`](@ref), and return `Y`. `field` is the
workspace [`VlasovMethods.LowRankFields1D1V`](@ref) of the same element type `T`, whose K-, L-
and S-fields evaluate the full-order field on the factors, without the `nx × nv` array. The
fields are autonomous, so `t` is accepted and not read. Each K-, L-, S- and Galerkin substep is
one step of the classical RK4 tableau over its interval.

`τ` is converted to `T` before its first use; a non-finite `τ` throws an `ArgumentError` before
any write. `τ = 0` keeps `X S Vᵀ` to round-off, except for `AugmentedBUG(tol; rmin, rmax)`,
whose truncation still applies and changes `X S Vᵀ` by up to `tol` in the weighted Frobenius
norm. The step replaces the factors `X`, `S` and
`V` of `Y`, counts every column of `V` in the rank, the fixed modes of `Y` included, and sets
`nfixed = 0`, because no reference method keeps the fixed modes in the span of `V`.

The entries of `Y` are not checked. An S-step or a Galerkin step that leaves a non-finite value
does not throw: the step skips the truncation and leaves the untruncated factors in `Y`, so the
status that a caller tests is `all(isfinite, Y.S)`.
"""
function lowrank_step!(Y::LowRankDistribution{T, <:GridFactorSpace}, method::BUGMethod,
        field::LowRankFields1D1V{T}, t, τ::Real) where {T}
    h = T(τ)
    isfinite(h) || throw(ArgumentError("the step τ of a BUG integrator is finite; got $h"))
    _check_rank_bounds(method, size(Y.S, 1))
    hx, hv = _weights(Y.space)
    Y.X, Y.S, Y.V = _bug_step(method, Y.X, Y.S, Y.V, field, hx, hv, h)
    Y.nfixed = 0
    return Y
end

"""
    LowRankIntegrator(Y, method, field, tspan, tstep)

The run of a dynamical low-rank `method` from the [`LowRankDistribution`](@ref) `Y` over the
span `tspan = (t₀, t₁)` with the step `tstep`, for any method with a [`lowrank_step!`](@ref) on
the workspace `field`. The fields `Y`, `method`, `field`, `tspan` and `tstep` hold them, the times
promoted to one type. [`run!`](@ref) advances `Y`.

Throws an `ArgumentError` for a step that is not finite and positive, and for a span that is not
finite or has `t₁ < t₀`.
"""
struct LowRankIntegrator{YT <: LowRankDistribution, MT, FT, TT <: Real} <: VlasovMethod
    Y::YT
    method::MT
    field::FT
    tspan::Tuple{TT, TT}
    tstep::TT

    function LowRankIntegrator(Y::YT, method::MT, field::FT, tspan::Tuple{Real, Real},
            tstep::Real) where {YT <: LowRankDistribution, MT, FT}
        t₀, t₁, Δt = promote(tspan..., tstep)
        isfinite(Δt) && Δt > 0 || throw(ArgumentError(
            "LowRankIntegrator takes a finite and positive step; got $tstep"))
        isfinite(t₀) && isfinite(t₁) && t₀ ≤ t₁ || throw(ArgumentError(
            "LowRankIntegrator takes a finite span (t₀, t₁) with t₀ ≤ t₁; got $tspan"))
        new{YT, MT, FT, typeof(Δt)}(Y, method, field, (t₀, t₁), Δt)
    end
end

# the columns of the stored factors: the largest rank that the method can reach
_stored_rank(m::AugmentedBUG{<:Real}, Y) = min(m.rmax, size(Y)...)
_stored_rank(method, Y) = size(Y.S, 1)

"""
    run!(int::LowRankIntegrator, h5file)

Advance the state `Y` of `int` over the times `t₀:tstep:t₁`, the initial time included, in place
with `lowrank_step!(Y, method, field, tₙ, tstep)`, write every state into the HDF5 file
`h5file`, and return `Y`. The file is created, or overwritten, and closed on return, also when a
step throws. A state whose rank is outside `rmin:rmax` of an `AugmentedBUG(tol; rmin, rmax)`
throws an `ArgumentError` before the file is opened, and so does, for [`BUG`](@ref),
[`AugmentedBUG`](@ref) and [`MidpointBUG`](@ref), a `tstep` that is zero or not finite in the
element type of `Y`.

Time `t₀` is the initial state, written before any step, and step `s` is the `s`-th call of
`lowrank_step!`, which writes the time `t₀ + s · tstep`. The file holds the datasets `t` of the
`nt` times, `rank` of the `nt` ranks, and `X` of size `(nx, R, nt)`, `S` of size `(R, R, nt)`
and `V` of size `(nv, R, nt)`, zero beyond the rank of each state, with `R` the rank of `Y` or,
for `AugmentedBUG(tol; rmin, rmax)`, `min(rmax, nx, nv)`.

A step that leaves `all(isfinite, Y.S)` false stops the run without an error: nothing is written
for it, its `rank` and factors stay zero, and the file attribute `status` is
`"failed at step <s>"`. A run that completes sets `status` to `"ok"`.
"""
function run!(int::LowRankIntegrator, h5file)
    (; Y, method, field, tstep) = int
    t₀, t₁ = int.tspan
    times = t₀:tstep:t₁
    nt = length(times)
    nx, nv = size(Y)
    _check_rank_bounds(method, size(Y.S, 1))
    _check_run_step(method, Y, tstep)
    R = _stored_rank(method, Y)
    DT = eltype(Y)

    h5open(h5file, "w") do h5
        h5t = create_dataset(h5, "t", eltype(times), (nt,); chunk = (1,))
        h5r = create_dataset(h5, "rank", Int, (nt,); chunk = (1,))
        h5X = create_dataset(h5, "X", DT, (nx, R, nt); chunk = (nx, R, 1))
        h5S = create_dataset(h5, "S", DT, (R, R, nt); chunk = (R, R, 1))
        h5V = create_dataset(h5, "V", DT, (nv, R, nt); chunk = (nv, R, 1))

        # the factors of Y in the first columns of R, zero beyond
        function store!(n)
            X, S, V = zeros(DT, nx, R), zeros(DT, R, R), zeros(DT, nv, R)
            X[:, axes(Y.X, 2)] .= Y.X
            S[axes(Y.S)...] .= Y.S
            V[:, axes(Y.V, 2)] .= Y.V
            h5t[n], h5r[n] = times[n], size(Y.S, 1)
            h5X[:, :, n], h5S[:, :, n], h5V[:, :, n] = X, S, V
            return nothing
        end

        store!(1)
        status = "ok"
        for s in 1:(nt - 1)
            lowrank_step!(Y, method, field, times[s], tstep)
            if !all(isfinite, Y.S)
                status = "failed at step $s"
                break
            end
            store!(s + 1)
        end
        attrs(h5)["status"] = status
    end

    return Y
end
