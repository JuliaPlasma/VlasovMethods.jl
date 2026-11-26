@doc raw"""
    ImplicitMidpointSubsteps(; projected = true, tol = nothing, max_iterations = 50)

The implicit-midpoint substeps of a [`ProjectorSplitting`](@ref), its keyword `substeps`: a
reference method for the energy. Each K-, S- and L-substep is one implicit-midpoint step of the
field of the bracket at a gradient of the Hamiltonian,

```math
\dot{K} = h_v J(K V^\mathsf{T}, P_V h) V , \qquad
\dot{S} = -h_x h_v X^\mathsf{T} J(X S V^\mathsf{T}, P_X P_V h) V , \qquad
\dot{L} = \left(h_x X^\mathsf{T} J(X L^\mathsf{T}, P_X h)\right)^\mathsf{T} ,
```

with ``J(f, g)`` the full-order field of [`_grid_1d1v_rhs!`](@ref) at the gradient ``g`` in
place of the Hamiltonian of the model, ``h = \phi \otimes 1 + 1 \otimes v^2/2`` the gradient
with ``\phi = G \rho`` of the state of the substep, and the weighted projections
``P_V g = h_v g V V^\mathsf{T}`` onto the span of `V` and ``P_X g = h_x X X^\mathsf{T} g`` onto
the span of `X`. With `projected = true`, the field of each substep is the skew-gradient field
of its own projected gradient, so each implicit-midpoint substep keeps the quadratic Hamiltonian
``H = h_x h_v \sum (v^2/2) f + (h_x/2) \rho^\mathsf{T} G \rho`` and the L² norm
``h_x h_v \sum f^2`` up to the tolerance of its solve. `projected = false` takes the plain
gradient ``h`` in every substep, which does not keep ``H``: the control of the projected
gradient.

Each substep solves ``y_1 = y_0 + τ F((y_0 + y_1)/2)`` with `GeometricSolvers.Newton` and
`GeometricSolvers.GMRES` under `EisenstatWalkerChoice2(T)`, through the analytic product
``dF(y) w = B(w, Q h(y)) + B(y, Q h_\mathrm{lin}(w))`` of the bilinear field ``B`` and the
projection ``Q`` of the affine gradient, with ``h_\mathrm{lin}(w) = G \rho(w) \otimes 1``, from the
explicit Euler predictor ``y_0 + τ F(y_0)``, with
`Options(T; f_abstol = tol * norm(y₀), f_reltol = tol, max_iterations)`. `tol = nothing` is
`32 * eps(T)` in the element type `T` of the state; a step throws an `ArgumentError` for a `tol`
outside `[eps(T), 1)` or not finite in `T`.

Throws an `ArgumentError` for `max_iterations < 1`.
"""
struct ImplicitMidpointSubsteps{TT <: Union{Nothing, Real}}
    projected::Bool
    tol::TT
    max_iterations::Int

    function ImplicitMidpointSubsteps(; projected::Bool = true,
            tol::Union{Nothing, Real} = nothing, max_iterations::Integer = 50)
        max_iterations ≥ 1 || throw(ArgumentError(
            "ImplicitMidpointSubsteps takes max_iterations ≥ 1; got $max_iterations"))
        new{typeof(tol)}(projected, tol, max_iterations)
    end
end

# The field of the bracket for a separable gradient h = a₁ ⊗ b₁ + a₂ ⊗ b₂, on the factors. With
# the offsets d = (dx, dv) of h and e = (ex, ev) of f in -1:1 per direction, Arakawa's stencil is
# [h, f]_I = factor Σ_{d,e} A[d, e] h_{I+d} f_{I+e}. For a term a ⊗ b, the v part of a factor
# enters through the 3 × 3 shifted grams M[ev, dv] = Σ_j V[j + ev, :] b[j + dv] V[j, :]ᵀ and the
# x part through N[ex, dx] = Σ_i X[i + ex, :] a[i + dx] X[i, :]ᵀ, each r × r, and the weights A
# are summed over r × r matrices, so no nx × nv array is formed. The helpers below add the term of
# one a ⊗ b without the factor.

# u[j] = y[j + s], periodic
function _shift!(u, y, s)
    n = length(y)
    for j in 1:n
        u[j] = y[_wrap(j + s, n)]
    end
    return u
end

# G[s + 2, d + 2] = Σ_j Y[j + s, a] u[j + d] Y[j, c], for s and d in -1:1, periodic in j
function _shifted_grams(Y, u)
    T, r = eltype(Y), size(Y, 2)
    ũ = similar(u)
    G = [Matrix{T}(undef, r, r) for _ in 1:3, _ in 1:3]
    for s in -1:1, d in -1:1

        _shifted_gram!(G[s + 2, d + 2], Y, s, _shift!(ũ, u, d))
    end
    return G
end

# K̇[i, c] += Σ_{ex,dx} a[i + dx] Σ_a K[i + ex, a] P[a, c],
# with P = Σ_{ev,dv} A[dx, dv, ex, ev] M[ev, dv]
function _kterm!(K̇, K, V, w, a, b)
    A = w.arakawa.A
    M = _shifted_grams(V, b)
    ã, P, Z = similar(a), similar(M[1]), zero(M[1])
    for ex in -1:1, dx in -1:1

        fill!(P, zero(eltype(P)))
        for ev in -1:1, dv in -1:1

            P .+= A[dx + 2, dv + 2, ex + 2, ev + 2] .* M[ev + 2, dv + 2]
        end
        _shifted_apply!(K̇, K, ex, _shift!(ã, a, dx), P, Z)
    end
    return K̇
end

# Ṡ += Σ_{ex,dx} N[ex, dx]ᵀ S P, with P as in `_kterm!`
function _sterm!(Ṡ, X, S, V, w, a, b)
    A = w.arakawa.A
    N, M = _shifted_grams(X, a), _shifted_grams(V, b)
    P, B = similar(S), similar(S)
    for ex in -1:1, dx in -1:1

        fill!(P, zero(eltype(P)))
        for ev in -1:1, dv in -1:1

            P .+= A[dx + 2, dv + 2, ex + 2, ev + 2] .* M[ev + 2, dv + 2]
        end
        mul!(B, S, P)
        mul!(Ṡ, N[ex + 2, dx + 2]', B, true, true)
    end
    return Ṡ
end

# L̇[j, c] += Σ_{ev,dv} b[j + dv] Σ_a L[j + ev, a] Q[a, c],
# with Q = Σ_{ex,dx} A[dx, dv, ex, ev] N[ex, dx]
function _lterm!(L̇, X, L, w, a, b)
    A = w.arakawa.A
    N = _shifted_grams(X, a)
    b̃, Q, Z = similar(b), similar(N[1]), zero(N[1])
    for ev in -1:1, dv in -1:1

        fill!(Q, zero(eltype(Q)))
        for ex in -1:1, dx in -1:1

            Q .+= A[dx + 2, dv + 2, ex + 2, ev + 2] .* N[ex + 2, dx + 2]
        end
        _shifted_apply!(L̇, L, ev, _shift!(b̃, b, dv), Q, Z)
    end
    return L̇
end

@doc raw"""
    _separable_kfield!(K̇, K, V, w::LowRankFields1D1V, a₁, b₁, a₂, b₂)
    _separable_sfield!(Ṡ, X, S, V, w::LowRankFields1D1V, a₁, b₁, a₂, b₂)
    _separable_lfield!(L̇, X, L, w::LowRankFields1D1V, a₁, b₁, a₂, b₂)

Write the K-, S- or L-field of the bracket of the workspace `w` at the separable gradient
``h = a_1 \otimes b_1 + a_2 \otimes b_2`` into the first argument and return it:

```math
\dot{K} = h_v J(K V^\mathsf{T}, h) V , \qquad
\dot{S} = -h_x h_v X^\mathsf{T} J(X S V^\mathsf{T}, h) V , \qquad
\dot{L} = \left(h_x X^\mathsf{T} J(X L^\mathsf{T}, h)\right)^\mathsf{T} ,
```

with ``J(f, h)`` the field `poisson_apply(w.arakawa, f, h) / (hx * hv)` of the bracket, which
is the right-hand side of [`_grid_1d1v_rhs!`](@ref) at the Hamiltonian ``h``, periodic in both
directions. The S-field carries the minus sign of the backward S-step of the projector
splitting. `a₁` and `a₂` are
vectors of length `nx`, `b₁` and `b₂` of length `nv`, `K` is `nx × r`, `X` is `nx × r`, `V` and
`L` are `nv × r`, and `S` is `r × r`. The field is evaluated on the factors, without the
``n_x \times n_v`` array, at a cost of ``O((n_x + n_v) r^2)`` per term. The sizes, the closure of
`w` and the entries are not checked.
"""
function _separable_kfield!(K̇, K, V, w::LowRankFields1D1V, a₁, b₁, a₂, b₂)
    fill!(K̇, zero(eltype(K̇)))
    _kterm!(K̇, K, V, w, a₁, b₁)
    _kterm!(K̇, K, V, w, a₂, b₂)
    return rmul!(K̇, w.arakawa.factor * w.arakawa.hv)
end

function _separable_sfield!(Ṡ, X, S, V, w::LowRankFields1D1V, a₁, b₁, a₂, b₂)
    fill!(Ṡ, zero(eltype(Ṡ)))
    _sterm!(Ṡ, X, S, V, w, a₁, b₁)
    _sterm!(Ṡ, X, S, V, w, a₂, b₂)
    return rmul!(Ṡ, -w.arakawa.factor * w.arakawa.hx * w.arakawa.hv)
end

function _separable_lfield!(L̇, X, L, w::LowRankFields1D1V, a₁, b₁, a₂, b₂)
    fill!(L̇, zero(eltype(L̇)))
    _lterm!(L̇, X, L, w, a₁, b₁)
    _lterm!(L̇, X, L, w, a₂, b₂)
    return rmul!(L̇, w.arakawa.factor * w.arakawa.hx)
end

# The gradient h = ϕ ⊗ 1 + 1 ⊗ v²/2 of each substep in the separable form, projected or plain:
# P_V h = ϕ ⊗ P_V 1 + 1 ⊗ P_V k, P_X P_V h = P_X ϕ ⊗ P_V 1 + P_X 1 ⊗ P_V k and
# P_X h = P_X ϕ ⊗ 1 + P_X 1 ⊗ k, with P_V b = hv V Vᵀ b, P_X a = hx X Xᵀ a and k = v²/2. The
# part h_lin(g) = G ρ(g) ⊗ 1 that is linear in the state g is projected the same way.
_project(Y, u, h) = h .* (Y * (Y' * u))
_project_v(V, b, w, projected) = projected ? _project(V, b, w.arakawa.hv) : b
_project_x(X, a, w, projected) = projected ? _project(X, a, w.arakawa.hx) : a

function _potential_of(K, V, w::LowRankFields1D1V{T}) where {T}
    _potential!(Vector{T}(undef, w.arakawa.nx), K, V, w)
end
function _potential_of(X, S, V, w::LowRankFields1D1V{T}) where {T}
    _potential!(Vector{T}(undef, w.arakawa.nx), X, S, V, w)
end

@doc raw"""
    _midpoint_kfield!(K̇, K, V, w, projected)
    _midpoint_kproduct!(K̇w, K, W, V, w, projected)
    _midpoint_sfield!(Ṡ, X, S, V, w, projected)
    _midpoint_sproduct!(Ṡw, X, S, W, V, w, projected)
    _midpoint_lfield!(L̇, X, L, w, projected)
    _midpoint_lproduct!(L̇w, X, L, W, w, projected)

The fields ``F`` of the K-, S- and L-substep of [`ImplicitMidpointSubsteps`](@ref), at the
gradient ``Q h`` of the state of the substep, and their products ``dF(y) W``, written into the
first argument, which is returned. ``Q`` is ``P_V``, ``P_X P_V`` and ``P_X`` for the K-, S- and
L-substep with `projected = true`, and the identity with `projected = false`. The fields are
those of [`_separable_kfield!`](@ref) and its S- and L-twins; the potential ``\phi = G \rho`` is
that of `K Vᵀ`, `X S Vᵀ` or `X Lᵀ`, from [`_potential!`](@ref). ``F`` is quadratic in the state,
so ``dF(y) W = B(W, Q h(y)) + B(y, Q h_\mathrm{lin}(W))`` is exact.
"""
function _midpoint_kfield!(K̇, K, V, w::LowRankFields1D1V{T}, projected) where {T}
    b₁ = _project_v(V, ones(T, w.arakawa.nv), w, projected)
    b₂ = _project_v(V, w.kinetic, w, projected)
    ϕ = _potential_of(K, V, w)
    return _separable_kfield!(K̇, K, V, w, ϕ, b₁, ones(T, w.arakawa.nx), b₂)
end

function _midpoint_kproduct!(K̇w, K, W, V, w::LowRankFields1D1V{T}, projected) where {T}
    b₁ = _project_v(V, ones(T, w.arakawa.nv), w, projected)
    b₂ = _project_v(V, w.kinetic, w, projected)
    ϕ, ϕw = _potential_of(K, V, w), _potential_of(W, V, w)
    fill!(K̇w, zero(T))
    _kterm!(K̇w, W, V, w, ϕ, b₁)
    _kterm!(K̇w, W, V, w, ones(T, w.arakawa.nx), b₂)
    _kterm!(K̇w, K, V, w, ϕw, b₁)
    return rmul!(K̇w, w.arakawa.factor * w.arakawa.hv)
end

function _midpoint_sfield!(Ṡ, X, S, V, w::LowRankFields1D1V{T}, projected) where {T}
    a₁ = _project_x(X, _potential_of(X, S, V, w), w, projected)
    a₂ = _project_x(X, ones(T, w.arakawa.nx), w, projected)
    b₁ = _project_v(V, ones(T, w.arakawa.nv), w, projected)
    b₂ = _project_v(V, w.kinetic, w, projected)
    return _separable_sfield!(Ṡ, X, S, V, w, a₁, b₁, a₂, b₂)
end

function _midpoint_sproduct!(Ṡw, X, S, W, V, w::LowRankFields1D1V{T}, projected) where {T}
    a₁ = _project_x(X, _potential_of(X, S, V, w), w, projected)
    a₂ = _project_x(X, ones(T, w.arakawa.nx), w, projected)
    aw = _project_x(X, _potential_of(X, W, V, w), w, projected)
    b₁ = _project_v(V, ones(T, w.arakawa.nv), w, projected)
    b₂ = _project_v(V, w.kinetic, w, projected)
    fill!(Ṡw, zero(T))
    _sterm!(Ṡw, X, W, V, w, a₁, b₁)
    _sterm!(Ṡw, X, W, V, w, a₂, b₂)
    _sterm!(Ṡw, X, S, V, w, aw, b₁)
    return rmul!(Ṡw, -w.arakawa.factor * w.arakawa.hx * w.arakawa.hv)
end

function _midpoint_lfield!(L̇, X, L, w::LowRankFields1D1V{T}, projected) where {T}
    a₁ = _project_x(X, _potential_of(X, L, w), w, projected)
    a₂ = _project_x(X, ones(T, w.arakawa.nx), w, projected)
    return _separable_lfield!(L̇, X, L, w, a₁, ones(T, w.arakawa.nv), a₂, w.kinetic)
end

function _midpoint_lproduct!(L̇w, X, L, W, w::LowRankFields1D1V{T}, projected) where {T}
    a₁ = _project_x(X, _potential_of(X, L, w), w, projected)
    a₂ = _project_x(X, ones(T, w.arakawa.nx), w, projected)
    aw = _project_x(X, _potential_of(X, W, w), w, projected)
    b₁ = ones(T, w.arakawa.nv)
    fill!(L̇w, zero(T))
    _lterm!(L̇w, X, W, w, a₁, b₁)
    _lterm!(L̇w, X, W, w, a₂, w.kinetic)
    _lterm!(L̇w, X, L, w, aw, b₁)
    return rmul!(L̇w, w.arakawa.factor * w.arakawa.hx)
end

# The residual R(x) = x − x₀ − τ F((y₀ + y₁)/2) of the implicit-midpoint step from y₀, with x the
# vector of y₁ and x₀ that of y₀, and its product dR(x) v = v − (τ/2) dF((y₀ + y₁)/2) v, in the
# forms `residual!(R, x, p)` and `jvp!(Jv, x, v, p)` of GeometricSolvers
function _midpoint_residual(F!, dF!, y₀, τ)
    sz = size(y₀)
    x₀ = vec(copy(y₀))
    ym, ẏ, ẏw = similar(y₀), similar(y₀), similar(y₀)
    function residual!(R, x, p)
        ym .= (y₀ .+ reshape(x, sz)) ./ 2
        F!(ẏ, ym)
        R .= x .- x₀ .- τ .* vec(ẏ)
        return R
    end
    function jvp!(Jv, x, v, p)
        ym .= (y₀ .+ reshape(x, sz)) ./ 2
        dF!(ẏw, ym, reshape(v, sz))
        Jv .= v .- (τ / 2) .* vec(ẏw)
        return Jv
    end
    return residual!, jvp!
end

# One implicit-midpoint step y₁ = y₀ + τ F((y₀ + y₁)/2) of the substep `name`, with `F!(ẏ, y)`
# and `dF!(ẏw, y, w)`, solved by Newton–GMRES from the explicit Euler predictor. A solve that
# does not return SUCCESS throws, before the state of the step is written.
function _midpoint_solve(F!, dF!, y₀::Matrix{T}, τ::T, tol::T, max_iterations,
        name) where {T}
    sz = size(y₀)
    x₀ = vec(copy(y₀))
    residual!, jvp! = _midpoint_residual(F!, dF!, y₀, τ)
    GS = GeometricSolvers
    method = GS.Newton(; linsolve = GS.GMRES(; forcing = GS.EisenstatWalkerChoice2(T)),
        jacobian = GS.JVPFunction(jvp!))
    y₀norm = norm(y₀)
    options = GS.Options(T; f_abstol = tol * y₀norm, f_reltol = tol, max_iterations)
    predictor = x₀ .+ τ .* vec(F!(similar(y₀), y₀))
    x, status = GS.solve(
        GS.NonlinearProblem(residual!), predictor, method, nothing; options)
    status.code === GS.SUCCESS || throw(ErrorException(
        "the $name-substep solve of ImplicitMidpointSubsteps returned $(status.code), with " *
        "‖R‖₂/‖y₀‖₂ = $(status.fnorm / y₀norm)"))
    return reshape(x, sz)
end

# the K-substep from K = X S, then X, S = orthonormalize_x(space, K)
function _midpoint_kstep(space, K, V, w, τ, tol, sub)
    F!(K̇, K) = _midpoint_kfield!(K̇, K, V, w, sub.projected)
    dF!(K̇w, K, W) = _midpoint_kproduct!(K̇w, K, W, V, w, sub.projected)
    return orthonormalize_x(space, _midpoint_solve(
        F!, dF!, K, τ, tol, sub.max_iterations, "K"))
end

# the backward S-substep, whose field carries the minus sign, over +τ
function _midpoint_sstep(X, S, V, w, τ, tol, sub)
    F!(Ṡ, S) = _midpoint_sfield!(Ṡ, X, S, V, w, sub.projected)
    dF!(Ṡw, S, W) = _midpoint_sproduct!(Ṡw, X, S, W, V, w, sub.projected)
    return _midpoint_solve(F!, dF!, S, τ, tol, sub.max_iterations, "S")
end

# the L-substep from L = V Sᵀ, then V, R = orthonormalize_v(space, L) and S = Rᵀ
function _midpoint_lstep(space, X, S, V, w, τ, tol, sub)
    F!(L̇, L) = _midpoint_lfield!(L̇, X, L, w, sub.projected)
    dF!(L̇w, L, W) = _midpoint_lproduct!(L̇w, X, L, W, w, sub.projected)
    L = _midpoint_solve(F!, dF!, V * S', τ, tol, sub.max_iterations, "L")
    V₁, R = orthonormalize_v(space, L)
    return permutedims(R), V₁
end

# τ and the solve tolerance in the element type T of Y, after the checks that run before any
# solve: τ finite and nonzero, tol in [eps(T), 1), finite factors, and the periodic closure
function _check_midpoint(Y, sub, w, τ)
    T = eltype(Y)
    τT = _finite_step(τ, T)
    iszero(τT) && throw(ArgumentError(
        "ImplicitMidpointSubsteps takes a nonzero step τ; got τ = $τ"))
    # the range refuses NaN and ±Inf too
    tol = sub.tol === nothing ? 32 * eps(T) : T(sub.tol)
    eps(T) ≤ tol < 1 || throw(ArgumentError(
        "ImplicitMidpointSubsteps takes a tol in [eps($T), 1) in the element type $T of the " *
        "state; got tol = $(sub.tol)"))
    for (name, A) in ((:X, Y.X), (:S, Y.S), (:V, Y.V))
        all(isfinite, A) || throw(ArgumentError(
            "ImplicitMidpointSubsteps takes a state with finite entries; the factor $name has " *
            "a non-finite entry"))
    end
    _bounded(w) && throw(ArgumentError(
        "ImplicitMidpointSubsteps takes a field of the PeriodicVelocityClosure; got the " *
        "BoundedVelocityClosure"))
    return τT, tol
end

@doc raw"""
    lowrank_step!(Y::LowRankDistribution,
                  method::ProjectorSplitting{C, <:ImplicitMidpointSubsteps},
                  field::LowRankFields1D1V, t, τ)

Advance the low-rank state `Y` by one step of length `τ` of the projector splitting with the
[`ImplicitMidpointSubsteps`](@ref) `method.substeps`, in place, and return `Y`. Each substep is
one implicit-midpoint step of its self-consistent field at its own state; after a K-substep,
`X, S = orthonormalize_x(Y.space, K)`, and after an L-substep,
`V, R = orthonormalize_v(Y.space, L)` and `S = Rᵀ`, as for the RK4 substeps of
[`ProjectorSplitting`](@ref).

  - `C = :lie` runs K(τ) S(τ) L(τ): the K-substep from ``K = X S``, the backward S-substep and the
    L-substep from ``L = V S^\mathsf{T}``, each over the step `τ`.
  - `C = :strang` runs the symmetric K(τ/2) S(τ/2) L(τ) S(τ/2) K(τ/2), so a step of `-τ` after a
    step of `τ` returns the state up to the tolerance of the solves.

The field is autonomous, so `t` is not read. `τ` is converted to the element type `T` of `Y`
and may be negative. Before any solve, the step throws the `ArgumentError` of the generic
[`lowrank_step!`](@ref) for a `τ` that is not finite in `T`, and an `ArgumentError` for `τ = 0`,
for a non-finite entry of `X`, `S` or `V`, which the message names, for a `tol` outside
`[eps(T), 1)`, and for a `field` with the [`BoundedVelocityClosure`](@ref). A solve that does not
return `SUCCESS` throws an `ErrorException` whose message names the substep, `K`, `S` or `L`, the
return code and ``\|R\|_2 / \|y_0\|_2``. Each throw leaves `Y` as it was. The step forms no
``n_x \times n_v`` array, and sets `Y.nfixed = 0`. A zero `S` stays zero.
"""
function lowrank_step!(Y::LowRankDistribution,
        method::ProjectorSplitting{:lie, <:ImplicitMidpointSubsteps},
        field::LowRankFields1D1V, t, τ)
    sub = method.substeps
    τ, tol = _check_midpoint(Y, sub, field, τ)
    X, S₁ = _midpoint_kstep(Y.space, Y.X * Y.S, Y.V, field, τ, tol, sub)
    S₂ = _midpoint_sstep(X, S₁, Y.V, field, τ, tol, sub)
    S₃, V = _midpoint_lstep(Y.space, X, S₂, Y.V, field, τ, tol, sub)
    Y.X, Y.S, Y.V, Y.nfixed = X, S₃, V, 0
    return Y
end

function lowrank_step!(Y::LowRankDistribution,
        method::ProjectorSplitting{:strang, <:ImplicitMidpointSubsteps},
        field::LowRankFields1D1V, t, τ)
    sub = method.substeps
    τ, tol = _check_midpoint(Y, sub, field, τ)
    τ½ = τ / 2
    X₁, S₁ = _midpoint_kstep(Y.space, Y.X * Y.S, Y.V, field, τ½, tol, sub)
    S₂ = _midpoint_sstep(X₁, S₁, Y.V, field, τ½, tol, sub)
    S₃, V₁ = _midpoint_lstep(Y.space, X₁, S₂, Y.V, field, τ, tol, sub)
    S₄ = _midpoint_sstep(X₁, S₃, V₁, field, τ½, tol, sub)
    X₂, S₅ = _midpoint_kstep(Y.space, X₁ * S₄, V₁, field, τ½, tol, sub)
    Y.X, Y.S, Y.V, Y.nfixed = X₂, S₅, V₁, 0
    return Y
end
