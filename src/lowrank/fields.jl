@doc raw"""
    LowRankFields1D1V(arakawa::Arakawa{T}, G::AbstractMatrix{T}, v::AbstractVector;
                      closure = PeriodicVelocityClosure())
    LowRankFields1D1V(space::GridFactorSpace, G; closure = PeriodicVelocityClosure())
    LowRankFields1D1V(dist::LowRankDistribution, G; closure = PeriodicVelocityClosure())

The workspace of the K-, L- and S-fields of 1D1V Vlasov–Poisson on an `nx × nv` grid, in the
element type `T`: the bracket `arakawa`, the ``n_x \times n_x`` Green's matrix `G` of the
x-grid, as [`_green_matrix`](@ref) builds it, the closure `closure` of the velocity grid, and the
kinetic energy `kinetic` of the velocity nodes `v`, in `T`. With the
[`PeriodicVelocityClosure`](@ref), `kinetic = collect(v) .^ 2 ./ 2`, and the fields are those of
the full-order right-hand side [`_grid_1d1v_rhs!`](@ref) on the same bracket, `G` and nodes.
With the [`BoundedVelocityClosure`](@ref), `kinetic` is zero at the end nodes ``j = 1`` and
``j = n_v``, and the fields are those of [`_grid_1d1v_bounded_rhs!`](@ref).

The first method throws a `DimensionMismatch` when `size(G) != (nx, nx)` or
`length(v) != nv`, and a `TypeError` for a `closure` other than
[`PeriodicVelocityClosure`](@ref) and [`BoundedVelocityClosure`](@ref). The second builds the bracket
`Arakawa(nx, nv, step(space.x), step(space.v))` of the grid of `space`, with
`(nx, nv) = size(space)`, and takes the nodes `space.v`. The third is the second on the space of
`dist`, so it throws a `MethodError` for a distribution on a space other than a
[`GridFactorSpace`](@ref), because the stencil and its weights belong to the grid.

The functions [`_kfield!`](@ref), [`_lfield!`](@ref), [`_sfield!`](@ref) and
[`_potential!`](@ref) take the workspace.
"""
struct LowRankFields1D1V{T, GT <: AbstractMatrix{T}, CT}
    arakawa::Arakawa{T}
    G::GT
    kinetic::Vector{T}
    closure::CT

    function LowRankFields1D1V(
            arakawa::Arakawa{T}, G::AbstractMatrix{T}, v::AbstractVector;
            closure = PeriodicVelocityClosure()) where {T}
        nx, nv = arakawa.nx, arakawa.nv
        size(G) == (nx, nx) || throw(DimensionMismatch(
            "LowRankFields1D1V on a $nx × $nv grid takes G of size ($nx, $nx); got $(size(G))"))
        length(v) == nv || throw(DimensionMismatch(
            "LowRankFields1D1V on a $nx × $nv grid takes $nv velocity nodes; got $(length(v))"))
        closure isa VelocityClosure ||
            throw(TypeError(:LowRankFields1D1V, "keyword closure", VelocityClosure, closure))
        _check_closure(:LowRankFields1D1V, closure)
        kinetic = Vector{T}(collect(v) .^ 2 ./ 2)
        if closure isa BoundedVelocityClosure
            kinetic[1] = kinetic[nv] = zero(T)
        end
        new{T, typeof(G), typeof(closure)}(arakawa, G, kinetic, closure)
    end
end

function LowRankFields1D1V(space::GridFactorSpace, G::AbstractMatrix;
        closure = PeriodicVelocityClosure())
    nx, nv = size(space)
    return LowRankFields1D1V(Arakawa(nx, nv, step(space.x), step(space.v)), G, space.v;
        closure)
end

function LowRankFields1D1V(dist::LowRankDistribution, G::AbstractMatrix;
        closure = PeriodicVelocityClosure())
    LowRankFields1D1V(dist.space, G; closure)
end

# The closure types are defined with the full-order model, which is included after this file,
# so the fields branch on the closure at run time; the branch is on a type parameter of the
# workspace, and the compiler folds it.
_bounded(w::LowRankFields1D1V) = w.closure isa BoundedVelocityClosure

# The factored form of Arakawa's Jacobian. With the offsets d (of h) and e (of f) in -1:1 per
# direction, [h, f]_I = factor Σ_{d,e} A[d, e] h_{I+d} f_{I+e}, periodic in both directions. For
# h = ϕ ⊗ 1 + 1 ⊗ k, the term of ϕ sums A over the v offset of d and the term of k over its x
# offset, so for each offset e of f the stencil sums ϕ̃_e[i] = Σ_dx c_dx ϕ[i + dx] and
# k̃_e[j] = Σ_dv c_dv k[j + dv] carry all of h. A factor of f shifted by e then enters through the
# r × r matrices of `_shifted_gram!`, and no nx × nv array is formed.

# the periodic neighbour i of 1:n, for an i in 0:(n + 1)
_wrap(i, n) = i < 1 ? i + n : (i > n ? i - n : i)

# the coefficients of ϕ[i + dx] and of k[j + dv], dx and dv in -1:1, for the offset (ex, ev) of f
function _ϕ_coefficients(A, ex, ev)
    ntuple(m -> sum(A[m, dv + 2, ex + 2, ev + 2] for dv in -1:1), 3)
end
function _k_coefficients(A, ex, ev)
    ntuple(m -> sum(A[dx + 2, m, ex + 2, ev + 2] for dx in -1:1), 3)
end

# u[i] = c[1] y[i - 1] + c[2] y[i] + c[3] y[i + 1], periodic
function _stencil_sum!(u, y, c)
    n = length(y)
    for i in 1:n
        u[i] = c[1] * y[_wrap(i - 1, n)] + c[2] * y[i] + c[3] * y[_wrap(i + 1, n)]
    end
    return u
end

# M[a, b] = Σ_i Y[i + s, a] Y[i, b], periodic in i
function _shifted_gram!(M, Y, s)
    n = size(Y, 1)
    for b in axes(M, 2), a in axes(M, 1)

        m = zero(eltype(M))
        for i in 1:n
            m += Y[_wrap(i + s, n), a] * Y[i, b]
        end
        M[a, b] = m
    end
    return M
end

# M[a, b] = Σ_i Y[i + s, a] u[i] Y[i, b], periodic in i
function _shifted_gram!(M, Y, s, u)
    n = size(Y, 1)
    for b in axes(M, 2), a in axes(M, 1)

        m = zero(eltype(M))
        for i in 1:n
            m += Y[_wrap(i + s, n), a] * u[i] * Y[i, b]
        end
        M[a, b] = m
    end
    return M
end

# out[i, b] += Σ_a Y[i + s, a] (p[i] P[a, b] + Q[a, b]), periodic in i
function _shifted_apply!(out, Y, s, p, P, Q)
    n = size(Y, 1)
    for b in axes(out, 2), a in axes(Y, 2)

        pab, qab = P[a, b], Q[a, b]
        for i in 1:n
            out[i, b] += Y[_wrap(i + s, n), a] * (p[i] * pab + qab)
        end
    end
    return out
end

# The bounded closure reads f and h as zero on and beyond the end nodes, so the velocity sums run
# over the unknowns j = 2 … nv-1 only, without a wrap, and Ph = ϕ ⊗ χ + 1 ⊗ (χ ∘ v²/2) has rank 2
# too. The term of ϕ then no longer sums A over the v offset dv of h alone: for each dv, the
# coefficient of ϕ[i + dx] is A[dx, dv, e] and the weight of the velocity sum is χ[j + dv], which
# restricts the sum to the j whose neighbour j + dv is an unknown. The term of k keeps the stencil
# sum of the cut kinetic energy. A row of a factor outside the unknowns is never read.

# the j in 2:(n - 1) with j + s and j + t in 2:(n - 1) too
_unknowns(n, s, t = 0) = max(2, 2 - s, 2 - t):min(n - 1, n - 1 - s, n - 1 - t)

# the coefficients of ϕ[i + dx], dx in -1:1, for the offset dv of h and the offset (ex, ev) of f
_ϕ_coefficients(A, ex, ev, dv) = ntuple(m -> A[m, dv + 2, ex + 2, ev + 2], 3)

# M[a, b] = Σ_{j ∈ js} Y[j + s, a] Y[j, b]
function _cut_gram!(M, Y, s, js)
    for b in axes(M, 2), a in axes(M, 1)

        m = zero(eltype(M))
        for j in js
            m += Y[j + s, a] * Y[j, b]
        end
        M[a, b] = m
    end
    return M
end

# M[a, b] = Σ_{j ∈ js} Y[j + s, a] u[j] Y[j, b]
function _cut_gram!(M, Y, s, js, u)
    for b in axes(M, 2), a in axes(M, 1)

        m = zero(eltype(M))
        for j in js
            m += Y[j + s, a] * u[j] * Y[j, b]
        end
        M[a, b] = m
    end
    return M
end

# out[j, b] += Σ_a Y[j + s, a] P[a, b], for j ∈ js
function _cut_apply!(out, Y, s, js, P)
    for b in axes(out, 2), a in axes(Y, 2)

        pab = P[a, b]
        for j in js
            out[j, b] += Y[j + s, a] * pab
        end
    end
    return out
end

# out[j, b] += Σ_a Y[j + s, a] u[j] P[a, b], for j ∈ js
function _cut_apply!(out, Y, s, js, u, P)
    for b in axes(out, 2), a in axes(Y, 2)

        pab = P[a, b]
        for j in js
            out[j, b] += Y[j + s, a] * u[j] * pab
        end
    end
    return out
end

# the sum of the rows of V that the density reads: all rows, or the unknowns of the bounded closure
function _velocity_sum(V, w)
    _bounded(w) || return vec(sum(V; dims = 1))
    return vec(sum(view(V, 2:(w.arakawa.nv - 1), :); dims = 1))
end

function _check_alias(name, outname, out, inputs::NamedTuple)
    for (label, x) in pairs(inputs)
        Base.mightalias(out, x) && throw(ArgumentError(
            "$name: the output $outname shares memory with the input $label"))
    end
    return nothing
end

_ϕ_wanted(ϕ, nx) = ϕ === nothing ? "" : ", and ϕ of length $nx"
_ϕ_got(ϕ) = ϕ === nothing ? "" : ", and ϕ of length $(length(ϕ))"

function _check_kfield(K̇, K, V, w, ϕ = nothing)
    nx, nv = w.arakawa.nx, w.arakawa.nv
    r = size(K, 2)
    size(K̇) == size(K) == (nx, r) && size(V) == (nv, r) &&
    (ϕ === nothing || length(ϕ) == nx) ||
        throw(DimensionMismatch(
            "_kfield! on a $nx × $nv grid takes K̇ and K of $nx rows and V of $nv rows, " *
            "with $r columns each$(_ϕ_wanted(ϕ, nx)); got K̇ of size $(size(K̇)), " *
            "K of size $(size(K)) and V of size $(size(V))$(_ϕ_got(ϕ))"))
    return nothing
end

function _check_lfield(L̇, X, L, w, ϕ = nothing)
    nx, nv = w.arakawa.nx, w.arakawa.nv
    r = size(L, 2)
    size(L̇) == size(L) == (nv, r) && size(X) == (nx, r) &&
    (ϕ === nothing || length(ϕ) == nx) ||
        throw(DimensionMismatch(
            "_lfield! on a $nx × $nv grid takes L̇ and L of $nv rows and X of $nx rows, " *
            "with $r columns each$(_ϕ_wanted(ϕ, nx)); got L̇ of size $(size(L̇)), " *
            "X of size $(size(X)) and L of size $(size(L))$(_ϕ_got(ϕ))"))
    return nothing
end

function _check_sfield(Ṡ, X, S, V, w, ϕ = nothing)
    nx, nv = w.arakawa.nx, w.arakawa.nv
    r = size(X, 2)
    size(X) == (nx, r) && size(V) == (nv, r) && size(Ṡ) == size(S) == (r, r) &&
    (ϕ === nothing || length(ϕ) == nx) ||
        throw(DimensionMismatch(
            "_sfield! on a $nx × $nv grid takes X of $nx rows, V of $nv rows, and Ṡ and S " *
            "of size r × r, with r the columns of X and V$(_ϕ_wanted(ϕ, nx)); got Ṡ of size " *
            "$(size(Ṡ)), X of size $(size(X)), S of size $(size(S)) and V of size " *
            "$(size(V))$(_ϕ_got(ϕ))"))
    return nothing
end

@doc raw"""
    _kfield!(K̇, K, V, w::LowRankFields1D1V)
    _kfield!(K̇, K, V, w::LowRankFields1D1V, ϕ)

Write the K-field of 1D1V Vlasov–Poisson into `K̇` and return it:

```math
\dot{K} = h_v \, F(K V^\mathsf{T}) \, V , \qquad F(f) = [h, f] ,
\qquad h = \phi \otimes 1 + 1 \otimes \tfrac{v^2}{2} ,
```

where ``F`` is the right-hand side of [`_grid_1d1v_rhs!`](@ref) on the bracket, `G` and the
velocity nodes of the workspace `w`, `K` is `nx × r`, `V` is `nv × r`, and `K̇` is `nx × r`.
Without `ϕ` the potential is that of ``f = K V^\mathsf{T}``, from [`_potential!`](@ref); with
`ϕ`, a `Vector` of length `nx`, the potential is the given `ϕ`, frozen. The factors are
`Matrix`es of the element type `T` of `w`; `V` need not be orthonormal.

On a workspace with the [`BoundedVelocityClosure`](@ref), ``F`` is the bounded field
``F_b(f) = \chi \circ [P h, \chi \circ f]`` of [`_grid_1d1v_bounded_rhs!`](@ref), with
``P h = \chi \circ (\phi \otimes 1 + 1 \otimes v^2/2)`` also at a frozen `ϕ`: rows ``1`` and
``n_v`` of `V` are read as zero, and the velocity sums run over the unknowns.

The field is evaluated on the factors, without the ``n_x \times n_v`` array: ``h`` has rank 2,
so each term of Arakawa's stencil projects onto `V` through ``r \times r`` matrices of shifted
columns, at a cost of ``O((n_x + n_v) r^2)`` beyond the product with `G`.

Throws a `DimensionMismatch` before any write when the sizes do not fit the grid of `w` and the
rank `r = size(K, 2)`, and an `ArgumentError` when `K̇` shares memory with `K`, `V` or `ϕ`. At
`r = 0` it returns the empty `K̇`. A non-finite entry propagates to the result, unchecked.
"""
function _kfield!(K̇::Matrix{T}, K::Matrix{T}, V::Matrix{T},
        w::LowRankFields1D1V{T}) where {T}
    _check_kfield(K̇, K, V, w)
    ϕ = _potential!(Vector{T}(undef, w.arakawa.nx), K, V, w)
    return _kfield!(K̇, K, V, w, ϕ)
end

function _kfield!(K̇::Matrix{T}, K::Matrix{T}, V::Matrix{T}, w::LowRankFields1D1V{T},
        ϕ::Vector{T}) where {T}
    _check_kfield(K̇, K, V, w, ϕ)
    _check_alias("_kfield!", "K̇", K̇, (K = K, V = V, ϕ = ϕ))
    _bounded(w) && return _bounded_kfield!(K̇, K, V, w, ϕ)
    (; arakawa, kinetic) = w
    A = arakawa.A
    r = size(K, 2)
    ϕ̃, k̃ = similar(ϕ), similar(kinetic)
    M0, Mk = Matrix{T}(undef, r, r), Matrix{T}(undef, r, r)

    fill!(K̇, zero(T))
    for ex in -1:1, ev in -1:1

        _stencil_sum!(ϕ̃, ϕ, _ϕ_coefficients(A, ex, ev))
        _stencil_sum!(k̃, kinetic, _k_coefficients(A, ex, ev))
        _shifted_gram!(M0, V, ev)
        _shifted_gram!(Mk, V, ev, k̃)
        _shifted_apply!(K̇, K, ex, ϕ̃, M0, Mk)
    end
    return rmul!(K̇, arakawa.factor * arakawa.hv)
end

@doc raw"""
    _lfield!(L̇, X, L, w::LowRankFields1D1V)
    _lfield!(L̇, X, L, w::LowRankFields1D1V, ϕ)

Write the L-field of 1D1V Vlasov–Poisson into `L̇` and return it:

```math
\dot{L} = h_x \, F(X L^\mathsf{T})^\mathsf{T} X ,
```

with ``F`` as in [`_kfield!`](@ref), `X` of size `nx × r`, and `L` and `L̇` of size `nv × r`.
Without `ϕ` the potential is that of ``f = X L^\mathsf{T}``, from
[`_potential!`](@ref)`(ϕ, X, L, w)`; with `ϕ`, the potential is the given `ϕ`, frozen. `X` need
not be orthonormal. The evaluation, the closure, the edges and the errors are those of
[`_kfield!`](@ref), with `L̇` in place of `K̇` and the rank `r = size(L, 2)`; on a workspace with
the [`BoundedVelocityClosure`](@ref), rows ``1`` and ``n_v`` of `L` are read as zero, and rows
``1`` and ``n_v`` of `L̇` are zero.
"""
function _lfield!(L̇::Matrix{T}, X::Matrix{T}, L::Matrix{T},
        w::LowRankFields1D1V{T}) where {T}
    _check_lfield(L̇, X, L, w)
    ϕ = _potential!(Vector{T}(undef, w.arakawa.nx), X, L, w)
    return _lfield!(L̇, X, L, w, ϕ)
end

function _lfield!(L̇::Matrix{T}, X::Matrix{T}, L::Matrix{T}, w::LowRankFields1D1V{T},
        ϕ::Vector{T}) where {T}
    _check_lfield(L̇, X, L, w, ϕ)
    _check_alias("_lfield!", "L̇", L̇, (X = X, L = L, ϕ = ϕ))
    _bounded(w) && return _bounded_lfield!(L̇, X, L, w, ϕ)
    (; arakawa, kinetic) = w
    A = arakawa.A
    r = size(L, 2)
    ϕ̃, k̃ = similar(ϕ), similar(kinetic)
    Nϕ, N0 = Matrix{T}(undef, r, r), Matrix{T}(undef, r, r)

    fill!(L̇, zero(T))
    for ex in -1:1, ev in -1:1

        _stencil_sum!(ϕ̃, ϕ, _ϕ_coefficients(A, ex, ev))
        _stencil_sum!(k̃, kinetic, _k_coefficients(A, ex, ev))
        _shifted_gram!(Nϕ, X, ex, ϕ̃)
        _shifted_gram!(N0, X, ex)
        _shifted_apply!(L̇, L, ev, k̃, N0, Nϕ)
    end
    return rmul!(L̇, arakawa.factor * arakawa.hx)
end

@doc raw"""
    _sfield!(Ṡ, X, S, V, w::LowRankFields1D1V)
    _sfield!(Ṡ, X, S, V, w::LowRankFields1D1V, ϕ)

Write the S-field of 1D1V Vlasov–Poisson into `Ṡ` and return it:

```math
\dot{S} = h_x h_v \, X^\mathsf{T} F(X S V^\mathsf{T}) \, V ,
```

with ``F`` as in [`_kfield!`](@ref), `X` of size `nx × r`, `V` of size `nv × r`, and `S` and
`Ṡ` of size `r × r`. Without `ϕ` the potential is that of ``f = X S V^\mathsf{T}``, from
[`_potential!`](@ref)`(ϕ, X, S, V, w)`; with `ϕ`, the potential is the given `ϕ`, frozen. `X`
and `V` need not be orthonormal. At `r = 1` the field vanishes, up to round-off, by the
antisymmetry of Arakawa's form. The evaluation, the closure, the edges and the errors are those
of [`_kfield!`](@ref), with `Ṡ` in place of `K̇` and the rank `r = size(X, 2)`.
"""
function _sfield!(Ṡ::Matrix{T}, X::Matrix{T}, S::Matrix{T}, V::Matrix{T},
        w::LowRankFields1D1V{T}) where {T}
    _check_sfield(Ṡ, X, S, V, w)
    ϕ = _potential!(Vector{T}(undef, w.arakawa.nx), X, S, V, w)
    return _sfield!(Ṡ, X, S, V, w, ϕ)
end

function _sfield!(Ṡ::Matrix{T}, X::Matrix{T}, S::Matrix{T}, V::Matrix{T},
        w::LowRankFields1D1V{T}, ϕ::Vector{T}) where {T}
    _check_sfield(Ṡ, X, S, V, w, ϕ)
    _check_alias("_sfield!", "Ṡ", Ṡ, (X = X, S = S, V = V, ϕ = ϕ))
    _bounded(w) && return _bounded_sfield!(Ṡ, X, S, V, w, ϕ)
    (; arakawa, kinetic) = w
    A = arakawa.A
    ϕ̃, k̃ = similar(ϕ), similar(kinetic)
    Nϕ, N0, M0, Mk = similar(S), similar(S), similar(S), similar(S)
    B = similar(S)

    fill!(Ṡ, zero(T))
    for ex in -1:1, ev in -1:1

        _stencil_sum!(ϕ̃, ϕ, _ϕ_coefficients(A, ex, ev))
        _stencil_sum!(k̃, kinetic, _k_coefficients(A, ex, ev))
        _shifted_gram!(Nϕ, X, ex, ϕ̃)
        _shifted_gram!(N0, X, ex)
        _shifted_gram!(M0, V, ev)
        _shifted_gram!(Mk, V, ev, k̃)
        mul!(B, S, M0)
        mul!(Ṡ, Nϕ', B, true, true)
        mul!(B, S, Mk)
        mul!(Ṡ, N0', B, true, true)
    end
    return rmul!(Ṡ, arakawa.factor * arakawa.hx * arakawa.hv)
end

# The K-, L- and S-fields of the bounded closure, after the checks of `_kfield!`, `_lfield!` and
# `_sfield!`: the sums of those fields with the velocity sums over the unknowns, the weight
# χ[j + dv] of the term of ϕ for each offset dv of h, and the cut kinetic energy of `w`.
function _bounded_kfield!(K̇, K, V, w, ϕ)
    (; arakawa, kinetic) = w
    A, nv = arakawa.A, arakawa.nv
    r = size(K, 2)
    T = eltype(K̇)
    ϕ̃, k̃ = similar(ϕ), similar(kinetic)
    M, Mk, Z = Matrix{T}(undef, r, r), Matrix{T}(undef, r, r), zeros(T, r, r)

    fill!(K̇, zero(T))
    for ex in -1:1, ev in -1:1

        _stencil_sum!(k̃, kinetic, _k_coefficients(A, ex, ev))
        _cut_gram!(Mk, V, ev, _unknowns(nv, ev), k̃)
        for dv in -1:1
            _stencil_sum!(ϕ̃, ϕ, _ϕ_coefficients(A, ex, ev, dv))
            _cut_gram!(M, V, ev, _unknowns(nv, ev, dv))
            _shifted_apply!(K̇, K, ex, ϕ̃, M, dv == 0 ? Mk : Z)
        end
    end
    return rmul!(K̇, arakawa.factor * arakawa.hv)
end

function _bounded_lfield!(L̇, X, L, w, ϕ)
    (; arakawa, kinetic) = w
    A, nv = arakawa.A, arakawa.nv
    r = size(L, 2)
    T = eltype(L̇)
    ϕ̃, k̃ = similar(ϕ), similar(kinetic)
    Nϕ, N0 = Matrix{T}(undef, r, r), Matrix{T}(undef, r, r)

    fill!(L̇, zero(T))
    for ex in -1:1, ev in -1:1

        _stencil_sum!(k̃, kinetic, _k_coefficients(A, ex, ev))
        _shifted_gram!(N0, X, ex)
        _cut_apply!(L̇, L, ev, _unknowns(nv, ev), k̃, N0)
        for dv in -1:1
            _stencil_sum!(ϕ̃, ϕ, _ϕ_coefficients(A, ex, ev, dv))
            _shifted_gram!(Nϕ, X, ex, ϕ̃)
            _cut_apply!(L̇, L, ev, _unknowns(nv, ev, dv), Nϕ)
        end
    end
    return rmul!(L̇, arakawa.factor * arakawa.hx)
end

function _bounded_sfield!(Ṡ, X, S, V, w, ϕ)
    (; arakawa, kinetic) = w
    A, nv = arakawa.A, arakawa.nv
    ϕ̃, k̃ = similar(ϕ), similar(kinetic)
    Nϕ, N0, M, Mk, B = similar(S), similar(S), similar(S), similar(S), similar(S)

    fill!(Ṡ, zero(eltype(Ṡ)))
    for ex in -1:1, ev in -1:1

        _stencil_sum!(k̃, kinetic, _k_coefficients(A, ex, ev))
        _shifted_gram!(N0, X, ex)
        _cut_gram!(Mk, V, ev, _unknowns(nv, ev), k̃)
        mul!(B, S, Mk)
        mul!(Ṡ, N0', B, true, true)
        for dv in -1:1
            _stencil_sum!(ϕ̃, ϕ, _ϕ_coefficients(A, ex, ev, dv))
            _shifted_gram!(Nϕ, X, ex, ϕ̃)
            _cut_gram!(M, V, ev, _unknowns(nv, ev, dv))
            mul!(B, S, M)
            mul!(Ṡ, Nϕ', B, true, true)
        end
    end
    return rmul!(Ṡ, arakawa.factor * arakawa.hx * arakawa.hv)
end

@doc raw"""
    _potential!(ϕ, K, V, w::LowRankFields1D1V)
    _potential!(ϕ, X, S, V, w::LowRankFields1D1V)

Write the potential ``\phi = G \rho`` of factored data into the `Vector` `ϕ` of length `nx` and
return it, with `G` of the workspace `w` and the density

```math
\rho = h_v K (V^\mathsf{T} 1) \quad \text{of } f = K V^\mathsf{T} , \qquad
\rho = h_v X S (V^\mathsf{T} 1) \quad \text{of } f = X S V^\mathsf{T} ,
```

the density ``\rho_i = h_v \sum_j f_{ij}`` of [`_grid_1d1v_rhs!`](@ref), without the
``n_x \times n_v`` array. On a workspace with the [`BoundedVelocityClosure`](@ref), ``1`` is
``\chi``, the indicator of the unknowns ``j = 2, …, n_v - 1``, which gives the density
``\rho_i = h_v \sum_{j=2}^{n_v-1} f_{ij}`` of [`_grid_1d1v_bounded_rhs!`](@ref). The potential of
``f = X L^\mathsf{T}`` is `_potential!(ϕ, X, L, w)`. The factors are `Matrix`es of the element
type `T` of `w`.

Throws a `DimensionMismatch` before any write when the sizes do not fit the grid of `w` and the
rank, and an `ArgumentError` when `ϕ` shares memory with a factor. At the rank 0 it writes zeros.
A non-finite entry propagates to the result, unchecked.
"""
function _potential!(ϕ::Vector{T}, K::Matrix{T}, V::Matrix{T},
        w::LowRankFields1D1V{T}) where {T}
    nx, nv = w.arakawa.nx, w.arakawa.nv
    length(ϕ) == size(K, 1) == nx && size(V, 1) == nv && size(K, 2) == size(V, 2) ||
        throw(DimensionMismatch(
            "_potential! on a $nx × $nv grid takes ϕ of length $nx, K of $nx rows and V of " *
            "$nv rows, with as many columns each; got ϕ of length $(length(ϕ)), K of size " *
            "$(size(K)) and V of size $(size(V))"))
    _check_alias("_potential!", "ϕ", ϕ, (K = K, V = V))
    return _green!(ϕ, K, _velocity_sum(V, w), w)
end

function _potential!(ϕ::Vector{T}, X::Matrix{T}, S::Matrix{T}, V::Matrix{T},
        w::LowRankFields1D1V{T}) where {T}
    nx, nv = w.arakawa.nx, w.arakawa.nv
    r = size(X, 2)
    length(ϕ) == nx && size(X) == (nx, r) && size(V) == (nv, r) && size(S) == (r, r) ||
        throw(DimensionMismatch(
            "_potential! on a $nx × $nv grid takes ϕ of length $nx, X of $nx rows, V of " *
            "$nv rows and S of size r × r, with r the columns of X and V; got ϕ of length " *
            "$(length(ϕ)), X of size $(size(X)), S of size $(size(S)) and V of size " *
            "$(size(V))"))
    _check_alias("_potential!", "ϕ", ϕ, (X = X, S = S, V = V))
    return _green!(ϕ, X, S * _velocity_sum(V, w), w)
end

# ϕ = G ρ with the density ρ = hv Y u
function _green!(ϕ, Y, u, w)
    ρ = Y * u
    ρ .*= w.arakawa.hv
    return mul!(ϕ, w.G, ρ)
end
