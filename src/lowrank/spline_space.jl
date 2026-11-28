@doc raw"""
    SplineFactorSpace(xbasis, vbasis)
    SplineFactorSpace(dist::TensorSplineDistribution)

The factor space of the coefficients of a tensor spline: the periodic B-spline basis `xbasis`
``\{\varphi_i\}`` in x and the B-spline basis `vbasis` ``\{\psi_j\}`` in v, the bases that a
[`TensorSplineDistribution`](@ref) takes. A factor is a matrix of coefficient vectors, an x
factor of `nx = nbasis(xbasis)` rows and a v factor of `nv = nbasis(vbasis)` rows, so
`size(space)` is `(nx, nv)`. The second form takes the bases of `dist`.

The space holds the bases, their default `SplineQuadrature`s `xquadrature` and `vquadrature`,
the mass matrices `Mx` and `Mv`, `mass_matrix` of the two quadratures, and their upper Cholesky
factors `Ux` and `Uv`, sparse and without a permutation, with ``M_x = U_x^\mathsf{T} U_x`` and
``M_v = U_v^\mathsf{T} U_v``. The inner products are those of ``L^2``,

```math
\langle a, b \rangle_x = a^\mathsf{T} M_x b , \qquad
\langle a, b \rangle_v = a^\mathsf{T} M_v b ,
```

which [`VlasovMethods.inner_x`](@ref) and [`VlasovMethods.inner_v`](@ref) evaluate, and
[`VlasovMethods.orthonormalize_x`](@ref) and [`VlasovMethods.orthonormalize_v`](@ref)
orthonormalise in them by the Householder QR of the Cholesky factor times the factor.

The constructor throws, as [`TensorSplineDistribution`](@ref) does:

  - `ArgumentError` unless `xbasis` is a `PeriodicBSplineBasis`;
  - `ArgumentError` unless `vbasis` reproduces ``v^2``, `polynomial_reproduction(vbasis) ≥ 2`: a
    Dirichlet, a periodic and a linear clamped basis do not;
  - `ArgumentError` unless the two bases have one element type.

The velocity domain may be any interval; only [`VlasovMethods.reflect_v`](@ref) needs it
symmetric.
"""
struct SplineFactorSpace{T, XB <: PeriodicBSplineBasis{T}, VB <: AbstractBSplineBasis{T},
    XQ <: SplineQuadrature{T}, VQ <: SplineQuadrature{T}, MT <: AbstractMatrix{T},
    UT <: AbstractMatrix{T}} <: AbstractFactorSpace{T}
    xbasis::XB
    vbasis::VB
    xquadrature::XQ
    vquadrature::VQ
    Mx::MT
    Mv::MT
    Ux::UT
    Uv::UT

    function SplineFactorSpace(xbasis::AbstractBSplineBasis, vbasis::AbstractBSplineBasis)
        xbasis isa PeriodicBSplineBasis || throw(ArgumentError(
            "the x-basis of a SplineFactorSpace is periodic, a PeriodicBSplineBasis; got a " *
            "$(nameof(typeof(xbasis)))"))
        polynomial_reproduction(vbasis) ≥ 2 || throw(ArgumentError(
            "the v-basis of a SplineFactorSpace reproduces v², so that the kinetic energy " *
            "v²/2 is in its span: a clamped (Free()) basis of degree 2 or more; got a basis " *
            "that reproduces polynomials of degree $(polynomial_reproduction(vbasis))"))
        eltype(xbasis) === eltype(vbasis) || throw(ArgumentError(
            "the bases of a SplineFactorSpace have one element type; got $(eltype(xbasis)) " *
            "and $(eltype(vbasis))"))
        xq, vq = SplineQuadrature(xbasis), SplineQuadrature(vbasis)
        Mx, Mv = mass_matrix(xq), mass_matrix(vq)
        Ux, Uv = _upper_cholesky(Mx), _upper_cholesky(Mv)
        new{eltype(xbasis), typeof(xbasis), typeof(vbasis), typeof(xq), typeof(vq),
            typeof(Mx), typeof(Ux)}(xbasis, vbasis, xq, vq, Mx, Mv, Ux, Uv)
    end
end

function SplineFactorSpace(dist::TensorSplineDistribution)
    SplineFactorSpace(dist.xbasis, dist.vbasis)
end

# The upper factor U of M = UᵀU. CHOLMOD permutes the rows and columns by default, and its factor
# then belongs to M[p, p]; the identity permutation keeps the factor of M.
function _upper_cholesky(M::AbstractMatrix)
    F = cholesky(Symmetric(M); perm = 1:size(M, 1))
    return copy(transpose(sparse(F.L)))
end

Base.size(space::SplineFactorSpace) = (nbasis(space.xbasis), nbasis(space.vbasis))

function _check_basis_rows(A, n, axis)
    size(A, 1) == n || throw(DimensionMismatch(
        "a factor of $(size(A, 1)) rows on the $axis basis of $n functions"))
    return nothing
end

# Q = U⁻¹ Q̃ and R from the thin Householder QR U A = Q̃ R, so A = Q R and QᵀM Q = I for M = UᵀU
function _cholesky_qr(A, U, n, axis)
    _check_basis_rows(A, n, axis)
    k = size(A, 2)
    k ≤ n || throw(ArgumentError(
        "cannot orthonormalise $k columns on the $axis basis of $n functions"))
    F = qr(convert(Matrix, U * A))
    return UpperTriangular(U) \ Matrix(F.Q), F.R
end

"""
    inner_x(space::SplineFactorSpace, A, B)

`Aᵀ Mx B` for the vectors or matrices `A` and `B` of `nx` rows, the ``L^2`` inner products of
the x splines of their columns. Throws a `DimensionMismatch` for another number of rows.
"""
function inner_x(space::SplineFactorSpace, A::AbstractVecOrMat, B::AbstractVecOrMat)
    nx = size(space)[1]
    _check_basis_rows(A, nx, "x")
    _check_basis_rows(B, nx, "x")
    return A' * (space.Mx * B)
end

"""
    inner_v(space::SplineFactorSpace, A, B)

`Aᵀ Mv B` for the vectors or matrices `A` and `B` of `nv` rows, the ``L^2`` inner products of
the v splines of their columns. Throws a `DimensionMismatch` for another number of rows.
"""
function inner_v(space::SplineFactorSpace, A::AbstractVecOrMat, B::AbstractVecOrMat)
    nv = size(space)[2]
    _check_basis_rows(A, nv, "v")
    _check_basis_rows(B, nv, "v")
    return A' * (space.Mv * B)
end

"""
    orthonormalize_x(space::SplineFactorSpace, A)

`(Q, R)` with `Q = Ux⁻¹ Q̃` and `R = R̃` from the thin Householder QR `Ux A = Q̃ R̃`, for `A` of
size `nx × k`: `A = Q R`, `Qᵀ Mx Q = I`, `Q` of size `nx × k` and `R` of size `k × k`. The QR of
a zero matrix gives an orthonormal `Q` and a zero `R`. Throws an `ArgumentError` for `k > nx`,
and a `DimensionMismatch` for another number of rows than `nx`. The entries of `A` are not
checked: a non-finite entry propagates into `Q` and `R`.
"""
function orthonormalize_x(space::SplineFactorSpace, A::AbstractMatrix)
    _cholesky_qr(A, space.Ux, size(space)[1], "x")
end

"""
    orthonormalize_v(space::SplineFactorSpace, A)

`(Q, R)` with `Q = Uv⁻¹ Q̃` and `R = R̃` from the thin Householder QR `Uv A = Q̃ R̃`, for `A` of
size `nv × k`: `A = Q R`, `Qᵀ Mv Q = I`, `Q` of size `nv × k` and `R` of size `k × k`. The QR of
a zero matrix gives an orthonormal `Q` and a zero `R`. Throws an `ArgumentError` for `k > nv`,
and a `DimensionMismatch` for another number of rows than `nv`. The entries of `A` are not
checked: a non-finite entry propagates into `Q` and `R`.
"""
function orthonormalize_v(space::SplineFactorSpace, A::AbstractMatrix)
    _cholesky_qr(A, space.Uv, size(space)[2], "v")
end

@doc raw"""
    reflect_v(space::SplineFactorSpace, A)

The reflection ``v \mapsto -v`` of the velocity coefficients `A`, a vector or a matrix of `nv`
rows: `A[nv:-1:1, :]`. On a v basis of a uniform mesh of a symmetric interval ``[-a, a]`` the
mirror of the function ``\psi_j`` is ``\psi_{n_v + 1 - j}``, ``\psi_j(-v) = \psi_{n_v+1-j}(v)``,
so the reversal of the coefficients is the reflection of the spline, and it commutes with `Mv`,
both up to the rounding of the breakpoints.
A coefficient vector `a` is even when `reflect_v(space, a) == a` and odd when it is `-a`.

Throws an `ArgumentError` unless the v basis lies on a `UniformMesh` whose domain has the ends
``-a`` and ``a``, compared exactly, and a `DimensionMismatch` for another number of rows than
`nv`.
"""
function reflect_v(space::SplineFactorSpace, A::AbstractVecOrMat)
    vbasis = space.vbasis
    a, b = leftendpoint(domain(vbasis)), rightendpoint(domain(vbasis))
    SimpleSplines.mesh(vbasis) isa UniformMesh && a == -b || throw(ArgumentError(
        "the reflection of the velocity coefficients needs a uniform velocity mesh on a " *
        "symmetric interval [-a, a]; got a $(nameof(typeof(SimpleSplines.mesh(vbasis)))) on " *
        "$a .. $b"))
    _check_basis_rows(A, size(space)[2], "v")
    return reverse(A; dims = 1)
end

@doc raw"""
    LowRankDistribution(space::SplineFactorSpace, F::AbstractMatrix, r; fixed)
    LowRankDistribution(dist::TensorSplineDistribution, r; fixed)

Factorise the `nx × nv` coefficient matrix `F` at the rank `r` on the [`SplineFactorSpace`](@ref)
`space`, in the ``L^2`` inner products of the space, with ``M_v = U_v^\mathsf{T} U_v`` and
``M_x = U_x^\mathsf{T} U_x``. `fixed` is an `nv × m` matrix of velocity coefficient vectors,
converted to the element type `T` of `space`, none by default:

 1. the fixed modes are orthonormalised first, by
    [`orthonormalize_v`](@ref VlasovMethods.orthonormalize_v), so the first `k` columns of `V`
    span the first `k` columns of `fixed`, and `nfixed = m`;
 2. the next `q` columns of `V` span the leading right singular vectors of
    ``U_x g \, U_v^\mathsf{T}``, mapped by ``U_v^{-1}``, with ``g = F - F M_v V_m V_m^\mathsf{T}``
    the part of `F` outside the fixed modes ``V_m``, where `q` is the smaller of `r - m` and the
    number of singular values above ``\max(n_x, n_v)\, \epsilon \, \| U_x F U_v^\mathsf{T} \|_F``,
    with ``\epsilon`` = `eps(T)`;
 3. the columns `m + q + 1` to `r` complete `V` deterministically, from the full orthogonal
    factor of the Householder QR of ``U_v [V_m \; V_q]``, mapped by ``U_v^{-1}``; no random
    numbers are drawn;
 4. `X, S = orthonormalize_x(space, F Mv V)`.

So ``X S V^\mathsf{T} = F M_v V V^\mathsf{T}`` is the ``L^2``-orthogonal projection of the spline
of `F` onto the span of `V`, with ``X^\mathsf{T} M_x X = I`` and ``V^\mathsf{T} M_v V = I``.
Without fixed modes it is the truncated SVD of rank `r` in the ``L^2`` norm
``\| U_x F U_v^\mathsf{T} \|_F``; with `m` fixed modes it keeps every velocity moment of `F` that
the fixed modes span, such as the density, the momentum and the kinetic energy for the
coefficients of `1`, `v` and `v²`. The factorisation forms dense `nx × nv` arrays, so it serves
the start of a run, tests and diagnostics.

The second form is the first on [`SplineFactorSpace`](@ref)`(dist)` and the coefficients of
`dist`.

Throws an `ArgumentError` for `r` outside `1:min(nx, nv)`, for `m > r`, for linearly dependent
fixed modes, and for a non-finite entry of `F` or of `fixed`; and a `DimensionMismatch` for `F`
of another size than `(nx, nv)` and for `fixed` of another number of rows than `nv`.
"""
function LowRankDistribution(space::SplineFactorSpace{T}, F::AbstractMatrix, r::Integer;
        fixed::AbstractMatrix = zeros(T, size(space)[2], 0)) where {T}
    nx, nv = size(space)
    _check_rank(r, nx, nv)
    size(F) == (nx, nv) || throw(DimensionMismatch(
        "a coefficient matrix of size $(size(F)) on a spline space of $nx × $nv functions"))
    size(fixed, 1) == nv || throw(DimensionMismatch(
        "fixed modes of $(size(fixed, 1)) rows on a velocity basis of $nv functions"))
    m = size(fixed, 2)
    m ≤ r || throw(ArgumentError("$m fixed modes exceed the rank $r"))
    all(isfinite, F) ||
        throw(ArgumentError("the coefficient matrix holds a non-finite entry"))
    all(isfinite, fixed) || throw(ArgumentError("the fixed modes hold a non-finite entry"))

    (; Mv, Ux, Uv) = space
    f = convert(Matrix{T}, F)
    Φ = convert(Matrix{T}, fixed)
    if m > 0
        σ = svdvals(convert(Matrix, Uv * Φ))
        last(σ) > nv * eps(T) * first(σ) ||
            throw(ArgumentError("the fixed modes are linearly dependent"))
    end
    Vm, _ = orthonormalize_v(space, Φ)

    # the leading right singular vectors of the part of f outside the fixed modes, in L²
    g = f - (f * (Mv * Vm)) * Vm'
    G = svd!(convert(Matrix, (Ux * g) * Uv'))
    q = min(r - m, count(>(max(nx, nv) * eps(T) * norm((Ux * f) * Uv')), G.S))

    # the first r columns of the full orthogonal factor: orthonormal to working precision, and
    # a deterministic completion below the rank
    B = qr(hcat(convert(Matrix, Uv * Vm), G.V[:, 1:q]))
    V = UpperTriangular(Uv) \ lmul!(B.Q, Matrix{T}(I, nv, r))

    X, S = orthonormalize_x(space, f * (Mv * V))
    return LowRankDistribution(space, X, S, V; nfixed = m)
end

function LowRankDistribution(dist::TensorSplineDistribution{T}, r::Integer;
        fixed::AbstractMatrix = zeros(T, nbasis(dist.vbasis), 0)) where {T}
    LowRankDistribution(SplineFactorSpace(dist), dist.coefficients, r; fixed)
end

"""
    TensorSplineDistribution(d::LowRankDistribution)

The [`TensorSplineDistribution`](@ref) of the coefficients `X S Vᵀ` on the bases of the
[`SplineFactorSpace`](@ref) of `d`. It forms the `nx × nv` array, so it serves tests and
diagnostics.
"""
function TensorSplineDistribution(d::LowRankDistribution{T, <:SplineFactorSpace}) where {T}
    TensorSplineDistribution(d.space.xbasis, d.space.vbasis, d.X * d.S * d.V')
end
