@doc raw"""
    LowRankDistribution{DT, ST, XT, SMT, VT}

A distribution function of a 1D1V phase space in the factored form

```math
f \approx X S V^\mathsf{T} ,
```

a `DistributionFunction{DT, 1, 1}`. Its fields:

  - `space`, the factor space, a concrete subtype of
    [`VlasovMethods.AbstractFactorSpace`](@ref)`{DT}`, with `nx, nv = size(space)`;
  - `X`, the `nx × r` matrix of the x factors;
  - `S`, the `r × r` coefficient matrix, so the rank is `r = size(S, 1)`;
  - `V`, the `nv × r` matrix of the v factors;
  - `nfixed`, the number of leading columns of `V` that span the fixed velocity modes.

The type is mutable: a step may assign to `X`, `S`, `V` and `nfixed` factors of another rank of
the same array types. The assignment is not checked against the sizes. `size(d)` is
`size(d.space)` and `eltype(d)` is `DT`.

    LowRankDistribution(space, X, S, V; nfixed = 0)

Wrap the factors `X`, `S` and `V`, of the element type `DT` of `space`, without a copy. Throws a
`DimensionMismatch` for factors of other sizes, and an `ArgumentError` for a rank `r` outside
`1:min(nx, nv)` or `nfixed` outside `0:r`. It does not check that `X` and `V` are orthonormal:
that is the caller's contract.

    LowRankDistribution(dist::GridDistribution, r; fixed)

Factorise the node values `f = dist.values` at the rank `r` on
[`GridFactorSpace`](@ref)`(dist)`. `fixed` is an `nv × m` matrix of velocity modes, converted to
the element type of `dist`, none by default. With ``h_x``, ``h_v`` the steps of the grid:

 1. the fixed modes are orthonormalised first, by the weighted QR of `fixed`, so the first `k`
    columns of `V` span the first `k` columns of `fixed`, and `nfixed = m`;
 2. the next `q` columns of `V` span the leading right singular vectors of
    ``\sqrt{h_x h_v}\, g``, with ``g = f - h_v f V_m V_m^\mathsf{T}`` the part of `f` outside the
    fixed modes ``V_m``, where `q` is the smaller of `r - m` and the number of singular values
    above ``\max(n_x, n_v)\, \epsilon \, \| \sqrt{h_x h_v} f \|_F``, with ``\epsilon`` = `eps(DT)`;
 3. the columns `m + q + 1` to `r` complete `V` deterministically, from the full orthogonal factor
    of the Householder QR of ``\sqrt{h_v} [V_m \; V_q]``, with ``V_q`` the `q` singular vectors of
    item 2, divided by ``\sqrt{h_v}``; no random numbers are drawn;
 4. `X, S = orthonormalize_x(space, h_v f V)`.

So ``X S V^\mathsf{T} = h_v f V V^\mathsf{T}`` is the weighted orthogonal projection of `f` onto
the span of `V`, with ``h_x X^\mathsf{T} X = I`` and ``h_v V^\mathsf{T} V = I``. Without fixed
modes it is the weighted truncated SVD of rank `r`; with `m` fixed modes it keeps every moment
of `f` that the fixed modes span, such as the density, the momentum and the energy for the modes
`[1 v v²]`.

Throws an `ArgumentError` for `r` outside `1:min(nx, nv)`, for `m > r`, for linearly dependent
fixed modes, and for a non-finite entry of `f` or of `fixed`; and a `DimensionMismatch` for
`fixed` of another number of rows than `nv`.
"""
mutable struct LowRankDistribution{
    DT, ST <: AbstractFactorSpace{DT}, XT <: AbstractMatrix{DT},
    SMT <: AbstractMatrix{DT}, VT <: AbstractMatrix{DT}} <: DistributionFunction{DT, 1, 1}
    const space::ST
    X::XT
    S::SMT
    V::VT
    nfixed::Int

    function LowRankDistribution(space::ST, X::XT, S::SMT, V::VT;
            nfixed::Integer = 0) where {DT, ST <: AbstractFactorSpace{DT},
            XT <: AbstractMatrix{DT}, SMT <: AbstractMatrix{DT}, VT <: AbstractMatrix{DT}}
        nx, nv = size(space)
        r = size(S, 1)
        _check_rank(r, nx, nv)
        size(X) == (nx, r) && size(S) == (r, r) && size(V) == (nv, r) ||
            throw(DimensionMismatch("factors of sizes $(size(X)), $(size(S)) and " *
                                    "$(size(V)) in a space of $nx × $nv nodes"))
        0 ≤ nfixed ≤ r || throw(ArgumentError(
            "the number of fixed modes $nfixed is not in 0:$r for the rank $r"))
        new{DT, ST, XT, SMT, VT}(space, X, S, V, nfixed)
    end
end

function _check_rank(r, nx, nv)
    1 ≤ r ≤ min(nx, nv) || throw(ArgumentError(
        "the rank $r is not in 1:$(min(nx, nv)) for a space of $nx × $nv nodes"))
    return nothing
end

function LowRankDistribution(dist::GridDistribution{DT}, r::Integer;
        fixed::AbstractMatrix = zeros(DT, length(dist.v), 0)) where {DT}
    space = GridFactorSpace(dist)
    nx, nv = size(space)
    f = dist.values
    _check_rank(r, nx, nv)
    size(fixed, 1) == nv || throw(DimensionMismatch(
        "fixed modes of $(size(fixed, 1)) rows on a velocity grid of $nv nodes"))
    m = size(fixed, 2)
    m ≤ r || throw(ArgumentError("$m fixed modes exceed the rank $r"))
    all(isfinite, f) ||
        throw(ArgumentError("the values of the distribution hold a non-finite entry"))
    all(isfinite, fixed) || throw(ArgumentError("the fixed modes hold a non-finite entry"))

    hx, hv = step(space.x), step(space.v)
    Φ = Matrix{DT}(fixed)
    if m > 0
        σ = svdvals(sqrt(hv) .* Φ)
        last(σ) > nv * eps(DT) * first(σ) ||
            throw(ArgumentError("the fixed modes are linearly dependent"))
    end
    Vm, _ = orthonormalize_v(space, Φ)

    # the leading right singular vectors of the part of f outside the fixed modes
    w = sqrt(hx * hv)
    g = f - (hv .* (f * Vm)) * Vm'
    F = svd!(w .* g)
    q = min(r - m, count(>(max(nx, nv) * eps(DT) * w * norm(f)), F.S))

    # the first r columns of the full orthogonal factor: orthonormal to working precision, and
    # a deterministic completion below the rank
    s = sqrt(hv)
    B = qr(s .* hcat(Vm, F.V[:, 1:q]))
    V = lmul!(B.Q, Matrix{DT}(I, nv, r)) ./ s

    X, S = orthonormalize_x(space, hv .* (f * V))
    return LowRankDistribution(space, X, S, V; nfixed = m)
end

Base.eltype(::LowRankDistribution{DT}) where {DT} = DT
Base.size(d::LowRankDistribution) = size(d.space)

"""
    GridDistribution(d::LowRankDistribution)

The [`GridDistribution`](@ref) of the values `X S Vᵀ` on the grid of the
[`GridFactorSpace`](@ref) of `d`. It forms the `nx × nv` array, so it serves tests and
diagnostics.
"""
function GridDistribution(d::LowRankDistribution{DT, <:GridFactorSpace}) where {DT}
    GridDistribution(d.space.x, d.space.v, d.X * d.S * d.V')
end

@doc raw"""
    velocity_moments(d::LowRankDistribution)

The velocity moments `(density, momentum, energy)` of `d` on a [`GridFactorSpace`](@ref), in the
convention of `velocity_moments` of a [`GridDistribution`](@ref), computed from the factors
without the `nx × nv` array:

```math
n = h_v X S (V^\mathsf{T} 1) , \qquad
p = h_v X S (V^\mathsf{T} v) , \qquad
\varepsilon = \tfrac{1}{2} h_v X S (V^\mathsf{T} v^2) .
```
"""
function velocity_moments(d::LowRankDistribution{DT, <:GridFactorSpace}) where {DT}
    v = collect(d.space.v)
    hv = step(d.space.v)
    XS = d.X * d.S
    n = XS * (hv .* vec(sum(d.V; dims = 1)))
    j = XS * (hv .* (d.V' * v))
    ε = XS * ((hv / 2) .* (d.V' * (v .^ 2)))
    return (density = n, momentum = j, energy = ε)
end

"""
    evaluate(d::LowRankDistribution, x, v)

The bilinear interpolant of the grid values `X S Vᵀ` of `d` on a [`GridFactorSpace`](@ref) at
`(only(x), only(v))`, with the index rule of the `GridDistribution`: periodic in `x`, zero
outside the velocity range, and `v` at its upper end in the last cell. It interpolates two rows
of `X` and two rows of `V` and contracts them with `S`, without the `nx × nv` array.
"""
function evaluate(d::LowRankDistribution{DT, <:GridFactorSpace}, x, v) where {DT}
    xs, vs = d.space.x, d.space.v
    nx, nv = size(d.space)
    first(vs) ≤ only(v) ≤ last(vs) || return zero(DT)
    η = (only(v) - first(vs)) / step(vs)
    ξ = (only(x) - first(xs)) / step(xs)

    ξ₀ = floor(ξ)
    t = ξ - ξ₀
    i = mod(Int(ξ₀), nx) + 1
    i⁺ = mod(i, nx) + 1

    # `v = v_b` falls in the last cell, at its upper end
    j = min(floor(Int, η), nv - 2) + 1
    s = η - (j - 1)

    X, S, V = d.X, d.S, d.V
    value = zero(promote_type(DT, typeof(t), typeof(s)))
    for b in axes(S, 2)
        vb = (1 - s) * V[j, b] + s * V[j + 1, b]
        for a in axes(S, 1)
            value += ((1 - t) * X[i, a] + t * X[i⁺, a]) * S[a, b] * vb
        end
    end
    return value
end
