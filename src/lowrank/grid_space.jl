"""
    AbstractFactorSpace{DT}

The supertype of the spaces that hold the factors of a [`LowRankDistribution`](@ref), with the
element type `DT`. A concrete space defines `size(space) = (nx, nv)`, the lengths of the x and
the v factors, and extends the four functions of its weighted inner products:
[`inner_x`](@ref VlasovMethods.inner_x), [`inner_v`](@ref VlasovMethods.inner_v),
[`orthonormalize_x`](@ref VlasovMethods.orthonormalize_x) and
[`orthonormalize_v`](@ref VlasovMethods.orthonormalize_v). `eltype(space)` is `DT`.
"""
abstract type AbstractFactorSpace{DT} end

Base.eltype(::Type{<:AbstractFactorSpace{DT}}) where {DT} = DT

"""
    inner_x(space, A, B)

The weighted inner product of the x factors `A` and `B` in `space`: a scalar for two vectors,
the matrix of the inner products of the columns for two matrices.
"""
function inner_x end

"""
    inner_v(space, A, B)

The weighted inner product of the v factors `A` and `B` in `space`: a scalar for two vectors,
the matrix of the inner products of the columns for two matrices.
"""
function inner_v end

"""
    orthonormalize_x(space, A)

The factors `(Q, R)` of `A = Q R`, where the columns of `Q` are orthonormal in the weighted
inner product [`inner_x`](@ref VlasovMethods.inner_x) of `space` and `R` is upper triangular.
"""
function orthonormalize_x end

"""
    orthonormalize_v(space, A)

The factors `(Q, R)` of `A = Q R`, where the columns of `Q` are orthonormal in the weighted
inner product [`inner_v`](@ref VlasovMethods.inner_v) of `space` and `R` is upper triangular.
"""
function orthonormalize_v end

@doc raw"""
    GridFactorSpace(x, v)
    GridFactorSpace(dist::GridDistribution)

The factor space of the nodes of a [`GridDistribution`](@ref): the periodic x grid `x` and the
bounded v grid `v`, two ranges of one element type `DT`. `x` holds at least one node and `v` at
least two. `size(space)` is `(length(x), length(v))`.

The inner products are the rectangle rules of the grid, with ``h_x = `` `step(x)` and
``h_v = `` `step(v)`:

```math
\langle a, b \rangle_x = h_x \sum_i a_i b_i , \qquad
\langle a, b \rangle_v = h_v \sum_j a_j b_j ,
```

so the rule in `v` is the rule of `velocity_moments` of a `GridDistribution`.
[`VlasovMethods.inner_x`](@ref) and [`VlasovMethods.inner_v`](@ref) evaluate them, and
[`VlasovMethods.orthonormalize_x`](@ref) and [`VlasovMethods.orthonormalize_v`](@ref)
orthonormalise in them by the Householder QR of `LinearAlgebra.qr`.
"""
struct GridFactorSpace{DT, XT <: AbstractRange{DT}, VT <: AbstractRange{DT}} <:
       AbstractFactorSpace{DT}
    x::XT
    v::VT

    function GridFactorSpace(x::XT,
            v::VT) where {DT, XT <: AbstractRange{DT}, VT <: AbstractRange{DT}}
        isempty(x) &&
            throw(ArgumentError("a periodic x grid needs at least one node, got 0"))
        length(v) ≥ 2 || throw(ArgumentError(
            "a bounded velocity grid needs at least two nodes, got $(length(v))"))
        new{DT, XT, VT}(x, v)
    end
end

GridFactorSpace(dist::GridDistribution) = GridFactorSpace(dist.x, dist.v)

Base.size(space::GridFactorSpace) = (length(space.x), length(space.v))

function _check_rows(A, n, axis)
    size(A, 1) == n || throw(DimensionMismatch(
        "a factor of $(size(A, 1)) rows on the $axis grid of $n nodes"))
    return nothing
end

function _weighted_inner(A, B, h, n, axis)
    _check_rows(A, n, axis)
    _check_rows(B, n, axis)
    return h * (A' * B)
end

# Q = Q̃ / √h and R = R̃ of the thin Householder QR √h A = Q̃ R̃, so A = Q R and h QᵀQ = I.
# The working matrix is dense: the QR of a sparse matrix permutes its rows and columns, which
# Q and R do not carry. `convert` copies nothing for the dense `s .* A` of a dense `A`.
function _weighted_qr(A, h, n, axis)
    _check_rows(A, n, axis)
    k = size(A, 2)
    k ≤ n || throw(ArgumentError(
        "cannot orthonormalise $k columns on the $axis grid of $n nodes"))
    s = sqrt(h)
    F = qr(convert(Matrix, s .* A))
    return Matrix(F.Q) ./ s, F.R
end

"""
    inner_x(space::GridFactorSpace, A, B)

`h_x AᵀB` for the vectors or matrices `A` and `B` of `nx` rows, with `h_x = step(space.x)`.
Throws a `DimensionMismatch` for another number of rows.
"""
function inner_x(space::GridFactorSpace, A::AbstractVecOrMat, B::AbstractVecOrMat)
    _weighted_inner(A, B, step(space.x), length(space.x), "x")
end

"""
    inner_v(space::GridFactorSpace, A, B)

`h_v AᵀB` for the vectors or matrices `A` and `B` of `nv` rows, with `h_v = step(space.v)`.
Throws a `DimensionMismatch` for another number of rows.
"""
function inner_v(space::GridFactorSpace, A::AbstractVecOrMat, B::AbstractVecOrMat)
    _weighted_inner(A, B, step(space.v), length(space.v), "v")
end

"""
    orthonormalize_x(space::GridFactorSpace, A)

`(Q, R)` with `Q = Q̃ / √h_x` and `R = R̃` from the thin Householder QR `√h_x A = Q̃ R̃`, for `A`
of size `nx × k`: `A = Q R`, `h_x QᵀQ = I`, `Q` of size `nx × k` and `R` of size `k × k`. Throws
an `ArgumentError` for `k > nx`, and a `DimensionMismatch` for another number of rows than `nx`.
The entries of `A` are not checked: a non-finite entry propagates into `Q` and `R`.
"""
function orthonormalize_x(space::GridFactorSpace, A::AbstractMatrix)
    _weighted_qr(A, step(space.x), length(space.x), "x")
end

"""
    orthonormalize_v(space::GridFactorSpace, A)

`(Q, R)` with `Q = Q̃ / √h_v` and `R = R̃` from the thin Householder QR `√h_v A = Q̃ R̃`, for `A`
of size `nv × k`: `A = Q R`, `h_v QᵀQ = I`, `Q` of size `nv × k` and `R` of size `k × k`. Throws
an `ArgumentError` for `k > nv`, and a `DimensionMismatch` for another number of rows than `nv`.
The entries of `A` are not checked: a non-finite entry propagates into `Q` and `R`.
"""
function orthonormalize_v(space::GridFactorSpace, A::AbstractMatrix)
    _weighted_qr(A, step(space.v), length(space.v), "v")
end
