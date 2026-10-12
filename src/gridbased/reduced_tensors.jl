### Potential Reduced Tensor

"""
    PotentialReducedTensor(tensor::GridTensor{DT, 3}, Pi, Pj, Pk)

The lazy `size(Pi, 2) × size(Pj, 2) × size(Pk, 2)` array that projects the first two indices of
the bracket tensor `tensor` onto the columns of `Pi` and `Pj`, and its third index onto the
columns of `Pk`. With `R = ReducedTensor(tensor, Pi, Pj)`,

    PotentialReducedTensor[i, j, k] = ∑_K R[i, j, K] * Pk[K, k] ,

so `Pk` has one row per grid node: `tensor` is a rank-3 `GridTensor`, `Pi` and `Pj` are matrices
with `size(tensor, 1)` and `size(tensor, 2)` rows, and `Pk` is a matrix with `size(tensor, 3)`
rows.

The constructor builds `ReducedTensor(tensor, Pi, Pj)` first, so the checks of that constructor
run here as well: `Pi` and `Pj` need their rows, or it throws a `DimensionMismatch`, their
element type has to be `DT`, or it throws an `ArgumentError`, and every direction of the grid
needs at least `2w + 1` nodes for the tensor's stencil width `w`, or it throws an
`ArgumentError`. `Pk` needs `size(tensor, 3)` rows, or the constructor throws a
`DimensionMismatch`, and the element type `DT`, or it throws an `ArgumentError`. An index off
the array throws a `BoundsError`.

The grid may have any dimension, because only the first two indices are reduced here. The nodes
are counted in the order of the grid's own linear indices, the first direction fastest, as in
`vec(dist.values)` of a `GridDistribution`.

The tensor's contract is `ReducedTensor`'s and is assumed, not checked: `tensor` is a bracket's
tensor on the grid that bracket was built on, nonzero only where its stencil says so. An
`Arakawa` on another grid breaks the contract silently, because nothing here can tell the two
grids apart.

```julia
julia> using MultiIndexArrays: GridTensor, MultiIndexLazyArray

julia> gt = GridTensor(MultiIndexLazyArray(Float64, (I, J, K) -> 1.0,
           (3, 3), (3, 3), (3, 3)), 1);   # the width-1 stencil covers the 3 × 3 grid

julia> prt = PotentialReducedTensor(gt, ones(9, 2), ones(9, 1), ones(9, 1));

julia> size(prt)
(2, 1, 1)
```

`getindex` allocates nothing, and neither does materialising the array into a preallocated
output.
"""
struct PotentialReducedTensor{DT, RT <: ReducedTensor{DT}, PM3} <: AbstractArray{DT, 3}
    reduced::RT
    projection_k::PM3

    function PotentialReducedTensor(tensor::GridTensor{DT, 3}, Pi::AbstractMatrix,
            Pj::AbstractMatrix, Pk::AbstractMatrix) where {DT}
        reduced = ReducedTensor(tensor, Pi, Pj)
        size(Pk, 1) == size(tensor, 3) || throw(DimensionMismatch(
            "Pk needs $(size(tensor, 3)) rows, one per grid node, got $(size(Pk, 1))"))
        eltype(Pk) == DT || throw(ArgumentError(
            "Pk needs the element type $DT of the tensor, got $(eltype(Pk))"))
        new{DT, typeof(reduced), typeof(Pk)}(reduced, Pk)
    end
end

function Base.size(rt::PotentialReducedTensor)
    (size(rt.reduced, 1), size(rt.reduced, 2), size(rt.projection_k, 2))
end

function Base.getindex(rt::PotentialReducedTensor{DT}, i::Int, j::Int, k::Int) where {DT}
    checkbounds(rt, i, j, k)

    local r = zero(DT)

    for K in axes(rt.reduced, 3)
        r += rt.reduced[i, j, K] * rt.projection_k[K, k]
    end

    return r
end

### Velocity Reduced Matrix

"""
    VelocityReducedMatrix(tensor::GridTensor{DT, 3}, Pi, Pj, v²)

The lazy `size(Pi, 2) × size(Pj, 2)` matrix that projects the first two indices of the bracket
tensor `tensor` onto the columns of `Pi` and `Pj`, and contracts its third index with `v²/2`.
With `R = ReducedTensor(tensor, Pi, Pj)` and `(n₁, n₂)` the grid indices of the node `K`,

    VelocityReducedMatrix[i, j] = ∑_K R[i, j, K] * v²[n₂] / 2 .

`v²` holds the **squared** velocity nodes of the grid, `dist.v .^ 2` of a `GridDistribution`, one
per velocity node, so it has `nv` entries and not one per phase-space node; the `½` is applied
inside the sum. The grid indices of a node are read on the grid that `tensor` is built on, the
first direction fastest, so `n₂` counts the velocity nodes in the order of `dist.v`.

The constructor builds `ReducedTensor(tensor, Pi, Pj)` first, so the checks of that constructor
run here as well: `Pi` and `Pj` need their rows, or it throws a `DimensionMismatch`, their
element type has to be `DT`, or it throws an `ArgumentError`, and every direction of the grid
needs at least `2w + 1` nodes for the tensor's stencil width `w`, or it throws an
`ArgumentError`. The grid needs two dimensions, because a velocity index is one per velocity
node, or the constructor throws an `ArgumentError`. `v²` needs `nv` entries, or it throws a
`DimensionMismatch`, and the element type `DT`, or it throws an `ArgumentError`. An index off
the array throws a `BoundsError`.

The tensor's contract is `ReducedTensor`'s and is assumed, not checked: `tensor` is a bracket's
tensor on the grid that bracket was built on, nonzero only where its stencil says so. An
`Arakawa` on another grid breaks the contract silently, because nothing here can tell the two
grids apart.

```julia
julia> using MultiIndexArrays: GridTensor, MultiIndexLazyArray

julia> gt = GridTensor(MultiIndexLazyArray(Float64, (I, J, K) -> 1.0,
           (3, 3), (3, 3), (3, 3)), 1);   # the width-1 stencil covers the 3 × 3 grid

julia> vrm = VelocityReducedMatrix(gt, ones(9, 2), ones(9, 1), [1.0, 2.0, 3.0]);

julia> size(vrm)
(2, 1)
```

`getindex` allocates nothing, and neither does materialising the array into a preallocated
output.
"""
struct VelocityReducedMatrix{DT, RT <: ReducedTensor{DT}, PV, AX <: MultiIndexAxis} <:
       AbstractArray{DT, 2}
    reduced::RT
    v²::PV
    axis::AX

    function VelocityReducedMatrix(tensor::GridTensor{DT, 3}, Pi::AbstractMatrix,
            Pj::AbstractMatrix, v²::AbstractVector) where {DT}
        reduced = ReducedTensor(tensor, Pi, Pj)
        ax = axes(parent(tensor), 1)
        ndims(ax) == 2 || throw(ArgumentError(
            "a velocity index is one per velocity node, which needs a two-dimensional grid, " *
            "and this grid has $(ndims(ax)) dimensions"))
        length(v²) == size(ax, 2) || throw(DimensionMismatch(
            "v² needs $(size(ax, 2)) entries, one per velocity node, got $(length(v²))"))
        eltype(v²) == DT || throw(ArgumentError(
            "v² needs the element type $DT of the tensor, got $(eltype(v²))"))
        new{DT, typeof(reduced), typeof(v²), typeof(ax)}(reduced, v², ax)
    end
end

Base.size(rt::VelocityReducedMatrix) = (size(rt.reduced, 1), size(rt.reduced, 2))

function Base.getindex(rt::VelocityReducedMatrix{DT}, i::Int, j::Int) where {DT}
    checkbounds(rt, i, j)

    local r = zero(DT)

    for K in axes(rt.reduced, 3)
        r += rt.reduced[i, j, K] * rt.v²[rt.axis[K][2]] / 2
    end

    return r
end

### Fully Reduced Tensor

"""
    FullyReducedTensor(tensor::GridTensor{DT, 3}, Pi, Pj, Px)

The lazy `size(Pi, 2) × size(Pj, 2) × size(Px, 2)` array that projects the first two indices of
the bracket tensor `tensor` onto the columns of `Pi` and `Pj`, and its third index onto the
columns of `Px`. With `R = ReducedTensor(tensor, Pi, Pj)` and `(n₁, n₂)` the grid indices of the
node `K`,

    FullyReducedTensor[i, j, k] = ∑_K R[i, j, K] * Px[n₁, k] ,

so `Px` has one row per `x` node: it is a matrix with `nx` rows and not one row per phase-space
node, and it is read at the `x` index of the node rather than at its linear index. The grid
indices are read on the grid that `tensor` is built on, the first direction fastest, so `n₁`
counts the `x` nodes in the order of `dist.x`.

The constructor builds `ReducedTensor(tensor, Pi, Pj)` first, so the checks of that constructor
run here as well: `Pi` and `Pj` need their rows, or it throws a `DimensionMismatch`, their
element type has to be `DT`, or it throws an `ArgumentError`, and every direction of the grid
needs at least `2w + 1` nodes for the tensor's stencil width `w`, or it throws an
`ArgumentError`. The grid needs two dimensions, because an `x` index is one per `x` node, or the
constructor throws an `ArgumentError`. `Px` needs `nx` rows, or it throws a
`DimensionMismatch`, and the element type `DT`, or it throws an `ArgumentError`. An index off
the array throws a `BoundsError`.

The tensor's contract is `ReducedTensor`'s and is assumed, not checked: `tensor` is a bracket's
tensor on the grid that bracket was built on, nonzero only where its stencil says so. An
`Arakawa` on another grid breaks the contract silently, because nothing here can tell the two
grids apart.

```julia
julia> using MultiIndexArrays: GridTensor, MultiIndexLazyArray

julia> gt = GridTensor(MultiIndexLazyArray(Float64, (I, J, K) -> 1.0,
           (3, 3), (3, 3), (3, 3)), 1);   # the width-1 stencil covers the 3 × 3 grid

julia> frt = FullyReducedTensor(gt, ones(9, 2), ones(9, 1), ones(3, 1));

julia> size(frt)
(2, 1, 1)
```

`getindex` allocates nothing, and neither does materialising the array into a preallocated
output.
"""
struct FullyReducedTensor{DT, RT <: ReducedTensor{DT}, PM3, AX <: MultiIndexAxis} <:
       AbstractArray{DT, 3}
    reduced::RT
    projection_x::PM3
    axis::AX

    function FullyReducedTensor(tensor::GridTensor{DT, 3}, Pi::AbstractMatrix,
            Pj::AbstractMatrix, Px::AbstractMatrix) where {DT}
        reduced = ReducedTensor(tensor, Pi, Pj)
        ax = axes(parent(tensor), 1)
        ndims(ax) == 2 || throw(ArgumentError(
            "an x index is one per x node, which needs a two-dimensional grid, and this grid " *
            "has $(ndims(ax)) dimensions"))
        size(Px, 1) == size(ax, 1) || throw(DimensionMismatch(
            "Px needs $(size(ax, 1)) rows, one per x node, got $(size(Px, 1))"))
        eltype(Px) == DT || throw(ArgumentError(
            "Px needs the element type $DT of the tensor, got $(eltype(Px))"))
        new{DT, typeof(reduced), typeof(Px), typeof(ax)}(reduced, Px, ax)
    end
end

function Base.size(rt::FullyReducedTensor)
    (size(rt.reduced, 1), size(rt.reduced, 2), size(rt.projection_x, 2))
end

function Base.getindex(rt::FullyReducedTensor{DT}, i::Int, j::Int, k::Int) where {DT}
    checkbounds(rt, i, j, k)

    local r = zero(DT)

    for K in axes(rt.reduced, 3)
        r += rt.reduced[i, j, K] * rt.projection_x[rt.axis[K][1], k]
    end

    return r
end
