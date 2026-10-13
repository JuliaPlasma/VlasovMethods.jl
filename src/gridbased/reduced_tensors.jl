### Bracket Tensor

"""
    bracket_tensor(a::Arakawa{DT})

The rank-3 `GridTensor` of the bracket `a` on `a`'s own `nx × nv` grid, with stencil width 1:
its entry at the nodes `(I, J, K)` is `a(I, J, K)`. The grid is read from `a`, so the tensor
cannot be on another grid than its bracket. This is the tensor that `PotentialReducedTensor`,
`VelocityReducedMatrix` and `FullyReducedTensor` take.

`bracket_tensor` is a bridge until VlasovMethods requires GeometricBrackets 0.3, whose
four-argument `PoissonTensor(DT, nx, nv, f)` builds the same tensor and checks the grid.

```julia
julia> using GeometricBrackets: Arakawa

julia> size(bracket_tensor(Arakawa(3, 4, 0.5, 0.5)))
(12, 12, 12)
```
"""
function bracket_tensor(a::Arakawa{DT}) where {DT}
    # GeometricBrackets 0.2 has no public accessor of the grid an `Arakawa` is built on
    g = MultiIndexAxis(a.nx, a.nv)
    # `MultiIndexLazyArray` refuses a callable struct, so the coefficients go through a closure
    return GridTensor(MultiIndexLazyArray(DT, (I, J, K) -> a(I, J, K), g, g, g), 1)
end

### Potential Reduced Tensor

"""
    PotentialReducedTensor(tensor::GridTensor{DT, 3}, Pi, Pj, Pk)

The lazy `size(Pi, 2) × size(Pj, 2) × size(Pk, 2)` array that projects the first two indices of
the bracket tensor `tensor` onto the columns of `Pi` and `Pj`, and its third index onto the
columns of `Pk`. With `R = ReducedTensor(tensor, Pi, Pj)`,

    PotentialReducedTensor[i, j, k] = ∑_K R[i, j, K] * Pk[K, k] .

`Pi`, `Pj` and `Pk` are matrices with one row per grid node. The nodes are counted in the order
of the grid's own linear indices, the first direction fastest, as in `vec(dist.values)` of a
`GridDistribution`. The grid may have any dimension.

The constructor builds `ReducedTensor(tensor, Pi, Pj)` first, so the checks of that constructor
run here as well: `Pi` and `Pj` need one row per node, or it throws a `DimensionMismatch`, their
element type has to be `DT`, or it throws an `ArgumentError`, and every direction of the grid
needs at least `2w + 1` nodes for the tensor's stencil width `w`, or it throws an
`ArgumentError`. `Pk` needs one row per node, or the constructor throws a `DimensionMismatch`,
and the element type `DT`, or it throws an `ArgumentError`. An index off the array throws a
`BoundsError`.

The tensor's contract is `ReducedTensor`'s and is assumed, not checked: `tensor` is a bracket's
tensor on the grid that bracket was built on, nonzero only where its stencil says so.
`bracket_tensor` builds such a tensor from an `Arakawa`. A tensor written by hand on another grid
than its bracket breaks the contract silently.

`getindex` allocates nothing, and neither does materialising the array into a preallocated
output.

```julia
julia> using GeometricBrackets: Arakawa

julia> gt = bracket_tensor(Arakawa(3, 3, 0.5, 0.5));

julia> size(PotentialReducedTensor(gt, ones(9, 2), ones(9, 1), ones(9, 4)))
(2, 1, 4)
```
"""
struct PotentialReducedTensor{DT, RT <: ReducedTensor{DT}, PK} <: AbstractArray{DT, 3}
    reduced::RT
    projection_k::PK

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

# The extents `(nx, nv)` of the grid of `tensor`. An `x` index and a velocity index are one per
# node of a direction, so they need a two-dimensional grid.
function _phase_space_grid(tensor::GridTensor)
    ax = axes(parent(tensor), 1)
    ndims(ax) == 2 || throw(ArgumentError(
        "an x index and a velocity index need a two-dimensional grid, and this grid has " *
        "$(ndims(ax)) dimensions"))
    return size(ax)
end

### Velocity Reduced Matrix

"""
    VelocityReducedMatrix(tensor::GridTensor{DT, 3}, Pi, Pj, v²)

The lazy `size(Pi, 2) × size(Pj, 2)` matrix that contracts the third index of the bracket tensor
`tensor` with `v²/2`. With `R = ReducedTensor(tensor, Pi, Pj)` and `(n₁, n₂)` the grid indices
of the node `K`,

    VelocityReducedMatrix[i, j] = ∑_K R[i, j, K] * v²[n₂] / 2 .

`v²` holds the squared velocity nodes, `dist.v .^ 2` of a `GridDistribution`: one entry per
velocity node, in the order of `dist.v`. The result is `view(prt, :, :, 1)` of the
`PotentialReducedTensor` `prt` whose one-column projection carries `v²[n₂] / 2` in the row of
each node.

The grid needs two dimensions, or the function throws an `ArgumentError`. `v²` needs `nv`
entries, or it throws a `DimensionMismatch`, and the element type `DT`, or it throws an
`ArgumentError`. `PotentialReducedTensor` gives the other checks, the node order, the tensor's
contract and the allocations.

```julia
julia> size(VelocityReducedMatrix(gt, ones(9, 2), ones(9, 1), [1.0, 4.0, 9.0]))   # gt as above
(2, 1)
```
"""
function VelocityReducedMatrix(tensor::GridTensor{DT, 3}, Pi::AbstractMatrix,
        Pj::AbstractMatrix, v²::AbstractVector) where {DT}
    nx, nv = _phase_space_grid(tensor)
    length(v²) == nv || throw(DimensionMismatch(
        "v² needs $nv entries, one per velocity node, got $(length(v²))"))
    eltype(v²) == DT || throw(ArgumentError(
        "v² needs the element type $DT of the tensor, got $(eltype(v²))"))
    # the node `(n₁, n₂)` carries `v²[n₂] / 2`, and `/ 2` keeps a `Float32` weight in `Float32`
    w = reshape(repeat(v² ./ 2, inner = nx), :, 1)
    return view(PotentialReducedTensor(tensor, Pi, Pj, w), :, :, 1)
end

### Fully Reduced Tensor

"""
    FullyReducedTensor(tensor::GridTensor{DT, 3}, Pi, Pj, Px)

The lazy `size(Pi, 2) × size(Pj, 2) × size(Px, 2)` array that projects the third index of the
bracket tensor `tensor` onto the columns of `Px`, read at the `x` index of the node. With
`R = ReducedTensor(tensor, Pi, Pj)` and `(n₁, n₂)` the grid indices of the node `K`,

    FullyReducedTensor[i, j, k] = ∑_K R[i, j, K] * Px[n₁, k] .

`Px` has one row per `x` node, in the order of `dist.x`, and not one per grid node. The result
is the `PotentialReducedTensor` whose projection has the row `Px[n₁, :]` for each node.

The grid needs two dimensions, or the function throws an `ArgumentError`. `Px` needs `nx` rows,
or it throws a `DimensionMismatch`, and the element type `DT`, or it throws an `ArgumentError`.
`PotentialReducedTensor` gives the other checks, the node order, the tensor's contract and the
allocations.

```julia
julia> size(FullyReducedTensor(gt, ones(9, 2), ones(9, 1), ones(3, 4)))   # gt as above
(2, 1, 4)
```
"""
function FullyReducedTensor(tensor::GridTensor{DT, 3}, Pi::AbstractMatrix,
        Pj::AbstractMatrix, Px::AbstractMatrix) where {DT}
    nx, nv = _phase_space_grid(tensor)
    size(Px, 1) == nx || throw(DimensionMismatch(
        "Px needs $nx rows, one per x node, got $(size(Px, 1))"))
    eltype(Px) == DT || throw(ArgumentError(
        "Px needs the element type $DT of the tensor, got $(eltype(Px))"))
    # the first grid index runs fastest, so the stacked copies give the node `(n₁, n₂)` the row
    # `Px[n₁, :]`
    return PotentialReducedTensor(tensor, Pi, Pj, repeat(Px, nv, 1))
end
