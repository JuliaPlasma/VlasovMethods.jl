### Potential Reduced Tensor
# A m × m × m₁ tensor:
# P̃[i,j,k] = ∑(l,m,n) P[l,m,n] Π¹[l,i] Π²[m,j] Π³[n,k]  

struct PotentialReducedTensor{DT, PT <: PoissonTensor{DT}, PM1, PM2, PM3} <:
       AbstractArray{DT, 3}
    tensor::PT
    projection_i::PM1
    projection_j::PM2
    projection_k::PM3

    function PotentialReducedTensor(
            tensor::PoissonTensor{DT}, Pi::PM1, Pj::PM2, Pk::PM3) where {DT, PM1, PM2, PM3}
        @assert size(Pi, 1) == size(tensor, 1)
        @assert size(Pj, 1) == size(tensor, 2)
        @assert size(Pk, 1) == size(tensor, 3)
        new{DT, typeof(tensor), PM1, PM2, PM3}(tensor, Pi, Pj, Pk)
    end
end

function Base.size(rt::PotentialReducedTensor)
    (size(rt.projection_i, 2), size(rt.projection_j, 2), size(rt.projection_k, 2))
end
Base.size(rt::PotentialReducedTensor, i) = size(rt)[i]
Base.axes(rt::PotentialReducedTensor, i) = Base.OneTo(size(rt, i))

_nx(t::PotentialReducedTensor) = _nx(t.tensor)
_nv(t::PotentialReducedTensor) = _nv(t.tensor)

function Base.getindex(rt::PotentialReducedTensor{DT}, i::Int, j::Int, k::Int) where {DT}
    @assert i ≥ 1 && i ≤ size(rt, 1)
    @assert j ≥ 1 && j ≤ size(rt, 2)
    @assert k ≥ 1 && k ≤ size(rt, 3)

    local nx = _nx(rt)
    local nv = _nv(rt)

    local r = zero(DT)

    # k1 here is the first index of the CartesianIndex K that describes the x,v space

    for l in 1:(nx * nv)
        nl = _stencil_indices(l, 1, nx, nv)
        for m in nl
            for n in nl
                r += rt.tensor[m, n, l] * rt.projection_i[m, i] * rt.projection_j[n, j] *
                     rt.projection_k[l, k]
            end
        end
    end

    return r
end

#function Base.materialize(rt::Union{ReducedTensor,PotentialReducedTensor})
#    [ rt[i,j,k] for i in axes(rt,1), j in axes(rt,2), k in axes(rt,3) ]
#end

### Velocity Reduced Matrix
# A m × m matrix where the first two indices are reduced using a projection on the phase space modes of f and the third dimension was contarcted with v²/2
# P̃[i,j] = ∑(l,m,(n₁,n₂)) P[l,m,n] Π¹[l,i] Π²[m,j] 1/2 v²[n₂(n)]

struct VelocityReducedMatrix{
    DT, PT <: PoissonTensor{DT}, PM1, PM2, PV <: AbstractVector{DT}} <: AbstractArray{DT, 2}
    tensor::PT
    projection_i::PM1
    projection_j::PM2
    v²::PV

    function VelocityReducedMatrix(tensor::PoissonTensor{DT}, Pi::PM1, Pj::PM2, v²::PV) where {
            DT, PM1, PM2, PV}
        @assert size(Pi, 1) == size(tensor, 1)
        @assert size(Pj, 1) == size(tensor, 2)
        new{DT, typeof(tensor), PM1, PM2, PV}(tensor, Pi, Pj, v²)
    end
end

Base.size(rt::VelocityReducedMatrix) = (size(rt.projection_i, 2), size(rt.projection_j, 2))
Base.size(rt::VelocityReducedMatrix, i) = size(rt)[i]
Base.axes(rt::VelocityReducedMatrix, i) = Base.OneTo(size(rt, i))

function Base.getindex(rt::VelocityReducedMatrix{DT}, i::Int, j::Int) where {DT}
    @assert i ≥ 1 && i ≤ size(rt, 1)
    @assert j ≥ 1 && j ≤ size(rt, 2)

    local nx = _nx(rt)
    local nv = _nv(rt)
    local r = zero(DT)

    for k in 1:(nx * nv)
        nk = _stencil_indices(k, 1, nx, nv)
        for m in nk
            for n in nk
                r += rt.tensor[m, n, k] * rt.projection_i[m, i] * rt.projection_j[n, j] *
                     rt.v²[k]
            end
        end
    end

    return r
end

_nx(t::VelocityReducedMatrix) = _nx(t.tensor)
_nv(t::VelocityReducedMatrix) = _nv(t.tensor)

### Fully Reduced Tensor
# A m × m × m₁ tensor where the first two indices are reduced using a projection on the phase space modes of f and the third on the modes of ϕ
# P̃[i,j,k] = ∑(l,m,(n₁,n₂)) P[l,m,n] Π¹[l,i] Π²[m,j] Π³[n₁(n),k]  

struct FullyReducedTensor{DT, PT <: PoissonTensor{DT}, PM1, PM2, PM3} <:
       AbstractArray{DT, 3}
    tensor::PT
    projection_i::PM1
    projection_j::PM2
    projection_k::PM3
    v::AbstractVector

    function FullyReducedTensor(tensor::PoissonTensor{DT}, Pi::PM1, Pj::PM2,
            Pα::PM3, v) where {DT, PM1, PM2, PM3}
        @assert size(Pi, 1) == size(tensor, 1)
        @assert size(Pj, 1) == size(tensor, 2)
        @assert size(Pk, 1) == size(tensor, 3)
        new{DT, typeof(tensor), PM1, PM2, PM3}(tensor, Pi, Pj, Pα, v)
    end
end

function Base.size(rt::FullyReducedTensor)
    (size(rt.projection_i, 2), size(rt.projection_j, 2), size(rt.projection_k, 2))
end
Base.size(rt::FullyReducedTensor, i) = size(rt)[i]
Base.axes(rt::FullyReducedTensor, i) = Base.OneTo(size(rt, i))

_nx(t::FullyReducedTensor) = _nx(t.tensor)
_nv(t::FullyReducedTensor) = _nv(t.tensor)

function Base.getindex(rt::FullyReducedTensor{DT}, i::Int, j::Int, k::Int) where {DT}
    @assert i ≥ 1 && i ≤ size(rt, 1)
    @assert j ≥ 1 && j ≤ size(rt, 2)
    @assert k ≥ 1 && k ≤ size(rt, 3)

    local nx = _nx(rt)
    local nv = _nv(rt)

    local r = zero(DT)

    # k1 here is the first index of the CartesianIndex K that describes the x,v space

    for k in 1:(nx * nv)
        nk = _stencil_indices(k, 1, nx, nv)
        k1 = Tuple(multiindex(k, nx, nv))[1]
        for m in nk
            for n in nk
                r += rt.tensor[m, n, k] * rt.projection_i[m, i] * rt.projection_j[n, j] *
                     rt.projection_k[k, α]
            end
        end
    end

    return r
end
