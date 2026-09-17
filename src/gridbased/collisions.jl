### Collision Tensor
# A N × N × N × tensor that discretizes the weak form of a Collision tensor at most quadratic in f: C[f] = C[f;f]
struct CollisionTensor{DT, FT} <: AbstractArray{DT, 3}
    nx::Int
    nv::Int
    f::FT

    function CollisionTensor(DT, nx, nv, f)
        new{DT, typeof(f)}(nx, nv, f)
    end
end

Base.size(ct::CollisionTensor) = tuple(ct.nx * ct.nv * ones(Int, 3)...)
Base.size(ct::CollisionTensor, i) = i ≥ 1 && i ≤ 3 ? ct.nx * ct.nv : 1

_nx(t::CollisionTensor) = t.nx
_nv(t::CollisionTensor) = t.nv

function Base.getindex(ct::CollisionTensor, I::CartesianIndex,
        J::CartesianIndex, K::CartesianIndex, L::CartesianIndex)
    @assert isvalid(I, ct.nx, ct.nv)
    @assert isvalid(J, ct.nx, ct.nv)
    @assert isvalid(K, ct.nx, ct.nv)

    ct.f(I, J, K)
end

function Base.getindex(ct::CollisionTensor, i::Int, j::Int, k::Int)
    I = multiindex(i, ct.nx, ct.nv)
    J = multiindex(j, ct.nx, ct.nv)
    K = multiindex(k, ct.nx, ct.nv)

    ct[I, J, K, L]
end

struct QuadraticCollisions{DT}
    nx::Int
    nv::Int
    hx::DT
    hv::DT
    v::AbstractVector
    factor::DT

    function QuadraticCollisions(
            nx::Int, nv::Int, hx::DT, hv::DT, v::AbstractVector) where {DT}
        new{DT}(nx, nv, hx, hv, v)
    end
end

function (quadraticcollisions::QuadraticCollisions{DT})(I, J, K) where {DT}
    i₁, j₁ = Tuple(I)
    i₂, j₂ = Tuple(J)
    i₃, j₃ = Tuple(K)

    # local in position space
    if i₁ != i₂ || i₁ != i₃
        return zero(DT)
    end

    v = quadraticcollisions.v

    j₋ = mod1(j₁-1, quadraticcollisions.nv)
    j₊ = mod1(j₁+1, quadraticcollisions.nv)

    if j₃ == j₁
        return -2 * v[j₂]^2
    elseif j₃ == j₊
        return v[j₂]^2 + v[j₊] * quadraticcollisions.hv / 2
    elseif j₃ == j₋
        return v[j₂]^2 - v[j₋] * quadraticcollisions.hv / 2
    else
        return zero(DT)
    end
end

struct ReducedCollisionTensor{DT, CT <: CollisionTensor{DT}, PM1, PM2, PM3} <:
       AbstractArray{DT, 3}
    tensor::CT
    projection_i::PM1
    projection_j::PM2
    projection_k::PM3

    function ReducedCollisionTensor(tensor::CollisionTensor{DT}, Pi::PM1,
            Pj::PM2, Pk::PM3) where {DT, PM1, PM2, PM3}
        @assert size(Pi, 1) == size(tensor, 1)
        @assert size(Pj, 1) == size(tensor, 2)
        @assert size(Pk, 1) == size(tensor, 3)
        new{DT, typeof(tensor), PM1, PM2, PM3}(tensor, Pi, Pj, Pk)
    end
end

function ReducedCollisionTensor(t::CollisionTensor{DT}, Pi::PM1) where {DT, PM1}
    ReducedCollisionTensor(t, Pi, Pi, Pi)
end

function Base.size(rt::ReducedCollisionTensor)
    (size(rt.projection_i, 2), size(rt.projection_j, 2), size(rt.projection_k, 2))
end
Base.size(rt::ReducedCollisionTensor, i) = size(rt)[i]
Base.axes(rt::ReducedCollisionTensor, i) = Base.OneTo(size(rt, i))

_nx(t::ReducedCollisionTensor) = _nx(t.tensor)
_nv(t::ReducedCollisionTensor) = _nv(t.tensor)

function _stencil_indices_v(i::Int, w::Int, nx, nv)
    indices = zeros(Int, 2*w + 1)
    ij = Tuple(multiindex(i, nx, nv))
    k = 1
    for d2 in (-w):w
        indices[k] = linearindex(CartesianIndex(mod1.(ij .+ (0, d2), (nx, nv))), nx, nv)
        k += 1
    end
    return indices
end

function Base.getindex(rt::ReducedCollisionTensor{DT}, i::Int, j::Int, k::Int) where {DT}
    @assert i ≥ 1 && i ≤ size(rt, 1)
    @assert j ≥ 1 && j ≤ size(rt, 2)
    @assert k ≥ 1 && k ≤ size(rt, 3)

    local nx = _nx(rt)
    local nv = _nv(rt)

    local r = zero(DT)

    for m1 in 1:nx
        for m2 in 1:nv
            m2₋ = mod1(m2-1, nv)
            m2₊ = mod1(m2+1, nv)
            for o2 in (m2₋, m2, m2₊)
                for n2 in 1:nv
                    M = CartesianIndex(m1, m2)
                    m = linearindex(M, nx, nv)
                    N = CartesianIndex(m1, n2)
                    n = linearindex(N, nx, nv)
                    O = CartesianIndex(m1, o2)
                    o = linearindex(O, nx, nv)

                    @inbounds r += rt.tensor[M, N, O] * rt.projection_i[m, i] *
                                   rt.projection_j[n, j] * rt.projection_k[o, k]
                end
            end
        end
    end

    return r
end

function _get_MC̃_cubic(V, ∫dv, ∫vdv, ∫v²dv, ci, li, h₁, h₂)
    local n₁, n₂ = size(ci)
    local n = n₁ * n₂
    local m = size(V, 2)

    @assert size(V, 1) == n

    MC̃ = zeros(m, m, m, m)

    ρ̂ = Matrix(∫dv * V)
    ρ̂u = Matrix(∫vdv * V)
    ρ̂ε = Matrix(∫v²dv * V)

    for ij in 1:n
        i, j = Tuple(multiindex(ij, n₁, n₂))
        j₊ = mod1(j+1, n₂)
        j₋ = mod1(j-1, n₂)
        for a in 1:m, b in 1:m, c in 1:m, o in 1:m
            @inbounds MC̃[o, a, b, c] += (1/h₂^2 *
                                         (ρ̂ε[i, a] * ρ̂[i, b] - ρ̂u[i, a] * ρ̂u[i, b]) *
                                         (V[li[i, j₋], c] - 2V[ij, c] + V[li[i, j₊], c])
                                         +
                                         1/(2h₂) *
                                         (v[j₊] * ρ̂[i, a] * ρ̂[i, b] - ρ̂[i, a] * ρ̂u[i, b]) *
                                         V[li[i, j₊], c]
                                         -
                                         1/(2h₂) *
                                         (v[j₋] * ρ̂[i, a] * ρ̂[i, b] - ρ̂[i, a] * ρ̂u[i, b]) *
                                         V[li[i, j₋], c]) * V[ij, o]
        end
    end

    return MC̃
end

function _get_MC̃_quadratic(V, ∫dv, ∫vdv, ∫v²dv, ci, li, h₁, h₂)
    local n₁, n₂ = size(ci)
    local n = n₁ * n₂
    local m = size(V, 2)

    @assert size(V, 1) == n

    MC̃ = zeros(m, m, m)
    ρ̂ε = Matrix(∫v²dv * V)
    ρ̂ = Matrix(∫dv * V)

    for ij in 1:n
        i, j = Tuple(multiindex(ij, n₁, n₂))
        j₊ = mod1(j+1, n₂)
        j₋ = mod1(j-1, n₂)
        for a in 1:m, c in 1:m, o in 1:m
            @inbounds MC̃[o, a, c] += (1/h₂^2 * ρ̂ε[i, a] *
                                      (V[li[i, j₋], c] - 2V[ij, c] + V[li[i, j₊], c])
                                      +
                                      1/(2h₂) * v[j₊] * ρ̂[i, a] * V[li[i, j₊], c]
                                      -
                                      1/(2h₂) * v[j₋] * ρ̂[i, a] * V[li[i, j₋], c]) *
                                     V[ij, o]
        end
    end

    return MC̃
end
