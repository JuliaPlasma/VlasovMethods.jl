# The three BUG integrators transcribed from their papers on the formed grid function
# `f = X S Vᵀ`, the oracle of `test/lowrank/bug.jl`. The field is the full-order right-hand side
# `F(f)` of the base's generic `_grid_1d1v_rhs!`, applied to the formed `f` and projected onto
# the other factor; the weighted orthonormal bases come from `LinearAlgebra.qr` of `√h A`; each
# substep is one step of the classical RK4 tableau.
#
#   - BUG: Ceruti & Lubich 2022, doi 10.1007/s10543-021-00873-0, arXiv 2010.02022 §3.1;
#   - augmented BUG: Ceruti, Kusch & Lubich 2022, doi 10.1007/s10543-021-00907-7,
#     arXiv 2104.05247 §2.1;
#   - midpoint BUG (4r): Ceruti, Einkemmer, Kusch & Lubich 2024, doi 10.1007/s10543-024-01032-x,
#     arXiv 2402.08607v1 §3, eqs. (10)–(13), truncated to the rank of its input.
#
# `d` is a named tuple `(p, hx, hv)` with `p = (arakawa, v, G)`, the parameters that the generic
# method of `_grid_1d1v_rhs!` reads.

using LinearAlgebra: Diagonal, norm, qr, svd
using VlasovMethods: VlasovMethods

# F(f) = [ϕ ⊗ 1 + 1 ⊗ v²/2, f] / (hx hv), with ϕ = G ρ(f), on the formed nx × nv array
function dense_rhs(d, f)
    ḟ = similar(vec(f))
    VlasovMethods._grid_1d1v_rhs!(ḟ, zero(eltype(f)), vec(f), d.p)
    return reshape(ḟ, size(f))
end

dense_kfield(d, K, V) = d.hv .* (dense_rhs(d, K * V') * V)
dense_lfield(d, X, L) = d.hx .* (dense_rhs(d, X * L')' * X)
dense_sfield(d, X, S, V) = (d.hx * d.hv) .* (X' * dense_rhs(d, X * S * V') * V)

# one step of the classical RK4 tableau, b = (1, 2, 2, 1)/6
function dense_rk4(g, y, τ)
    k1 = g(y)
    k2 = g(y .+ (τ / 2) .* k1)
    k3 = g(y .+ (τ / 2) .* k2)
    k4 = g(y .+ τ .* k3)
    return y .+ (τ / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

# the columns of Q / √h of the Householder QR of √h A: h QᵀQ = I, min(size(A)...) columns
dense_basis(A, h) = Matrix(qr(sqrt(h) .* A).Q) ./ sqrt(h)

# fixed-rank BUG
function dense_bug(d, X, S, V, τ)
    K = dense_rk4(K -> dense_kfield(d, K, V), X * S, τ)
    L = dense_rk4(L -> dense_lfield(d, X, L), V * S', τ)
    X₁, V₁ = dense_basis(K, d.hx), dense_basis(L, d.hv)
    M, N = d.hx .* (X₁' * X), d.hv .* (V₁' * V)
    S₁ = dense_rk4(S -> dense_sfield(d, X₁, S, V₁), M * S * N', τ)
    return X₁, S₁, V₁
end

# the augmented BUG step, untruncated
function dense_augmented(d, X, S, V, τ)
    K = dense_rk4(K -> dense_kfield(d, K, V), X * S, τ)
    L = dense_rk4(L -> dense_lfield(d, X, L), V * S', τ)
    X̂, V̂ = dense_basis(hcat(K, X), d.hx), dense_basis(hcat(L, V), d.hv)
    Ŝ₀ = (d.hx .* (X̂' * X)) * S * (d.hv .* (V' * V̂))
    return X̂, dense_rk4(S -> dense_sfield(d, X̂, S, V̂), Ŝ₀, τ), V̂
end

# the least r with ‖σ[r+1:end]‖ ≤ tol, clamped to [rmin, min(rmax, length(σ))]
function dense_rank(σ, tol, rmin, rmax)
    k = length(σ)
    r = findfirst(r -> norm(σ[(r + 1):k]) ≤ tol, 0:k) - 1
    return clamp(r, rmin, min(rmax, k))
end

function dense_truncate(X̂, Ŝ, V̂, r₁)
    F = svd(Ŝ)
    return X̂ * F.U[:, 1:r₁], Matrix(Diagonal(F.S[1:r₁])), V̂ * F.V[:, 1:r₁]
end

# the augmented BUG, truncated to the rank of its input (`tol === nothing`) or by `tol`
function dense_augmented_bug(d, X, S, V, τ; tol = nothing, rmin = 1, rmax = size(S, 1))
    X̂, Ŝ, V̂ = dense_augmented(d, X, S, V, τ)
    σ = svd(Ŝ).S
    r₁ = tol === nothing ? min(size(S, 1), length(σ)) : dense_rank(σ, tol, rmin, rmax)
    return dense_truncate(X̂, Ŝ, V̂, r₁)
end

# the midpoint BUG (4r)
function dense_midpoint_bug(d, X, S, V, τ)
    Û, Ŝ½, V̂ = dense_augmented(d, X, S, V, τ / 2)
    Ŷ½ = Û * Ŝ½ * V̂'
    FV = d.hv .* (dense_rhs(d, Ŷ½) * V̂)
    FU = d.hx .* (dense_rhs(d, Ŷ½)' * Û)
    Ū, V̄ = dense_basis(hcat(Û, τ .* FV), d.hx), dense_basis(hcat(V̂, τ .* FU), d.hv)
    S̄₀ = (d.hx .* (Ū' * X)) * S * (d.hv .* (V' * V̄))
    S̄ = dense_rk4(S -> dense_sfield(d, Ū, S, V̄), S̄₀, τ)
    return dense_truncate(Ū, S̄, V̄, min(size(S, 1), minimum(size(S̄))))
end
