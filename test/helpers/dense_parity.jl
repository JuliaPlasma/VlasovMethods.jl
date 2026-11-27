# The low-rank parity splitting transcribed on the formed grid function, the oracle of
# `test/lowrank/parity.jl` and the explicit control of `test/verification/lowrank_parity.jl`. The
# field is the grid model's bracket term `P(f, h) = poisson_apply(Arakawa(nx, nv, hx, hv), f, h) /
# (hx hv)` with `h = G ρ ⊗ 1 + 1 ⊗ v²/2`, applied to the formed `nx × nv` arrays; every projection,
# augmentation and truncation acts on the factors. The state is `f = 1 ⊗ f0 + f⁺ + f⁻`, with
# `f⁺ = X⁺ S⁺ V⁺ᵀ` (every column of `V⁺` even, the first two spanning `{1, v²}`) and
# `f⁻ = X⁻ S⁻ V⁻ᵀ` (every column odd, the first one `∝ v`). One step of `τ` is the `f⁻` substep
# of `τ/2` (the exact low-rank increment, then the truncation that keeps the fixed modes), the
# `f⁺` step of `τ` (the augmented BUG step with `f⁻` frozen, its S-step on the closed space
# `[X̂, 1, G X̂]`, every substep an implicit midpoint step by `iters` fixed-point iterations), and
# the `f⁻` substep of `τ/2` again.
#
# An object is a named tuple `(X, S, V, m)` with `m` fixed modes.

using GeometricBrackets: Arakawa, poisson_apply
using LinearAlgebra: I, norm, qr, svd

# the data of the case: the grid of `space`, the Green's matrix `G` and the background `f0`
function dense_parity_case(space, G, f0)
    T = eltype(f0)
    nx, nv = length(space.x), length(space.v)
    hx, hv = step(space.x), step(space.v)
    v = collect(space.v)
    a = Arakawa(nx, nv, hx, hv)
    P(f, h) = reshape(poisson_apply(a, vec(f), vec(h)), nx, nv) ./ (hx * hv)
    kin = v .^ 2 ./ 2
    ρ(g) = hv .* vec(sum(g; dims = 2))
    φ(g) = G * ρ(g)
    hof(ϕ) = ϕ .+ kin'
    F0 = ones(T, nx) .* f0'
    return (; T, nx, nv, hx, hv, v, kin, G, P, ρ, φ, hof, f0, F0)
end

dense_reflectv(V) = V[end:-1:1, :]

# the thin weighted QR in x; as many columns as nodes or more give the whole space
function dense_orth_x(d, A)
    size(A, 2) ≥ d.nx && return Matrix{d.T}(I, d.nx, d.nx) ./ sqrt(d.hx)
    return Matrix(qr(sqrt(d.hx) .* A).Q) ./ sqrt(d.hx)
end

dense_qr_x(d, Y) = (F = qr(sqrt(d.hx) .* Y); (Matrix(F.Q) ./ sqrt(d.hx), Matrix(F.R)))

# [Vfix Vnew]: an orthonormal basis of span [Vfix, A] in the subspace of the given parity (+1
# even, -1 odd), by the weighted Householder QR on the half grid, so that every column has its
# parity bit for bit; Vfix is kept as the leading columns
function dense_orth_v(d, Vfix, A, parity)
    T, nv, hv = d.T, d.nv, d.hv
    nh = parity > 0 ? cld(nv, 2) : fld(nv, 2)
    w = [(parity > 0 && isodd(nv) && j == nh) ? one(T) : T(2) for j in 1:nh]
    C = hcat(Vfix, (A .+ parity .* dense_reflectv(A)) ./ 2)
    size(C, 2) == 0 && return C
    s = sqrt.(hv .* w)
    Q = Matrix(qr(s .* C[1:nh, :]).Q) ./ s
    V = zeros(T, nv, size(Q, 2))
    for j in 1:nh
        V[j, :] .= Q[j, :]
        nv + 1 - j != j && (V[nv + 1 - j, :] .= parity .* Q[j, :])
    end
    return hcat(Vfix, V[:, (size(Vfix, 2) + 1):end])
end

# the least r with ‖σ[r+1:end]‖ ≤ abstol, clamped to [rlo, rcap] and to the number of values
function dense_trunc_rank(σ, abstol, rlo, rcap)
    n = length(σ)
    tail = sqrt.(reverse(cumsum(reverse(σ .^ 2))))
    rr = findfirst(j -> (j + 1 > n ? zero(eltype(σ)) : tail[j + 1]) ≤ abstol, 0:n)
    rleast = rr === nothing ? n : rr - 1
    return clamp(rleast, min(max(rlo, 0), n), min(rcap, n))
end

# f = K Vᵀ with the fixed modes V[:, 1:m]: K[:, 1:m] V[:, 1:m]ᵀ is kept, and the rest is truncated
# by its weighted SVD to r - m columns, r = clamp(r*, rmin, rmax), with the tolerance relative to
# ‖f‖ = √hx ‖K‖_F
function dense_truncate(d, K, V, m; tol, rmax, rmin = rmax)
    nrm = sqrt(d.hx) * norm(K)
    KM, K⊥ = K[:, 1:m], K[:, (m + 1):end]
    if size(K⊥, 2) == 0
        X, S = dense_qr_x(d, KM)
        return (X = X, S = S, V = V[:, 1:m], m = m)
    end
    F = svd(sqrt(d.hx) .* K⊥)
    r = dense_trunc_rank(F.S, tol * nrm, rmin - m, rmax - m)
    X, S = dense_qr_x(d, hcat(KM, K⊥ * F.V[:, 1:r]))
    return (X = X, S = S, V = hcat(V[:, 1:m], V[:, (m + 1):end] * F.V[:, 1:r]), m = m)
end

dense_object(o) = o.X * o.S * o.V'

# the implicit midpoint step y1 = y0 + τ R((y0 + y1)/2) by `iters` fixed-point iterations from
# y1 = y0; `iters = 1` is explicit Euler
function dense_midpoint(R, y0, τ; iters = 2)
    y1 = copy(y0)
    for _ in 1:iters
        y1 = y0 .+ τ .* R((y0 .+ y1) ./ 2)
    end
    return y1
end

# the f⁺ step of τ, f⁻ frozen
function dense_fplus_step(d, op, om, τ; tol, rmax, rmin = rmax, iters = 2)
    fm = dense_object(om)
    X0, S0, V0, m = op.X, op.S, op.V, op.m
    field(fp) = d.P(fm, d.hof(d.φ(d.F0 .+ fp)))
    K1 = dense_midpoint(K -> field(K * V0') * (d.hv .* V0), X0 * S0, τ; iters)
    L1 = dense_midpoint(L -> field(X0 * L')' * (d.hx .* X0), V0 * S0', τ; iters)
    X̂ = dense_orth_x(d, hcat(K1, X0))
    V̂ = dense_orth_v(d, V0[:, 1:m], hcat(V0[:, (m + 1):end], L1), +1)
    Xa = dense_orth_x(d, hcat(X̂, ones(d.T, d.nx), d.G * X̂))
    Ŝ0 = (d.hx .* (Xa' * X0)) * S0 * (d.hv .* (V0' * V̂))
    S1 = dense_midpoint(S -> d.hx * d.hv .* (Xa' * field(Xa * S * V̂') * V̂), Ŝ0, τ; iters)
    return dense_truncate(d, Xa * S1, V̂, m; tol, rmax, rmin)
end

# the v-range of P(g, h) for g = Σ X_α ⊗ V_α and h = ϕ ⊗ 1 + 1 ⊗ v²/2: the periodic shifts S_t V
# and (S_q k) ∘ (S_t V), q, t ∈ {-1, 0, 1}
function dense_vrange(d, Vg)
    sh(A, t) = circshift(A, (t, 0))
    return reduce(hcat,
        [q == 2 ? sh(Vg, t) : circshift(d.kin, q) .* sh(Vg, t) for t in -1:1 for q in -1:2])
end

# the f⁻ substep of length τ: f⁻ + τ P(F0 + f⁺, h(F0 + f⁺)), the increment in the odd part of
# the v-range, then truncated
function dense_fminus_step(d, op, om, τ; tol, rmax, rmin = rmax)
    g = d.F0 .+ dense_object(op)
    Vinc = dense_orth_v(d, zeros(d.T, d.nv, 0), dense_vrange(d, hcat(d.f0, op.V)), -1)
    Kinc = d.P(g, d.hof(d.φ(g))) * (d.hv .* Vinc)
    V = dense_orth_v(d, om.V[:, 1:(om.m)], hcat(om.V[:, (om.m + 1):end], Vinc), -1)
    K = om.X * om.S * (d.hv .* (om.V' * V)) .+ τ .* Kinc * (d.hv .* (Vinc' * V))
    return dense_truncate(d, K, V, om.m; tol, rmax, rmin)
end

# one step of τ, with the ranks `rmax`, least ranks `rmin` and tolerances `tol` of (f⁺, f⁻)
function dense_parity_step(d, op, om, τ; tol, rmax, rmin = rmax, iters = 2)
    om = dense_fminus_step(d, op, om, τ / 2; tol = tol[2], rmax = rmax[2], rmin = rmin[2])
    op = dense_fplus_step(d, op, om, τ;
        tol = tol[1], rmax = rmax[1], rmin = rmin[1], iters)
    om = dense_fminus_step(d, op, om, τ / 2; tol = tol[2], rmax = rmax[2], rmin = rmin[2])
    return op, om
end

# the objects of a `LowRankParityDistribution`
function dense_objects(Y)
    ((X = copy(Y.plus.X), S = copy(Y.plus.S), V = copy(Y.plus.V), m = 2),
        (X = copy(Y.minus.X), S = copy(Y.minus.S), V = copy(Y.minus.V), m = 1))
end
