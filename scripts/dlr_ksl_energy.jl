# Probe: projector-splitting (Lie, K–S–L) on the Arakawa grid, each substep by the implicit
# midpoint rule. Variant A uses the projected Hamiltonian gradient Q h in each substep
# (a skew-gradient system on the substep's linear space); variant B uses h itself (the
# Galerkin field of standard DLRA). Claim: A keeps H and the L² norm to the solver tolerance
# for any basis; B does not keep H. Also: the Arakawa form T(a,b,c) = ⟨a, J(b,c)⟩ is
# totally antisymmetric.
# Run, from the repository root:
#   julia --startup-file=no --project=scripts -e 'r = include("scripts/dlr_ksl_energy.jl"); println(r)'
using LinearAlgebra, Random
include("dlr_energy_tangency.jl")          # arakawa(f, h, hx, hv), orth(A, h)

function setup(; nx = 16, nv = 32, α = 0.3, k = 0.5)
    L = 2π / k
    hx = L / nx
    x = (0:(nx - 1)) .* hx
    v = collect(range(-6, 6; length = nv))
    hv = v[2] - v[1]
    # Green's matrix of -φ'' = ρ - mean(ρ), spectral, symmetric: the conjugate pairs of modes
    # with the weight 2/nx, and for even nx the Nyquist mode with the weight 1/nx
    G = zeros(nx, nx)
    for a in 1:nx, b in 1:nx, m in 1:(nx ÷ 2)
        w = 2m == nx ? 1 / nx : 2 / nx
        G[a, b] += w * cos(m * k * hx * (a - b)) / (m * k)^2
    end
    f0 = (1 .+ α .* cos.(k .* x)) .* (exp.(-v .^ 2 ./ 2) ./ sqrt(2π))'
    (; nx, nv, hx, hv, v, G, f0)
end

ρ(d, f) = d.hv .* vec(sum(f; dims = 2))
hgrad(d, f) = (d.G * ρ(d, f)) .+ (d.v .^ 2 ./ 2)'
function energy(d, f)
    d.hx * d.hv * sum((d.v .^ 2 ./ 2)' .* f) + d.hx / 2 * dot(ρ(d, f), d.G * ρ(d, f))
end
l2(d, f) = d.hx * d.hv * sum(abs2, f)
J(d, f, h) = -arakawa(f, h, d.hx, d.hv)    # ḟ = {h, f}, sign immaterial for the claims

# implicit midpoint on y' = F(y) by fixed-point iteration to round-off
function midpoint(F, y0, τ; tol = 1e-15, maxit = 200)
    y1 = y0 + τ * F(y0)
    for _ in 1:maxit
        y = y0 + τ * F((y0 + y1) / 2)
        norm(y - y1) ≤ tol * norm(y) && return y
        y1 = y
    end
    error("midpoint iteration did not converge")
end

PV(d, V, g) = g * (d.hv .* V * V')         # projector onto L²_x ⊗ span V
PX(d, X, g) = (d.hx .* X * X') * g         # projector onto span X ⊗ L²_v

function ksl_step(d, X, S, V, τ; projected)
    K = midpoint(K -> begin
            f = K * V'
            h = hgrad(d, f)
            J(d, f, projected ? PV(d, V, h) : h) * (d.hv .* V)
        end, X * S, τ)
    X = orth(K, d.hx)
    S = d.hx .* X' * K
    S = midpoint(S -> begin
            f = X * S * V'
            h = hgrad(d, f)
            hq = projected ? PX(d, X, PV(d, V, h)) : h
            d.hx * d.hv .* X' * J(d, f, hq) * V
        end, S, -τ)
    Lm = midpoint(Lm -> begin
            f = X * Lm'
            h = hgrad(d, f)
            (d.hx .* X' * J(d, f, projected ? PX(d, X, h) : h))'
        end, V * S', τ)
    V = orth(Lm, d.hv)
    S = (d.hv .* V' * Lm)'
    return X, S, V
end

function run(; r = 4, steps = 40, τ = 0.05, projected)
    d = setup()
    F = svd(sqrt(d.hx * d.hv) .* d.f0)
    X = F.U[:, 1:r] ./ sqrt(d.hx)
    V = F.V[:, 1:r] ./ sqrt(d.hv)
    # complete the rank-one initial state with random directions
    rng = Xoshiro(3)
    X = orth(hcat(X[:, 1:2], randn(rng, d.nx, r - 2)), d.hx)
    V = orth(hcat(V[:, 1:1], randn(rng, d.nv, r - 1)), d.hv)
    S = d.hx * d.hv .* X' * d.f0 * V
    f = X * S * V'
    H0, C0 = energy(d, f), l2(d, f)
    dH = dC = 0.0
    for _ in 1:steps
        X, S, V = ksl_step(d, X, S, V, τ; projected)
        f = X * S * V'
        dH = max(dH, abs(energy(d, f) - H0) / H0)
        dC = max(dC, abs(l2(d, f) - C0) / C0)
    end
    (energy_drift = dH, l2_drift = dC)
end

function antisymmetry(; seed = 7)
    d = setup()
    rng = Xoshiro(seed)
    a, b, c = (randn(rng, d.nx, d.nv) for _ in 1:3)
    T(a, b, c) = d.hx * d.hv * sum(a .* J(d, b, c))
    s = abs(T(a, b, c))
    (T_abc = s, swap_bc = abs(T(a, b, c) + T(a, c, b)) / s,
        swap_ab = abs(T(a, b, c) + T(b, a, c)) / s)
end

(antisymmetry = antisymmetry(), A_projected = run(projected = true),
    B_plain = run(projected = false))
