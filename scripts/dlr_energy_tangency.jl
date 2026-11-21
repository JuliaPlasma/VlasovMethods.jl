# Probe: on a doubly periodic Arakawa grid, the DLR field P(f) J(f, h) keeps the energy
# H = hx hv Σ (v²/2) f + ½ ρᵀ G ρ hx when 1 and v² lie in span(V), and the L² norm always.
# Control: without 1, v² in span(V) the energy rate is nonzero.
# Run, from the repository root:
#   julia --startup-file=no --project=scripts -e 'r = include("scripts/dlr_energy_tangency.jl"); println(r)'
using LinearAlgebra, Random

function arakawa(f, h, hx, hv)
    nx, nv = size(f)
    J = zeros(nx, nv)
    ip(i) = mod1(i + 1, nx)
    im(i) = mod1(i - 1, nx)
    jp(j) = mod1(j + 1, nv)
    jm(j) = mod1(j - 1, nv)
    for i in 1:nx, j in 1:nv

        J1 = (f[ip(i), j] - f[im(i), j]) * (h[i, jp(j)] - h[i, jm(j)]) -
             (f[i, jp(j)] - f[i, jm(j)]) * (h[ip(i), j] - h[im(i), j])
        J2 = f[ip(i), j] * (h[ip(i), jp(j)] - h[ip(i), jm(j)]) -
             f[im(i), j] * (h[im(i), jp(j)] - h[im(i), jm(j)]) -
             f[i, jp(j)] * (h[ip(i), jp(j)] - h[im(i), jp(j)]) +
             f[i, jm(j)] * (h[ip(i), jm(j)] - h[im(i), jm(j)])
        J3 = f[ip(i), jp(j)] * (h[i, jp(j)] - h[ip(i), j]) -
             f[im(i), jm(j)] * (h[im(i), j] - h[i, jm(j)]) -
             f[im(i), jp(j)] * (h[i, jp(j)] - h[im(i), j]) +
             f[ip(i), jm(j)] * (h[ip(i), j] - h[i, jm(j)])
        J[i, j] = (J1 + J2 + J3) / (12hx * hv)
    end
    return J
end

orth(A, h) = Matrix(qr(sqrt(h) .* A).Q)[:, 1:size(A, 2)] ./ sqrt(h)

function rates(; nx = 16, nv = 24, r = 4, fixed = true, seed = 1)
    rng = Xoshiro(seed)
    L = 4π
    hx = L / nx
    v = range(-6, 6; length = nv)
    hv = step(v)
    k = 2π / L
    # Green's matrix of -φ'' = ρ - mean(ρ), spectral, symmetric
    G = zeros(nx, nx)
    for a in 1:nx, b in 1:nx, m in 1:(nx ÷ 2 - 1)
        G[a, b] += 2 / nx * cos(m * k * hx * (a - b)) / (m * k)^2
    end
    X = orth(randn(rng, nx, r), hx)
    Vr = randn(rng, nv, r)
    fixed && (Vr[:, 1] .= 1; Vr[:, 2] .= v .^ 2)
    V = orth(Vr, hv)
    S = randn(rng, r, r)
    f = X * S * V'
    ρ = hv .* vec(sum(f; dims = 2))
    φ = G * ρ
    h = φ .+ (v .^ 2 ./ 2)'                # δH/δf on the grid, in the weighted product
    F = arakawa(f, h, hx, hv)              # ḟ = J(f, h) up to sign
    PX = X * (hx .* X')                    # orthogonal projectors in the grid products
    PV = V * (hv .* V')
    PF = PX * F + F * PV' - PX * F * PV'
    ip(a, b) = hx * hv * sum(a .* b)
    (energy_rate = ip(h, PF), l2_rate = ip(f, PF), mass_rate = ip(ones(nx, nv), PF),
        scale = sqrt(ip(h, h) * ip(PF, PF)), full_energy_rate = ip(h, F))
end

with = rates(fixed = true)
without = rates(fixed = false)
(with = map(x -> round(x; sigdigits = 3), with),
    without = map(x -> round(x; sigdigits = 3), without))
