# Probe: nonlinear Vlasov–Poisson in parity variables f = f0 + f⁺ + f⁻ (f0 even and stationary)
# on the Arakawa grid with a velocity grid symmetric about 0.
#  1. Grading: J(a, b) has parity −p_a p_b under v ↦ −v.
#  2. H depends on f⁺ only: ρ(f⁻) = 0 and Σ v² f⁻ = 0.
#  3. Parity Strang splitting: f⁻ += τ/2 J(f0 + f⁺, h(f⁺)) (exact: f⁺ frozen), then the implicit
#     midpoint step of ḟ⁺ = J(f⁻, h(f⁺)) with f⁻ frozen (a skew-gradient system in f⁺), then f⁻
#     again. Claim: H is kept to the solver tolerance; the error is of second order in τ.
#  Control: the same splitting with an explicit Euler f⁺ step does not keep H.
# Run, from the repository root:
#   julia --startup-file=no --project=scripts -e 'r = include("scripts/parity_splitting.jl"); println(r)'
using LinearAlgebra
include("dlr_ksl_energy.jl")    # setup, ρ, hgrad, energy, J, midpoint

reflect(g) = g[:, end:-1:1]
even(g) = (g + reflect(g)) / 2
odd(g) = (g - reflect(g)) / 2

function grading(d; seed = 11)
    rng = Xoshiro(seed)
    a⁺, b⁺ = even(randn(rng, d.nx, d.nv)), even(randn(rng, d.nx, d.nv))
    a⁻ = odd(randn(rng, d.nx, d.nv))
    n(g) = norm(g) / norm(J(d, a⁺, b⁺))
    (even_even_has_even_part = n(even(J(d, a⁺, b⁺))),
        even_odd_has_odd_part = n(odd(J(d, a⁺, a⁻))))
end

function parity_step(d, f⁺, f⁻, τ; implicit = true)
    f0 = d.f0h
    f⁻ = f⁻ + τ / 2 * J(d, f0 + f⁺, hgrad(d, f0 + f⁺))
    F(g) = J(d, f⁻, hgrad(d, f0 + g))
    f⁺ = implicit ? midpoint(F, f⁺, τ) : f⁺ + τ * F(f⁺)
    f⁻ = f⁻ + τ / 2 * J(d, f0 + f⁺, hgrad(d, f0 + f⁺))
    return f⁺, f⁻
end

rk4step(d, f, τ) = begin
    F(g) = J(d, g, hgrad(d, g))
    k1 = F(f)
    k2 = F(f + τ / 2 * k1)
    k3 = F(f + τ / 2 * k2)
    k4 = F(f + τ * k3)
    f + τ / 6 * (k1 + 2k2 + 2k3 + k4)
end

function run_parity(d, T, τ; implicit = true)
    f⁺ = d.f0 - d.f0h
    f⁻ = zero(f⁺)      # initial perturbation α cos kx f_M is even
    H0 = energy(d, d.f0h + f⁺ + f⁻)
    dH = 0.0
    for _ in 1:round(Int, T / τ)
        f⁺, f⁻ = parity_step(d, f⁺, f⁻, τ; implicit)
        dH = max(dH, abs(energy(d, d.f0h + f⁺ + f⁻) - H0) / H0)
    end
    (f = d.f0h + f⁺ + f⁻, dH = dH)
end

d0 = setup(; nx = 16, nv = 32, α = 0.2)
# the homogeneous even background f0(v) = mean over x of the initial state
d = merge(d0, (f0h = repeat(sum(d0.f0; dims = 1) ./ d0.nx, d0.nx, 1),))
T = 2.0
ref = let f = d.f0, τ = 1e-3
    for _ in 1:round(Int, T / τ)
        f = rk4step(d, f, τ)
    end
    f
end
errs = map((0.1, 0.05, 0.025)) do τ
    r = run_parity(d, T, τ)
    (τ = τ, err = norm(r.f - ref) / norm(ref - d.f0h), dH = r.dH)
end
ctrl = run_parity(d, T, 0.05; implicit = false)
(grading = grading(d),
    f⁻_moments = let g = odd(randn(Xoshiro(5), d.nx, d.nv))
        (ρ = norm(ρ(d, g)), v² = abs(sum((d.v .^ 2)' .* g)))
    end,
    runs = errs, orders = [log2(errs[i].err / errs[i + 1].err) for i in 1:2],
    control_explicit_dH = ctrl.dH)
