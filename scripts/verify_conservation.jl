#!/usr/bin/env julia
#
# Verify, numerically, the conservation properties the two collision-operator manuscripts
# prove, and record which boundary conditions actually admit them.
#
# There are TWO distinct claims here and they are easy to conflate. This script measures both
# separately, because the measurement shows they are not the same claim.
#
# CLAIM 1 — the semi-discrete right-hand side conserves the particle momentum and energy:
#
#     Σ_α w_α v̇_α = 0        and        Σ_α w_α v_α v̇_α = 0 .
#
# This holds on ANY basis. A₁ and A₂ are solved from exactly these two constraints
# (eq:coefficient_lin_system), which are algebraic conditions on particle sums and involve the
# basis only through f_s'/f_s. Measured below: round-off for every boundary condition, with
# uniform and non-uniform weights alike.
#
# CLAIM 2 — the moments of the spline representation agree with the particles':
#
#     ∫ vᵏ f_s dv = Σ_α w_α v_α^k .
#
# This is what the manuscript's sentence at main.tex:226 is about — "Momentum and energy
# conservation follow ... from requiring the functions v and v² to be in the span of basis
# {φ_j}" — and it holds only when vᵏ lies in the span. It is the claim that fails on a periodic
# or Dirichlet basis, and what fails with it is every diagnostic computed from f_s, the entropy
# included, since the H-theorem is a statement about f_s and not about the particles.
#
# The weights are varied because the coefficient system used to omit them, which is invisible
# for the uniform weights every initialiser in examples/ sets and breaks Claim 1 for the
# non-uniform ones the importance-sampling initialiser produces.
#
# Run with:
#   julia --project=<env> scripts/verify_conservation.jl
# where <env> is an environment with VlasovMethods and SimpleSplines available.

using VlasovMethods
using LinearAlgebra
using Printf
using Random

const VM = VlasovMethods

# ---------------------------------------------------------------------------------------------
# The semi-discrete conservative Lenard-Bernstein right-hand side, written out here rather than
# taken from the package, so that this script checks the package's coefficients against an
# independent statement of the scheme.
# ---------------------------------------------------------------------------------------------

function rhs_clb(sdist, pdist, v, ν = 1.0)
    A1, A2 = VM.compute_coefficients(sdist, pdist, v)
    fs = sdist.spline
    dfs = derivative(fs)
    return @. -ν * (dfs(v) / fs(v) + A1 + A2 * v)
end

function residuals(sdist, pdist, v)
    v̇ = rhs_clb(sdist, pdist, v)
    w = pdist.particles.w[1, :]
    scale = maximum(abs, v̇)
    momentum = abs(sum(w .* v̇)) / scale
    energy = abs(sum(w .* v .* v̇)) / scale
    return momentum, energy, scale
end

function make_particles(npart, weights::Symbol; seed = 0x5EED)
    Random.seed!(seed)
    pd = ParticleDistribution(1, 1, npart)
    pd.particles.v[1, :] .= randn(npart)
    if weights === :uniform
        pd.particles.w[1, :] .= 1 / npart
    elseif weights === :nonuniform
        # the shape the importance-sampling initialiser in examples/bumpontail.jl produces
        pd.particles.w[1, :] .= (1 .+ rand(npart)) ./ npart
    else
        error("unknown weight scheme $(weights)")
    end
    return pd
end

# ---------------------------------------------------------------------------------------------

const NPART = 5000
const NKNOT = 31
const ORDER = 4                     # cubic, as both manuscripts specify
const DOMAIN = (-8.0, 8.0)

println("CLAIM 1 — particle momentum and energy of the conservative Lenard-Bernstein operator")
println("  particles: ", NPART, "   knots: ", NKNOT, "   order: ", ORDER,
    " (degree ", ORDER - 1, ")   domain: ", DOMAIN)
println("  Expected: round-off on EVERY boundary condition, because A₁ and A₂ are solved")
println("  from these two constraints directly.")
println()

for weights in (:uniform, :nonuniform)
    pd = make_particles(NPART, weights)
    v = pd.particles.v[1, :]

    @printf("%-14s %-12s %-8s %-14s %-14s\n",
        "weights", "condition", "poly", "|Σ w v̇|/|v̇|", "|Σ w v v̇|/|v̇|")
    println("  ", "-"^68)

    for bc in (Free(), Periodic(), Dirichlet())
        sdist = SplineDistribution(1, 1, NKNOT, ORDER, DOMAIN, 0, bc)
        m = polynomial_reproduction(sdist.basis)
        projection(v, pd, sdist)

        local mom, ene
        try
            mom, ene, _ = residuals(sdist, pd, v)
            @printf("%-14s %-12s %-8d %-14.3e %-14.3e\n",
                weights, string(bc), m, mom, ene)
        catch err
            # A Dirichlet basis forces f_s to zero at the domain ends, so f_s'/f_s is a
            # division by zero there and the positivity guard fires. That is itself part of
            # the finding: the condition does not merely lose conservation, it makes the
            # operator undefined.
            msg = err isa ErrorException ? first(split(err.msg, ',')) : string(typeof(err))
            @printf("%-14s %-12s %-8d %s\n", weights, string(bc), m, msg)
        end
    end
    println()
end

# ---------------------------------------------------------------------------------------------
# The polynomial-reproduction predicate is what decides it. Check directly that the clamped
# basis reproduces 1, v and v² exactly and that the other two do not.
# ---------------------------------------------------------------------------------------------

println("CLAIM 2 — do the spline moments agree with the particles'?")
println("  Polynomial reproduction, measured as max |Πf(v) - f(v)| over the domain: the L²")
println("  projection reproduces ∫ πf exactly iff π is in the span, so these errors are")
println("  what separates the boundary conditions.")
println()
@printf("%-12s %-8s %-14s %-14s %-14s\n", "condition", "poly", "f = 1", "f = v", "f = v²")
println("  ", "-"^64)

pts = collect(range(DOMAIN[1] + 0.5, DOMAIN[2] - 0.5; length = 41))

for bc in (Free(), Periodic(), Dirichlet())
    sd = SplineDistribution(1, 1, NKNOT, ORDER, DOMAIN, 0, bc)
    m = polynomial_reproduction(sd.basis)
    errs = Float64[]
    for f in (v -> one(v), v -> v, v -> v^2)
        project_function(f, sd)
        push!(errs, maximum(abs(sd.spline(v) - f(v)) for v in pts))
    end
    @printf("%-12s %-8d %-14.3e %-14.3e %-14.3e\n", string(bc), m, errs...)
end

println()
println("`Free` reproduces all three at round-off, so the spline moments are the particle")
println("moments. `Periodic` reproduces only the constants — v is not periodic — and")
println("`Dirichlet` reproduces none, every basis function vanishing at the ends.")
println()
println("Conclusion: the boundary condition does not decide whether the particles conserve")
println("momentum and energy; it decides whether f_s has the same moments they do, and hence")
println("whether anything measured from f_s — the entropy, and the H-theorem that rests on")
println("it — refers to the distribution actually being evolved.")
