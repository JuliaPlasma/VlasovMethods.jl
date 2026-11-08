# The stalled steps of the rate check of `lenard_bernstein_metriplectic.jl` at the shifted
# normal's parameters (1000 particles, M = 41, h = 8e-4, t = 1, ν = 1), for each bracket.
#
#   julia --startup-file=no --project=scripts scripts/lenard_bernstein_metriplectic_rate_stalls.jl
#
# It repeats the loop of `relax_frozen` step by step and, for each step whose damped fixed-point
# solve stops at the iteration cap, prints the step, the norm of the residual, the particle with the
# largest residual component and its velocity (last iterate x), the mean and the standard deviation
# of x, and the count of particles of x within ±0.5 of that velocity (the particle itself included).
# It also calls the driver's own `relax_frozen` and checks that its `stalled` and `residual` agree.
# `KNOWN_ISSUES.md` K95 cites these figures.

module MPLB
include(joinpath(@__DIR__, "lenard_bernstein_metriplectic.jl"))
end

using LinearAlgebra, Printf
using .MPLB: CASES, BRACKETS, RATE_TEMPERATURE, RATE_SOLVER_ITERATIONS,
             RATE_SOLVER_TOLERANCE, RATE_SOLVER_DAMPING, frozen_rhs!, temperature,
             relax_frozen, initial_distribution, spline_distribution, nsteps,
             MetriplecticLenardBernstein, CollisionEntropy

const CASE = CASES["shifted-normal"]
const OPTIONS = CASE.paper

function probe(initial, bracket, options)
    dist = deepcopy(initial)
    w = dist.particles.w[1, :]
    v₀ = dist.particles.v[1, :]
    n₀ = sum(w)
    u₀ = dot(w, v₀) / n₀
    T₀ = RATE_TEMPERATURE
    model = MetriplecticLenardBernstein(dist,
        CollisionEntropy(spline_distribution(options.nbasis)); ν = options.nu, bracket)
    n = nsteps(options.final_time, options.tstep)
    h = options.tstep
    v = u₀ .+ (v₀ .- u₀) .* sqrt(T₀ / (2 * temperature(v₀, w)))
    v̇ = similar(v)
    r = similar(v)
    stalled, residual = 0, 0.0
    for k in 1:n
        x = copy(v)
        converged = false
        for _ in 1:RATE_SOLVER_ITERATIONS
            frozen_rhs!(v̇, (x .+ v) ./ 2, model, n₀, u₀, T₀ + u₀^2)
            r .= v .+ h .* v̇ .- x
            converged = norm(r) ≤ RATE_SOLVER_TOLERANCE * norm(x)
            converged && break
            x .+= RATE_SOLVER_DAMPING .* r
        end
        if !converged
            stalled += 1
            residual = max(residual, norm(r))
            i = argmax(abs.(r))
            μ = dot(w, x) / sum(w)
            σ = sqrt(temperature(x, w))
            near = count(y -> abs(y - x[i]) ≤ 0.5, x)
            @printf("  step %d (t = %.4f): ‖r‖ = %.4e, largest |rᵢ| = %.3e at particle %d, v = %.4f (vⁿ %.4f); mean %.4f, std %.4f; %d particles within ±0.5\n",
                k, k * h, norm(r), abs(r[i]), i, x[i], v[i], μ, σ, near)
        end
        v = x
    end
    return (; stalled, residual, n)
end

for name in ("dougherty", "jeyakumar")
    bracket = BRACKETS[name]
    initial = initial_distribution(CASE.example, OPTIONS.npart)
    println(name)
    p = probe(initial, bracket, OPTIONS)
    d = relax_frozen(initial, bracket, OPTIONS)
    @printf("  probe: %d of %d steps stalled, largest residual %.4e\n", p.stalled, p.n,
        p.residual)
    @printf("  relax_frozen: %d stalled, largest residual %.4e, agrees: %s\n",
        d.stalled, d.residual, d.stalled == p.stalled && d.residual == p.residual)
    d.stalled == p.stalled && d.residual == p.residual ||
        error("$name: the probe and relax_frozen disagree")
end
