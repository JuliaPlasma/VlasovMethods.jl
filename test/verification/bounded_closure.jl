# Runs of the grid 1D1V Vlasov–Poisson model with the bounded velocity closure, through the
# model's `GeometricIntegrator` and its default `ImplicitMidpoint()`.
#
# Run A needs `f` at the boundary of the velocity domain, `[-3.5, 3.5]`: the bounded closure differs
# from the periodic one only where `f` reaches the end nodes, and a run whose `f` vanishes there
# shows neither the boundary flux nor the failure of the control. Its discrete Hamiltonian
# `H = hx hv Σ_i Σ_{j=2}^{nv-1} f_ij v_j²/2 + (hx/2) ρᵀ G ρ` is quadratic, so the implicit midpoint
# keeps it to the Newton solve, and the mass change of each step is `τ` times the boundary flux
# `−hx hv Σ_i Σ_{j ∈ {1, nv}} [Ph, χ ∘ f]_ij` at the midpoint state. The control brackets with the
# Hamiltonian from its formula on every node, the uncut `h`. Run B keeps `f` round-off at `±vmax`,
# where the bounded and the periodic closure give the same run. H, the mass and the flux are
# computed here and not in `src`.

using GeometricBrackets: poisson_apply
using GeometricEquations: GeometricEquations
# an explicit only-list, and not `using GeometricIntegrators`: both it and `VlasovMethods` export
# `GeometricIntegrator`, and the only-list binds it to this one
using GeometricIntegrators: ImplicitMidpoint, integrate, ntime
using LinearAlgebra: dot, norm
using PoissonSolvers: FFTWBasis, Potential
using Test
using VlasovMethods

# the model on [0, 4π) × [-vlim, vlim] at f = (1 + α cos(x/2)) e^{-v²/2}/√(2π), with zero end
# columns
function model_of(T, nx, nv, vlim, α)
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, nx, nv), xd, (T(-vlim), T(vlim)))
    dist.values .= (1 .+ T(α) .* cos.(dist.x ./ 2)) .* exp.(-dist.v' .^ 2 ./ 2) ./
                   sqrt(T(2π))
    dist.values[:, [1, nv]] .= 0
    return VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
end

# the states of a run of `nsteps` steps of `τ`
function states(int, nsteps)
    nx, nv = size(int.model.distribution.values)
    sol = integrate(int.integrator)
    Fs = [reshape(Vector(sol.q[n]), nx, nv) for n in 0:ntime(int.equation)]
    @test length(Fs) == nsteps + 1
    return Fs
end

function bounded_run(model, τ, nsteps)
    T = typeof(τ)
    int = GeometricIntegrator(model, (zero(T), T(nsteps) * τ), τ;
        closure = BoundedVelocityClosure())
    return int.equation.parameters, states(int, nsteps)
end

function periodic_run(model, τ, nsteps)
    T = typeof(τ)
    int = GeometricIntegrator(model, (zero(T), T(nsteps) * τ), τ)
    return states(int, nsteps)
end

function cut(a::AbstractMatrix)
    b = copy(a)
    b[:, [1, end]] .= 0
    return b
end

# [a, b] on the periodic grid, from GeometricBrackets, which gives hx hv [c, û]
function bracket(p, a, b)
    A = p.arakawa
    return reshape(poisson_apply(A, vec(b), vec(a)) ./ (A.hx * A.hv), A.nx, A.nv)
end

density(p, F) = p.arakawa.hv .* vec(sum(cut(F); dims = 2))
cut_h(p, F) = cut(p.G * density(p, F) .+ p.kinetic')
function hamiltonian(p, F)
    A, ρ = p.arakawa, density(p, F)
    return A.hx * A.hv * sum(cut(F) .* p.kinetic') + (A.hx / 2) * dot(ρ, p.G * ρ)
end
mass(p, F) = p.arakawa.hx * p.arakawa.hv * sum(cut(F))
function boundary_flux(p, F)
    B = bracket(p, cut_h(p, F), cut(F))
    return -p.arakawa.hx * p.arakawa.hv * (sum(B[:, 1]) + sum(B[:, end]))
end

function drifts(p, Fs)
    [abs(hamiltonian(p, F) - hamiltonian(p, Fs[1])) / abs(hamiltonian(p, Fs[1]))
     for F in Fs]
end
# |M_{n+1} − M_n − τ flux(f_{n+½})| / |M_0| of each step
function balances(p, Fs, τ)
    M = [mass(p, F) for F in Fs]
    return [abs(M[n + 1] - M[n] - τ * boundary_flux(p, (Fs[n] .+ Fs[n + 1]) ./ 2)) /
            abs(M[1])
            for n in 1:(length(Fs) - 1)]
end

# The control: the test's own right-hand side with the uncut h = ϕ ⊗ 1 + 1 ⊗ v²/2 on every node
# and f cut to the unknowns, generic in the element type for the Newton Jacobian's Dual numbers.
function uncut_rhs!(ḟ, t, f, params)
    (; p) = params
    A = p.arakawa
    F = cut(reshape(collect(f), A.nx, A.nv))
    h = p.G * (A.hv .* vec(sum(F; dims = 2))) .+ p.kinetic'
    copyto!(ḟ, vec(cut(bracket(p, h, F))))
    return ḟ
end

function uncut_run(model, p, τ, nsteps)
    T = typeof(τ)
    q₀ = copy(vec(model.distribution.values))
    prob = GeometricEquations.ODEProblem(uncut_rhs!, (zero(T), T(nsteps) * τ), τ, q₀;
        parameters = (p = p,))
    sol = integrate(prob, ImplicitMidpoint())
    return [reshape(Vector(sol.q[n]), p.arakawa.nx, p.arakawa.nv) for n in 0:ntime(prob)]
end

const RUN_A = (nx = 16, nv = 32, vlim = 3.5, α = 0.5, τ = 0.05)

@testset "run A keeps the energy and balances the mass, $T" for T in (Float64, Float32)
    nsteps = 100
    (; nx, nv, vlim, α, τ) = RUN_A
    model = model_of(T, nx, nv, vlim, α)
    p, Fs = bounded_run(model, T(τ), nsteps)
    @test eltype(Fs[end]) === T

    @test maximum(drifts(p, Fs)) ≤ nsteps * 2eps(T)
    M = [mass(p, F) for F in Fs]
    @test abs(M[end] - M[1]) ≥ 1e-4 * abs(M[1])
    @test maximum(balances(p, Fs, T(τ))) ≤ 16eps(T)
    @test all(F -> all(==(0), F[:, [1, nv]]), Fs)

    # the control drifts the energy by at least ten times the bound
    Fu = uncut_run(model, p, T(τ), nsteps)
    @test maximum(drifts(p, Fu)) ≥ 10 * nsteps * 2eps(T)
end

@testset "run B agrees with the periodic run, $T" for T in (Float64, Float32)
    # the base's Float64 conservation row, whose f stays round-off at ±vmax
    nx, nv, τ, nsteps = 16, 64, T(0.05), 40
    model = model_of(T, nx, nv, 10, 0.01)
    p, Fb = bounded_run(model, τ, nsteps)
    Fp = periodic_run(model, τ, nsteps)
    @test eltype(Fb[end]) === eltype(Fp[end]) === T

    edge(F) = max(maximum(abs, F[:, 2]), maximum(abs, F[:, nv - 1])) / maximum(abs, F)
    @test maximum(edge, Fb) ≤ eps(T)
    @test maximum(norm(Fb[n] .- Fp[n]) / norm(Fp[n]) for n in eachindex(Fb)) ≤ 64eps(T)
    @test maximum(drifts(p, Fb)) ≤ nsteps * 2eps(T)
end

@testset "run A over 10⁴ steps, $T" for T in (Float64, Float32)
    nsteps = 10_000
    (; nx, nv, vlim, α, τ) = RUN_A
    model = model_of(T, nx, nv, vlim, α)
    p, Fs = bounded_run(model, T(τ), nsteps)
    @test maximum(drifts(p, Fs)) ≤ (T === Float64 ? 1e-12 : 1e-4)
    @test maximum(balances(p, Fs, T(τ))) ≤ 16eps(T)
end
