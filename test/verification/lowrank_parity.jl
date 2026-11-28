# The long runs of `LowRankParitySplitting` on the grid of x ∈ [0, 4π), v ∈ [-7.5, 7.5], from
# f = (1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π) with the default background: the energy and the mass over
# 10⁴ steps of τ = 0.05 on 16 × 32, an explicit control, the order against a Float64 RK4 reference
# of the full model, and the error against the full-rank `ParitySplitting` at four ranks on
# 64 × 128.
#
# `H = hx hv Σ f v²/2 + (hx/2) ρᵀ G ρ` is computed here from `GridDistribution(Y)`. The substeps
# are implicit midpoint steps solved exactly by two fixed-point iterations, so the energy drifts by
# round-off, not by a power of the step: the bounds are `1e-12` in Float64 and
# `8 · 3 · 10⁴ eps(Float32) = 2.9e-2` in Float32, three substeps per step, against at most 3.0e-14
# and 3.2e-4 measured at the ranks (6, 6) and (4, 4) on Julia 1.13.4. The mass bound is
# `1000 eps(T)`, against at most 91 eps measured.

using GeometricIntegrators: RK4, integrate
using HDF5: h5read
using LinearAlgebra: dot, norm
using PoissonSolvers: FFTWBasis, Potential
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V

include("../helpers/dense_parity.jl")

maxwellian(v) = exp(-v^2 / 2) / sqrt(2 * oftype(v, π))

function grid_model(T; nx = 16, nv = 32)
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, nx, nv), xd, (T(-7.5), T(7.5)))
    dist.values .= [(1 + T(0.2) * cos(x / 2)) * maxwellian(v) for x in dist.x, v in dist.v]
    return VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
end

function energy(dist, G, f)
    hx, hv = step(dist.x), step(dist.v)
    ρ = hv .* vec(sum(f; dims = 2))
    return hx * hv * sum(f .* (collect(dist.v)' .^ 2 ./ 2)) + (hx / 2) * dot(ρ, G * ρ)
end

mass(dist, f) = step(dist.x) * step(dist.v) * sum(f)
maxdrift(q) = maximum(abs, q .- q[1]) / abs(q[1])

# `nsteps` steps of τ at `rank` from the model's state; with `track`, the energy and the mass of
# every state
function lowrank_run(T, rank, nsteps, τ; nx = 16, nv = 32, track = true)
    model = grid_model(T; nx, nv)
    dist = model.distribution
    Y = LowRankParityDistribution(dist, rank)
    field = LowRankFields1D1V(Y.plus.space, model.work)
    method = LowRankParitySplitting(rank)
    H, m = T[], T[]
    record!(f) = track && (push!(H, energy(dist, model.work, f)); push!(m, mass(dist, f)))
    record!(GridDistribution(Y).values)
    for n in 1:nsteps
        lowrank_step!(Y, method, field, τ * (n - 1), τ)
        record!(GridDistribution(Y).values)
    end
    return (; model, Y, H, m, f = GridDistribution(Y).values)
end

# the runs of 10⁴ steps, made by the first testset that reads them
const LONG = Dict{Any, Any}()
long(T, r) = get!(() -> lowrank_run(T, (r, r), 10_000, 0.05), LONG, (T, r))

@testset "energy over 10⁴ steps of τ = 0.05, $T, rank ($r, $r)" for T in (Float32, Float64),
    r in (6, 4)

    run = long(T, r)
    @test length(run.H) == 10_001
    @test eltype(run.H) == T && eltype(run.Y.plus.S) == T
    @test maximum(abs, run.f .- run.model.distribution.values) >
          1000 * eps(T) * maximum(abs, run.model.distribution.values)
    @test maxdrift(run.H) ≤ (T == Float64 ? 1e-12 : 8 * 3 * 10^4 * eps(Float32))
end

@testset "mass over 10⁴ steps of τ = 0.05, $T, rank ($r, $r)" for T in (Float32, Float64),
    r in (6, 4)

    run = long(T, r)
    @test eltype(run.m) == T
    @test maxdrift(run.m) ≤ 1000 * eps(T)
end

@testset "control: explicit substeps drift, $T" for T in (Float32, Float64)
    model = grid_model(T)
    dist = model.distribution
    Y = LowRankParityDistribution(dist, (6, 6))
    d = dense_parity_case(Y.plus.space, model.work, Y.f0)
    op, om = dense_objects(Y)
    H₀ = energy(dist, model.work, GridDistribution(Y).values)
    drift = zero(T)
    for _ in 1:400
        op, om = dense_parity_step(d, op, om, T(0.05);
            tol = (zero(T), zero(T)), rmax = (6, 6), iters = 1)
        f = parity_merge(Y.f0, dense_object(op), dense_object(om))
        drift = max(drift, abs(energy(dist, model.work, f) - H₀) / abs(H₀))
    end
    @test eltype(op.S) == T
    # above the pass bound of the energy over 400 steps in test/lowrank/parity.jl; the explicit
    # substeps drift by 5.6e-3 in both precisions on Julia 1.13.4
    @test drift > (T == Float64 ? 1e-12 : 8 * 3 * 400 * eps(Float32))
end

# The Float64 RK4 solution of the model's right-hand side at τ = 10⁻³, at t = 2. Its own error, of
# order 10⁻¹², is far below the smallest error that it measures, of order 10⁻⁴.
const REFERENCE = let model = grid_model(Float64)
    int = GeometricIntegrator(model, (0.0, 2.0), 1e-3, RK4())
    reshape(collect(integrate(int.integrator).q[end]), 16, 32)
end

@testset "order: second order at rank (8, 8) against an RK4 reference, $T" for T in (Float32, Float64)
    f₀ = grid_model(Float64).distribution.values
    τs = (0.1, 0.05, 0.025)
    errs = map(τs) do τ
        run = lowrank_run(T, (8, 8), round(Int, 2 / τ), τ; track = false)
        @test eltype(run.f) == T
        norm(Float64.(run.f) .- REFERENCE) / norm(REFERENCE .- f₀)
    end
    lτ, le = log.(collect(τs)), log.(collect(errs))
    slope = sum((le .- sum(le) / 3) .* (lτ .- sum(lτ) / 3)) / sum(abs2, lτ .- sum(lτ) / 3)
    @test abs(slope - 2) ≤ 0.2
end

@testset "the error against the full rank at r ∈ {4, 8, 16, 32}, $T" for T in (Float32, Float64)
    model = grid_model(T; nx = 64, nv = 128)
    f₀ = copy(model.distribution.values)
    z = mktempdir() do dir
        file = joinpath(dir, "run.h5")
        run!(ParitySplitting(model, (0, 5), 0.05), file)
        h5read(file, "z")
    end
    f = reshape(z[:, end], 64, 128)
    @test size(z, 2) == 101
    errs = map((4, 8, 16, 32)) do r
        run = lowrank_run(T, (r, r), 100, 0.05; nx = 64, nv = 128, track = false)
        @test eltype(run.f) == T
        norm(run.f - f) / norm(f - f₀)
    end
    bounds = T == Float64 ? (0.5, 1e-3, 1e-10, 1e-11) : (0.5, 2e-3, 1e-3, 1e-3)
    @test all(errs .≤ bounds)
    @test errs[1] > errs[2] > errs[3]
end
