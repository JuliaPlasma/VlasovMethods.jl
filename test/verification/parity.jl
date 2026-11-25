# The long runs of `ParitySplitting` on the 16 × 32 grid of x ∈ [0, 4π), v ∈ [-7.5, 7.5], from
# f = (1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π): the energy and the mass over 10⁴ steps of τ = 0.05, the
# order against a Float64 RK4 reference of the full model, and an explicit control.
#
# `H = hx hv Σ f v²/2 + (hx/2) ρᵀ ϕ`, with `ϕ` the backend's field for `ρ = hv Σⱼ f`, is computed
# here and not in `src`. The f⁺ step solves its residual to `‖F‖₂ ≤ 4 eps(T) ‖f0 + f⁺‖₂`, so the
# energy drifts by the solver's round-off, not by a power of the step: the bounds are `1e-12` in
# Float64 and `1e-4` in Float32, against 2.9e-14 and 1.1e-6 measured on Julia 1.13.4 with
# GeometricSolvers 0.1.0. The mass bound is `100 eps(T)`, against 3.9 eps and 2.9 eps measured.

using GeometricBrackets: Arakawa, poisson_apply
# an explicit only-list, and not `using GeometricIntegrators`: both it and `VlasovMethods` export
# `GeometricIntegrator`, and the only-list binds it to this one
using GeometricIntegrators: RK4, integrate
using HDF5: h5read
using LinearAlgebra: dot, norm
using PoissonSolvers: FFTWBasis, Potential, solve
using Test
using VlasovMethods

function grid_model(T)
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, 16, 32), xd, (-T(7.5), T(7.5)))
    dist.values .= (1 .+ T(0.2) .* cos.(dist.x ./ 2)) .* exp.(-dist.v' .^ 2 ./ 2) ./
                   sqrt(T(2π))
    return VlasovPoisson(dist, Potential(FFTWBasis(xd, 16)))
end

function energy(model, F)
    dist = model.distribution
    hx, hv = step(dist.x), step(dist.v)
    ρ = hv .* vec(sum(F; dims = 2))
    return hx * hv * sum(F .* (dist.v' .^ 2 ./ 2)) +
           (hx / 2) * dot(ρ, solve(model.potential.solver, ρ))
end

mass(model, F) = step(model.distribution.x) * step(model.distribution.v) * sum(F)

# the states of a run of `ParitySplitting(grid_model(T), tspan, τ)`, as nx × nv matrices
function states(T, tspan, τ)
    model = grid_model(T)
    nx, nv = size(model.distribution.values)
    z = mktempdir() do dir
        file = joinpath(dir, "run.h5")
        run!(ParitySplitting(model, tspan, τ), file)
        h5read(file, "z")
    end
    return model, [reshape(z[:, n], nx, nv) for n in axes(z, 2)]
end

maxdrift(q) = maximum(abs, q .- q[1]) / abs(q[1])

# the run of 10⁴ steps of τ = 0.05 in `T`, made by the first testset that reads it, so that a run
# that throws fails that testset and not the file
const LONG = Dict{DataType, Any}()
long(T) = get!(() -> states(T, (0, 500), 0.05), LONG, T)

@testset "energy over 10⁴ steps of τ = 0.05, $T" for T in (Float32, Float64)
    model, Fs = long(T)
    @test length(Fs) == 10_001
    @test eltype(Fs[end]) == T
    @test maximum(abs, Fs[end] .- Fs[1]) > 1000 * eps(T) * maximum(abs, Fs[1])
    H = [energy(model, F) for F in Fs]
    @test eltype(H) == T
    @test maxdrift(H) ≤ (T == Float64 ? 1e-12 : 1e-4)
end

@testset "mass over 10⁴ steps of τ = 0.05, $T" for T in (Float32, Float64)
    model, Fs = long(T)
    m = [mass(model, F) for F in Fs]
    @test eltype(m) == T
    @test maxdrift(m) ≤ 100 * eps(T)
end

# The Float64 RK4 solution of the base's right-hand side `_grid_1d1v_rhs!` at τ = 10⁻³, at t = 2.
# Its own error, of order 10⁻¹², is far below the smallest error that it measures, of order 10⁻⁴.
const REFERENCE = let model = grid_model(Float64)
    int = GeometricIntegrator(model, (0.0, 2.0), 1e-3, RK4())
    reshape(collect(integrate(int.integrator).q[end]), 16, 32)
end

@testset "order: second order against an RK4 reference, $T" for T in (Float32, Float64)
    f₀ = grid_model(Float64).distribution.values
    τs = (0.1, 0.05, 0.025)
    errs = map(τs) do τ
        _, Fs = states(T, (0, 2), τ)
        @test eltype(Fs[end]) == T
        norm(Float64.(Fs[end]) .- REFERENCE) / norm(REFERENCE .- f₀)
    end
    lτ, le = log.(collect(τs)), log.(collect(errs))
    slope = sum((le .- sum(le) / 3) .* (lτ .- sum(lτ) / 3)) / sum(abs2, lτ .- sum(lτ) / 3)
    @test abs(slope - 2) ≤ 0.2
end

# The test's own parity splitting with the explicit f⁺ step f⁺ + τ P(f⁻, h(F0 + f⁺)), on the
# base's bracket term and Green's matrix: its energy drifts, so the energy testset can fail.
@testset "control: an explicit f⁺ step drifts, $T" for T in (Float32, Float64)
    model = grid_model(T)
    dist = model.distribution
    f = copy(dist.values)
    nx, nv = size(f)
    hx, hv = step(dist.x), step(dist.v)
    a = Arakawa(nx, nv, hx, hv)
    G = model.work
    P(g, h) = reshape(poisson_apply(a, vec(g), vec(h)), nx, nv) ./ (hx * hv)
    h(g) = G * (hv .* vec(sum(g; dims = 2))) .+ (dist.v' .^ 2 ./ 2)
    f0 = vec(sum((f .+ f[:, end:-1:1]) ./ 2; dims = 1)) ./ nx
    F0 = ones(T, nx) .* f0'
    fp, fm = parity_split(f, f0)
    τ = T(0.05)
    H0 = energy(model, f)
    drift = zero(T)
    for _ in 1:40
        fm = fm .+ τ / 2 .* P(F0 .+ fp, h(F0 .+ fp))
        fp = fp .+ τ .* P(fm, h(F0 .+ fp))
        fm = fm .+ τ / 2 .* P(F0 .+ fp, h(F0 .+ fp))
        drift = max(drift, abs(energy(model, parity_merge(f0, fp, fm)) - H0) / abs(H0))
    end
    @test eltype(fp) == eltype(fm) == T
    @test drift ≥ 1e-4
end
