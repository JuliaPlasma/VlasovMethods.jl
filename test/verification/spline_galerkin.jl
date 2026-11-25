using Test
using HDF5: h5read
using PoissonSolvers: Potential
using SimpleSplines: evaluate
using VlasovMethods

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

# the least-squares slope of y over t
function slope(t, y)
    tm, ym = sum(t) / length(t), sum(y) / length(y)
    return sum((t .- tm) .* (y .- ym)) / sum(abs2, t .- tm)
end

function run_states(model, nsteps, τ; kwargs...)
    file = tempname() * ".hdf5"
    run!(GalerkinMidpoint(model, (0, nsteps * τ), τ; kwargs...), file)
    return h5read(file, "z"), h5read(file, "t")
end

# The energy and the mass over 10⁴ steps of the Landau state on 16 × 100 splines of
# [0, 4π) × [−25, 25], a domain on which f(±V) stays at round-off, with the default tolerance.
@testset "energy over 10⁴ steps ($T)" for T in (Float64, Float32)
    xb, vb = spline_bases(T)
    t = spline_tables(xb, vb)
    model = VlasovPoisson(TensorSplineDistribution(xb, vb, landau_state(T)), Potential(xb))
    A = copy(model.distribution.coefficients)
    z, _ = run_states(model, 10_000, T(0.05))
    @test size(z, 2) == 10_001 && eltype(z) === T
    H₀, M₀ = spline_H(t, A), spline_M(t, A)
    dH = maximum(abs(spline_H(t, reshape(c, 16, 100)) - H₀) for c in eachcol(z)) / abs(H₀)
    dM = maximum(abs(spline_M(t, reshape(c, 16, 100)) - M₀) for c in eachcol(z)) / abs(M₀)
    @test dH ≤ (T === Float64 ? 1e-12 : 1e-4)
    @test dM ≤ (T === Float64 ? 1e-13 : 1e-5)
end

# E₁ = k (2/nx) |Σᵢ φ(xᵢ) exp(−i k xᵢ)| at the nx breakpoints xᵢ = (i − 1) 4π/nx
function log_e1(t, z, nx, nv, k)
    xs = [(i - 1) * 4π / nx for i in 1:nx]
    B = [evaluate(t.xb, j, eltype(z)(x)) for x in xs, j in 1:nx]
    return map(eachcol(z)) do c
        φ = B * spline_phi(t, reshape(c, nx, nv))
        log(k * (2 / nx) * abs(sum(φ[i] * exp(-im * k * xs[i]) for i in 1:nx)))
    end
end

# The rates of linear Landau damping and of the two-stream instability at k = 0.5 on 64 × 128
# splines of [0, 4π) × [−10, 10] at τ = 0.05, against the roots of the dispersion relation,
# ω = ±1.415661889 − 0.1533594669 i and ω = 0.2592499317 i (RCM book, ch. 38).
@testset "rates ($T)" for T in (Float64, Float32)
    nx, nv, k, τ = 64, 128, 0.5, T(0.05)
    xb, vb = spline_bases(T; nx, nv, V = 10)
    t = spline_tables(xb, vb)
    tol = max(T(1e-12), 32eps(T))

    landau = VlasovPoisson(TensorSplineDistribution(xb, vb, landau_state(T, 0.01)), Potential(xb))
    z, times = run_states(landau, 800, τ; tol)
    y = log_e1(t, z, nx, nv, k)
    tt = Float64.(times)
    peaks = [n for n in 2:(length(y) - 1) if tt[n] ≤ 40 && y[n - 1] < y[n] ≥ y[n + 1]]
    γ = slope(tt[peaks], y[peaks])
    @test γ < 0
    @test abs(γ / -0.15336 - 1) ≤ 0.02

    two_stream = (x, v) -> (1 + T(1e-3) * cos(x / 2)) * v^2 * exp(-v^2 / 2) / sqrt(T(2π))
    unstable = VlasovPoisson(TensorSplineDistribution(xb, vb, two_stream), Potential(xb))
    z, times = run_states(unstable, 500, τ; tol)
    y = log_e1(t, z, nx, nv, k)
    tt = Float64.(times)
    window = findall(s -> 5 ≤ s ≤ 25, tt)
    γ = slope(tt[window], y[window])
    @test γ > 0
    @test abs(γ / 0.25925 - 1) ≤ 0.05
end
