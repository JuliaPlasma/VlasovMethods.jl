# The independent side of the tests of the spline Galerkin model: the tables of a tensor
# spline space, built here from SimpleSplines' quadratures with dense mass matrices, and the
# quantities the tests compare with the model's: the plain Galerkin form, the potential, the
# energy, the mass, the boundary flux, the energy gradient and the plain field. Nothing here
# calls a function of the model except the long runs at the end, which four test files share:
# `test/verification/spline_galerkin_energy.jl` (Float64), `spline_galerkin_energy_float32.jl`,
# `spline_galerkin_landau_rates.jl` and `spline_galerkin_two_stream_rates.jl`. They are files of
# their own so that the CI job of their group can run them on separate workers.
using HDF5: h5read
using LinearAlgebra: Symmetric, cholesky, dot, norm
using SimpleSplines: SplineQuadrature, basis_values, quadrature_nodes, quadrature_weights,
                     PeriodicBSplineBasis, BSplineBasis, UniformMesh, Free, (..)
using PoissonSolvers: Potential, PoissonSolverSpline, solve
using Test
using VlasovMethods
using VlasovMethods: evaluate

# 16 periodic cubic x-splines of [0, L) and nv clamped cubic v-splines of [-V, V]
function spline_bases(T; nx = 16, nv = 100, L = 4π, V = 25)
    xb = PeriodicBSplineBasis(UniformMesh(nx, zero(T) .. T(L)), 3)
    vb = BSplineBasis(UniformMesh(nv - 3, T(-V) .. T(V)), 3, Free())
    return xb, vb
end

# the tables of the tensor rule with nq points per cell, or the default nq where nothing
function spline_tables(xb, vb; nq = nothing)
    qx = nq === nothing ? SplineQuadrature(xb) : SplineQuadrature(xb; nq)
    qv = nq === nothing ? SplineQuadrature(vb) : SplineQuadrature(vb; nq)
    X0, X1 = Matrix(basis_values(qx, 0)), Matrix(basis_values(qx, 1))
    V0, V1 = Matrix(basis_values(qv, 0)), Matrix(basis_values(qv, 1))
    xq, wx = quadrature_nodes(qx), quadrature_weights(qx)
    vq, wv = quadrature_nodes(qv), quadrature_weights(qv)
    Mx, Mv = X0 * (wx .* X0'), V0 * (wv .* V0')
    Mx, Mv = (Mx + Mx') / 2, (Mv + Mv') / 2
    return (; xb, vb, X0, X1, V0, V1, xq, vq, wx, wv, W = wx * wv', Mx, Mv,
        Cx = cholesky(Symmetric(Mx)), Cv = cholesky(Symmetric(Mv)),
        ix = X0 * wx, iv = V0 * wv, kv = V0 * (wv .* vq .^ 2 ./ 2),
        solver = PoissonSolverSpline(xb))
end

# values and first derivatives of a coefficient matrix at the nodes of the tensor rule
grid_values(t, A) = (t.X0' * A * t.V0, t.X1' * A * t.V0, t.X0' * A * t.V1)

# T(a, b, c) = ∫∫ a (∂ₓb ∂ᵥc − ∂ᵥb ∂ₓc) dx dv
function plain_form(t, a, b, c)
    a0, _, _ = grid_values(t, a)
    _, bx, bv = grid_values(t, b)
    _, cx, cv = grid_values(t, c)
    return sum(t.W .* a0 .* (bx .* cv .- bv .* cx))
end

minv(t, R) = (t.Cx \ R) / t.Cv

# the coefficients of the mean-free solution of -φ'' = ∫ f dv
spline_phi(t, A) = solve(t.solver, t.Mx * (A * t.iv))

# H = Σ (∫Bᵢ) A (∫Cⱼ v²/2) + ½ ∫ (φ')² dx, the field term by the x-quadrature of φ'
function spline_H(t, A)
    φx = t.X1' * spline_phi(t, A)
    return dot(t.ix, A * t.kv) + sum(t.wx .* φx .^ 2) / 2
end

spline_M(t, A) = dot(t.ix, A * t.iv)

# (1/3) ∫ φ'(x) (f(x, V) − f(x, −V)) dx by the x-quadrature
function spline_flux(t, A)
    φx = t.X1' * spline_phi(t, A)
    fb = t.X0' * (A[:, end] .- A[:, 1])
    return sum(t.wx .* φx .* fb) / 3
end

# the coefficients ĥ of h = φ + v²/2, from the projections of 1 and v²/2
function spline_hhat(t, A)
    one_x, one_v, κ = t.Cx \ t.ix, t.Cv \ t.iv, t.Cv \ t.kv
    return spline_phi(t, A) * one_v' .+ one_x * κ'
end

energy_gradient(t, A) = t.Mx * spline_hhat(t, A) * t.Mv

# Mx⁻¹ [T(Φᵢⱼ, h, A)] Mv⁻¹ for h = φ̂ + v²/2, or h = φ̂ alone where `kinetic` is false
function plain_bracket(t, φ̂, kinetic, A)
    hx = (t.X1' * φ̂) .+ zero(t.vq)'
    hv = zero(t.X0' * φ̂) .+ (kinetic ? t.vq' : zero(t.vq)')
    _, fx, fv = grid_values(t, A)
    return minv(t, t.X0 * (t.W .* (hx .* fv .- hv .* fx)) * t.V0')
end

# the plain field Mx⁻¹ [T(Φᵢⱼ, h, f)] Mv⁻¹ with h = φ(f) + v²/2 at the nodes
plain_field(t, A) = plain_bracket(t, spline_phi(t, A), true, A)

# the L² projection of f(x, v), by the tables
spline_project(t, f) = minv(t, t.X0 * (t.W .* f.(t.xq, t.vq')) * t.V0')

reflect(A) = A[:, end:-1:1]

landau_state(T, α = 0.2) = (x, v) -> (1 + T(α) * cos(x / 2)) * exp(-v^2 / 2) / sqrt(T(2π))
function pedestal_state(T)
    (x, v) -> landau_state(T)(x, v) +
              T(0.01) * (1 + T(0.5) * sin(x / 2)) * (1 + v / 15)
end

# The long runs. `z` and the times of `nsteps` steps of `τ` of `GalerkinMidpoint` on `model`
function galerkin_states(model, nsteps, τ; kwargs...)
    file = tempname() * ".hdf5"
    run!(GalerkinMidpoint(model, (0, nsteps * τ), τ; kwargs...), file)
    return h5read(file, "z"), h5read(file, "t")
end

# The energy and the mass over 10⁴ steps of the Landau state on 16 × 100 splines of
# [0, 4π) × [−25, 25], a domain on which f(±V) stays at round-off, with the default tolerance.
function energy_drift_test(T)
    xb, vb = spline_bases(T)
    t = spline_tables(xb, vb)
    model = VlasovPoisson(TensorSplineDistribution(xb, vb, landau_state(T)), Potential(xb))
    A = copy(model.distribution.coefficients)
    z, _ = galerkin_states(model, 10_000, T(0.05))
    @test size(z, 2) == 10_001 && eltype(z) === T
    H₀, M₀ = spline_H(t, A), spline_M(t, A)
    dH = maximum(abs(spline_H(t, reshape(c, 16, 100)) - H₀) for c in eachcol(z)) / abs(H₀)
    dM = maximum(abs(spline_M(t, reshape(c, 16, 100)) - M₀) for c in eachcol(z)) / abs(M₀)
    @test dH ≤ (T === Float64 ? 1e-12 : 1e-4)
    @test dM ≤ (T === Float64 ? 1e-13 : 1e-5)
end

# the least-squares slope of y over t
function slope(t, y)
    tm, ym = sum(t) / length(t), sum(y) / length(y)
    return sum((t .- tm) .* (y .- ym)) / sum(abs2, t .- tm)
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
# ω = ±1.415661889 − 0.1533594669 i and ω = 0.2592499317 i (RCM book, ch. 38). `rate_run`
# returns log E₁ and the times of `nsteps` steps from the projection of `f`.
function rate_run(T, f, nsteps)
    nx, nv, k, τ = 64, 128, 0.5, T(0.05)
    xb, vb = spline_bases(T; nx, nv, V = 10)
    t = spline_tables(xb, vb)
    model = VlasovPoisson(TensorSplineDistribution(xb, vb, f), Potential(xb))
    z, times = galerkin_states(model, nsteps, τ; tol = max(T(1e-12), 32eps(T)))
    return log_e1(t, z, nx, nv, k), Float64.(times)
end

function landau_rate_test(T)
    y, tt = rate_run(T, landau_state(T, 0.01), 800)
    peaks = [n for n in 2:(length(y) - 1) if tt[n] ≤ 40 && y[n - 1] < y[n] ≥ y[n + 1]]
    γ = slope(tt[peaks], y[peaks])
    @test γ < 0
    @test abs(γ / -0.15336 - 1) ≤ 0.02
end

function two_stream_rate_test(T)
    two_stream = (x, v) -> (1 + T(1e-3) * cos(x / 2)) * v^2 * exp(-v^2 / 2) / sqrt(T(2π))
    y, tt = rate_run(T, two_stream, 500)
    window = findall(s -> 5 ≤ s ≤ 25, tt)
    γ = slope(tt[window], y[window])
    @test γ > 0
    @test abs(γ / 0.25925 - 1) ≤ 0.05
end
