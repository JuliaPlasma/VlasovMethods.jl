# The independent side of the tests of the spline Galerkin model: the tables of a tensor
# spline space, built here from SimpleSplines' quadratures with dense mass matrices, and the
# quantities the tests compare with the model's: the plain Galerkin form, the potential, the
# energy, the mass, the boundary flux, the energy gradient and the plain field. Nothing here
# calls a function of the model.
using LinearAlgebra: Symmetric, cholesky, dot, norm
using SimpleSplines: SplineQuadrature, basis_values, quadrature_nodes, quadrature_weights,
                     PeriodicBSplineBasis, BSplineBasis, UniformMesh, Free, (..)
using PoissonSolvers: PoissonSolverSpline, solve

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

# the plain field Mx⁻¹ [T(Φᵢⱼ, h, f)] Mv⁻¹ with h = φ(f) + v²/2 at the nodes
function plain_field(t, A)
    φ̂ = spline_phi(t, A)
    hx = (t.X1' * φ̂) .+ zero(t.vq)'
    hv = zero(t.X0' * φ̂) .+ t.vq'
    _, fx, fv = grid_values(t, A)
    return minv(t, t.X0 * (t.W .* (hx .* fv .- hv .* fx)) * t.V0')
end

# the L² projection of f(x, v), by the tables
spline_project(t, f) = minv(t, t.X0 * (t.W .* f.(t.xq, t.vq')) * t.V0')

reflect(A) = A[:, end:-1:1]

landau_state(T, α = 0.2) = (x, v) -> (1 + T(α) * cos(x / 2)) * exp(-v^2 / 2) / sqrt(T(2π))
function pedestal_state(T)
    (x, v) -> landau_state(T)(x, v) +
              T(0.01) * (1 + T(0.5) * sin(x / 2)) * (1 + v / 15)
end
