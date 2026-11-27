@doc raw"""
    VlasovPoisson(dist::TensorSplineDistribution, potential::Potential)

The spline Galerkin model of 1D1V Vlasov–Poisson on the tensor spline space of `dist`,
``f = \sum_{ij} A_{ij} \Phi_{ij}`` with ``\Phi_{ij}(x, v) = B_i(x) C_j(v)``. Its field is the
Galerkin field of the cyclically averaged form [`_averaged_form`](@ref),

```math
\mathbb{M}_x \dot{A} \, \mathbb{M}_v = r , \qquad r_{ij} = \bar{T}(\Phi_{ij}, h, f) , \qquad
h = \phi(x) + \frac{v^2}{2} ,
```

with ``\phi = \sum_k \hat\phi_k B_k(x)``, ``\hat\phi`` the solve of `potential.solver` for the
load vector ``\mathbb{M}_x \hat\rho - \bar\rho \, \iota_x`` of the density
``\hat\rho = A \iota_v``, ``(\iota_v)_j = \int C_j \, dv``, ``(\iota_x)_k = \int B_k \, dx``:
the solution of ``-\phi'' = \rho - \bar\rho`` with ``\int \phi \, dx = 0``, for
``\rho = \int f \, dv`` and its mean ``\bar\rho = \int \rho \, dx / L`` on the x-domain of length
``L``, on a uniform or a nonuniform x-mesh. The averaged form is totally antisymmetric, so the
field keeps the quadratic energy [`_spline_energy`](@ref). The mass [`_spline_mass`](@ref) is not
kept: its rate is the boundary flux ``\frac{1}{3} \int \phi'(x) (f(x, V) - f(x, -V)) \, dx`` on
the velocity domain ``[-V, V]``, which vanishes where ``f`` does. [`GalerkinMidpoint`](@ref)
steps the model.

The model's `work` is the matrix of the weights of the tensor quadrature of the two bases. The
constructor throws an `ArgumentError` unless `potential.solver` is a `PoissonSolverSpline`, and
an `ArgumentError` unless the potential's basis is the x-basis of `dist`: a `PeriodicBSplineBasis`
of the same element type, number of functions, degree and breakpoints.
"""
VlasovPoisson(dist::TensorSplineDistribution, potential::Potential)

function _work(dist::TensorSplineDistribution, potential::Potential)
    solver = potential.solver
    solver isa PoissonSolverSpline || throw(ArgumentError(
        "the potential of the spline Vlasov–Poisson model solves -φ'' = ρ by a " *
        "PoissonSolverSpline on the x-basis of the distribution; got a " *
        "$(nameof(typeof(solver)))"))
    b, xb = basis(potential), dist.xbasis
    same = b isa PeriodicBSplineBasis && eltype(b) === eltype(xb) &&
           nbasis(b) == nbasis(xb) &&
           degree(b) == degree(xb) && breakpoints(b) == breakpoints(xb)
    same || throw(ArgumentError(
        "the potential of the spline Vlasov–Poisson model is on the x-basis of the " *
        "distribution, $(nbasis(xb)) functions of degree $(degree(xb)) in $(eltype(xb)) on " *
        "$(first(breakpoints(xb))) .. $(last(breakpoints(xb))); got $(nbasis(b)) functions of " *
        "degree $(degree(b)) in $(eltype(b)) on $(first(breakpoints(b))) .. " *
        "$(last(breakpoints(b)))"))
    return quadrature_weights(dist.xquadrature) * quadrature_weights(dist.vquadrature)'
end

function _check_coefficients(dist::TensorSplineDistribution, A::AbstractMatrix)
    size(A) == size(dist.coefficients) || throw(DimensionMismatch(
        "a coefficient matrix of the distribution has the size $(size(dist.coefficients)); " *
        "got $(size(A))"))
    return nothing
end

# the values of the spline of A and of its two first derivatives at the nodes of the tensor
# quadrature: (a, ∂ₓa, ∂ᵥa)
function _tensor_values(dist::TensorSplineDistribution, A::AbstractMatrix)
    X0, X1 = basis_values(dist.xquadrature, 0), basis_values(dist.xquadrature, 1)
    V0, V1 = basis_values(dist.vquadrature, 0), basis_values(dist.vquadrature, 1)
    XA = transpose(X0) * A
    return (XA * V0, (transpose(X1) * A) * V0, XA * V1)
end

@doc raw"""
    _galerkin_form(dist::TensorSplineDistribution, a, b, c)

The Galerkin trilinear form of the canonical bracket on the coefficient matrices `a`, `b` and
`c` of the bases of `dist`,

```math
T(a, b, c) = \int\!\!\int a \, (\partial_x b \, \partial_v c - \partial_v b \, \partial_x c)
\, dx \, dv ,
```

by the tensor rule of the two default `SplineQuadrature`s, whose ``\lceil 3p/2 \rceil`` Gauß
points per cell integrate the product exactly for the degree ``p``. On a clamped velocity basis
``T`` is antisymmetric in ``(b, c)``, and in ``(a, b)`` only up to the boundary terms at the ends
of the velocity domain. Throws a `DimensionMismatch` unless each matrix has the size of the
coefficients of `dist`.
"""
function _galerkin_form(
        dist::TensorSplineDistribution, a::AbstractMatrix, b::AbstractMatrix,
        c::AbstractMatrix)
    foreach(m -> _check_coefficients(dist, m), (a, b, c))
    W = quadrature_weights(dist.xquadrature) * quadrature_weights(dist.vquadrature)'
    a0, _, _ = _tensor_values(dist, a)
    _, bx, bv = _tensor_values(dist, b)
    _, cx, cv = _tensor_values(dist, c)
    return sum(W .* a0 .* (bx .* cv .- bv .* cx))
end

@doc raw"""
    _averaged_form(dist::TensorSplineDistribution, a, b, c)

The cyclic average of [`_galerkin_form`](@ref),

```math
\bar{T}(a, b, c) = \frac{1}{3} \left( T(a, b, c) + T(b, c, a) + T(c, a, b) \right) ,
```

which is totally antisymmetric on any bases, the Galerkin analogue of Arakawa's bracket. With the
energy gradient ``h``, ``\bar{T}(h, h, f) = 0``.
"""
function _averaged_form(
        dist::TensorSplineDistribution, a::AbstractMatrix, b::AbstractMatrix,
        c::AbstractMatrix)
    return (_galerkin_form(dist, a, b, c) + _galerkin_form(dist, b, c, a) +
            _galerkin_form(dist, c, a, b)) / 3
end

# the load vector b = Mx ρ̂ of the density ρ̂ = A ι_v, (ι_v)ⱼ = ∫ Cⱼ dv, and the coefficients φ̂
# of the potential of A, or of the load b: the solution of -φ'' = ρ - ρ̄ for the mean
# ρ̄ = ∫ρ dx / L, whose load b - ρ̄ ι_x, (ι_x)ₖ = ∫ Bₖ dx, sums to zero, taken with ∫φ dx = 0, so
# that φ̂ is the gradient of the field energy ½ φ̂ᵀ b in b. On a uniform mesh ι_x is constant, and
# both steps change the solver's mean-free solution by round-off only.
function _spline_load(model, A)
    dist = model.distribution
    mass_matrix(dist.xquadrature) * (A * basis_integrals(dist.vquadrature))
end
_spline_potential(model, A) = _load_potential(model, _spline_load(model, A))
function _load_potential(model, b)
    ι = basis_integrals(model.distribution.xquadrature)
    L = sum(ι)
    φ̂ = PoissonSolvers.solve(model.potential.solver, b .- (sum(b) / L) .* ι)
    return φ̂ .- dot(ι, φ̂) / L
end

# h = φ + v²/2, or φ alone, at the nodes of the tensor quadrature, with its two derivatives:
# h is a matrix, ∂ₓh a column and ∂ᵥh a row, which the products of the bracket broadcast
function _hamiltonian_values(dist, φ̂, kinetic::Bool)
    T = eltype(φ̂)
    X0, X1 = basis_values(dist.xquadrature, 0), basis_values(dist.xquadrature, 1)
    v = permutedims(quadrature_nodes(dist.vquadrature))
    k = kinetic ? one(T) : zero(T)
    return (transpose(X0) * φ̂ .+ k .* v .^ 2 ./ 2, transpose(X1) * φ̂, k .* v)
end

# r[i, j] = T̄(Φᵢⱼ, h, f) for the h of `_hamiltonian_values` and the f of A, the sum of
#   T(Φ, h, f) = ∫ Φ (hₓ f_v − h_v fₓ),  T(h, f, Φ) = ∫ h (fₓ Φ_v − f_v Φₓ),
#   T(f, Φ, h) = ∫ f (Φₓ h_v − Φ_v hₓ),
# collected by the factor of the test function: G00 multiplies Φ, G01 Φ_v and G10 Φₓ
function _averaged_bracket(dist, W, (h, hx, hv), A)
    f, fx, fv = _tensor_values(dist, A)
    X0, X1 = basis_values(dist.xquadrature, 0), basis_values(dist.xquadrature, 1)
    V0, V1 = basis_values(dist.vquadrature, 0), basis_values(dist.vquadrature, 1)
    G00 = W .* (hx .* fv .- hv .* fx)
    G01 = W .* (h .* fx .- f .* hx)
    G10 = W .* (f .* hv .- h .* fv)
    r = X0 * (G00 * transpose(V0) .+ G01 * transpose(V1)) .+ X1 * (G10 * transpose(V0))
    return r ./ 3
end

@doc raw"""
    _spline_1d1v_rhs!(Ȧ, model::VlasovPoisson{1, 1, <:TensorSplineDistribution}, A)

Write the field of the spline Galerkin model at the coefficients `A` into `Ȧ`, and return `Ȧ`:

```math
\dot{A} = \mathbb{M}_x^{-1} r \, \mathbb{M}_v^{-1} , \qquad
r_{ij} = \bar{T}(\Phi_{ij}, h, f) , \qquad h = \phi(f) + \frac{v^2}{2} ,
```

with ``f`` the spline of `A` and ``\bar{T}`` the form [`_averaged_form`](@ref). The field is in
the element type of `A`, the element type of the model. Throws a `DimensionMismatch` unless `Ȧ`
and `A` have the size of the coefficients of the model's distribution.
"""
function _spline_1d1v_rhs!(Ȧ::AbstractMatrix,
        model::VlasovPoisson{1, 1, <:TensorSplineDistribution{T}},
        A::AbstractMatrix{T}) where {T}
    dist = model.distribution
    _check_coefficients(dist, Ȧ)
    _check_coefficients(dist, A)
    hg = _hamiltonian_values(dist, _spline_potential(model, A), true)
    r = _averaged_bracket(dist, model.work, hg, A)
    Ȧ .= _tensor_mass_solve(dist.xquadrature, dist.vquadrature, r)
    return Ȧ
end

@doc raw"""
    _spline_energy(model::VlasovPoisson{1, 1, <:TensorSplineDistribution}, A)

The energy of the spline Galerkin model at the coefficients `A`,

```math
H = \sum_{ij} \Big( \int B_i \, dx \Big) A_{ij} \Big( \int C_j \frac{v^2}{2} \, dv \Big)
+ \frac{1}{2} \hat\phi^\mathsf{T} (\mathbb{M}_x \hat\rho) ,
```

the kinetic energy and the field energy ``\frac{1}{2} \int (\phi')^2 \, dx``. It is quadratic in
`A`, with the gradient ``\nabla H = \mathbb{M}_x \hat{h} \, \mathbb{M}_v`` for the coefficients
``\hat{h}`` of ``h = \phi + v^2/2``, so ``\nabla H \cdot \dot{A} = \bar{T}(h, h, f) = 0``.
"""
function _spline_energy(model::VlasovPoisson{1, 1, <:TensorSplineDistribution{T}},
        A::AbstractMatrix{T}) where {T}
    dist = model.distribution
    _check_coefficients(dist, A)
    qv = dist.vquadrature
    kinetic = basis_values(qv, 0) *
              (quadrature_weights(qv) .* quadrature_nodes(qv) .^ 2 ./ 2)
    load = _spline_load(model, A)
    φ̂ = _load_potential(model, load)
    return dot(basis_integrals(dist.xquadrature), A * kinetic) + dot(φ̂, load) / 2
end

@doc raw"""
    _spline_mass(model::VlasovPoisson{1, 1, <:TensorSplineDistribution}, A)

The mass ``\sum_{ij} (\int B_i \, dx) A_{ij} (\int C_j \, dv)`` of the coefficients `A`. Its rate
under the field is the boundary flux of [`VlasovPoisson`](@ref) on a tensor spline distribution.
"""
function _spline_mass(model::VlasovPoisson{1, 1, <:TensorSplineDistribution{T}},
        A::AbstractMatrix{T}) where {T}
    dist = model.distribution
    _check_coefficients(dist, A)
    return dot(basis_integrals(dist.xquadrature), A * basis_integrals(dist.vquadrature))
end

@doc raw"""
    GalerkinMidpoint(model::VlasovPoisson{1, 1, <:TensorSplineDistribution}, tspan::Tuple,
        tstep::Real; tol = 4eps(T), max_iterations = 10)

The implicit midpoint rule for the spline Galerkin model `model` over the time span `tspan` with
the step `tstep`, both converted to the element type `T` of the distribution. One step of ``\tau``
from ``A`` is ``A \leftarrow Y``, the solution of

```math
F(Y) = Y - A - \tau \, \mathbb{M}_x^{-1} r(\bar{f}) \, \mathbb{M}_v^{-1} = 0 ,
```

with ``\bar{f}`` the spline of ``(A + Y)/2`` and ``r`` the matrix of [`_spline_1d1v_rhs!`](@ref)
with ``h = \phi(\bar{f}) + v^2/2``. The field is the gradient of the quadratic energy
[`_spline_energy`](@ref) under a totally antisymmetric form, so the step keeps the energy to the
tolerance of the solve.

The solve is GeometricSolvers' Newton method on the vector `vec(Y)`, from `vec(A)`, with GMRES
of one unrestarted Krylov cycle of dimension 100 and the constant forcing ``10^{-3}``, and the
exact Jacobian-vector product

```math
W - \frac{\tau}{2} \mathbb{M}_x^{-1} \left( r_{\phi(W)}(\bar{f}) + r_{h(\bar{f})}(W) \right)
\mathbb{M}_v^{-1} ,
```

with ``r_h(f)`` the matrix ``r`` for the Hamiltonian ``h`` and the state ``f``; the first term
takes ``h = \phi(W)``, without the kinetic energy. It stops at
``\|F\|_2 \le \mathrm{tol} \, \|A\|_2``, with the Euclidean norm of the entries, after at most
`max_iterations` Newton iterations.

The method holds a copy of the coefficients of the distribution at construction, which every
[`run!`](@ref) starts from. The constructor throws an `ArgumentError` when a coefficient is not
finite, when `tstep` is not finite and positive in `T`, when `tol` is negative or `NaN`, and when
`max_iterations` is negative.
"""
struct GalerkinMidpoint{T, MT <: VlasovPoisson, AT <: AbstractMatrix{T}} <: VlasovMethod
    model::MT
    tspan::Tuple{T, T}
    tstep::T
    initial::AT
    tol::T
    max_iterations::Int
end

function GalerkinMidpoint(model::VlasovPoisson{1, 1, <:TensorSplineDistribution{T}},
        tspan::Tuple{Real, Real}, tstep::Real; tol::Real = 4eps(T),
        max_iterations::Integer = 10) where {T}
    A = model.distribution.coefficients
    all(isfinite, A) || throw(ArgumentError(
        "the coefficients of the distribution of GalerkinMidpoint are finite; got " *
        "$(A[findfirst(!isfinite, A)]) at $(Tuple(findfirst(!isfinite, A)))"))
    τ = T(tstep)
    isfinite(τ) && τ > 0 || throw(ArgumentError(
        "the step of GalerkinMidpoint is finite and positive in the element type $T of the " *
        "distribution; got tstep = $tstep"))
    tol ≥ 0 || throw(ArgumentError(
        "the tolerance of GalerkinMidpoint is not negative and not NaN; got tol = $tol"))
    max_iterations ≥ 0 || throw(ArgumentError(
        "the iteration count of GalerkinMidpoint is not negative; got " *
        "max_iterations = $max_iterations"))
    initial = copy(A)
    return GalerkinMidpoint{T, typeof(model), typeof(initial)}(
        model, (T(tspan[1]), T(tspan[2])), τ, initial, T(tol), Int(max_iterations))
end

# The residual F(y) = y − a − τ Mx⁻¹ r(f̄) Mv⁻¹ of the implicit midpoint step of τ from the
# coefficients A, a = vec(A), at the iterate y = vec(Y), with f̄ the spline of (A + Y)/2.
function _midpoint_residual!(r, model, A, τ, y)
    dist = model.distribution
    f̄ = (A .+ reshape(y, size(A))) ./ 2
    h̄ = _hamiltonian_values(dist, _spline_potential(model, f̄), true)
    ṙ = _averaged_bracket(dist, model.work, h̄, f̄)
    r .= y .- vec(A) .- τ .* vec(_tensor_mass_solve(dist.xquadrature, dist.vquadrature, ṙ))
    return r
end

# The derivative of `_midpoint_residual!` at y in the direction w. The field is bilinear in
# (h, f), and φ(W) enters h without the kinetic term, so the derivative is exact.
function _midpoint_jvp!(Jv, model, A, τ, y, w)
    dist = model.distribution
    f̄ = (A .+ reshape(y, size(A))) ./ 2
    W = reshape(w, size(A))
    hW = _hamiltonian_values(dist, _spline_potential(model, W), false)
    hf = _hamiltonian_values(dist, _spline_potential(model, f̄), true)
    ṙ = _averaged_bracket(dist, model.work, hW, f̄) .+
        _averaged_bracket(dist, model.work, hf, W)
    Jv .= w .- (τ / 2) .* vec(_tensor_mass_solve(dist.xquadrature, dist.vquadrature, ṙ))
    return Jv
end

# One implicit midpoint step of the method from A, the step `s` of a run; an error unless the
# solve ends with SUCCESS.
function _midpoint_step(method::GalerkinMidpoint, A::AbstractMatrix, s::Integer)
    model, τ = method.model, method.tstep
    a = vec(A)
    residual!(r, y, p) = _midpoint_residual!(r, model, A, τ, y)
    jvp!(Jv, y, w, p) = _midpoint_jvp!(Jv, model, A, τ, y, w)

    newton = GeometricSolvers.Newton(;
        linsolve = GeometricSolvers.GMRES(;
            maxiter = 100, forcing = GeometricSolvers.FixedForcing(1e-3)),
        jacobian = GeometricSolvers.JVPFunction(jvp!))
    target = method.tol * norm(A)
    options = GeometricSolvers.Options(eltype(A); f_abstol = target, f_reltol = 0,
        x_reltol = 0, max_iterations = method.max_iterations)
    y, status = GeometricSolvers.solve(
        GeometricSolvers.NonlinearProblem(residual!), copy(a), newton; options)
    status.code === GeometricSolvers.SUCCESS || error(
        "the implicit midpoint solve of GalerkinMidpoint at step $s returned $(status.code): " *
        "‖F‖ = $(status.fnorm) against the target $target after $(status.iterations) " *
        "iterations")
    return reshape(y, size(A))
end

"""
    run!(method::GalerkinMidpoint, h5file)

Advance `method` over its time span from the copy of the initial coefficients that it holds,
write every state into the HDF5 file `h5file`, write the last state into the coefficients of the
model's distribution, and return the distribution.

The time grid is `GeometricIntegrators.TimeSeries(t₀, t₁, τ)`, one time per state, the initial
one included, as for [`run!`](@ref)`(::GeometricIntegrator, h5file)`. The file holds the dataset
`z` of the size `(nx * nv, ntime)`, one column `vec(A)` per state, and the dataset `t` of the
`ntime` times, both in the element type of the distribution. The file is created, or
overwritten, and closed on return, also when a step throws. A solve that does not converge raises
an `ErrorException` that names its return code, `‖F‖` and the number of iterations, and leaves the
distribution unchanged.
"""
function run!(method::GalerkinMidpoint, h5file)
    dist = method.model.distribution
    T = eltype(dist)
    t₀, t₁ = method.tspan
    times = GeometricIntegrators.TimeSeries(t₀, t₁, method.tstep)
    nt = length(times)
    A = copy(method.initial)
    n = length(A)

    h5 = h5open(h5file, "w")
    try
        h5z = create_dataset(h5, "z", T, ((n, nt), (n, -1)), chunk = (n, 1))
        h5t = create_dataset(h5, "t", T, ((nt,), (-1,)), chunk = (1,))
        h5z[:, 1] = vec(A)
        h5t[1] = t₀
        for s in 1:(nt - 1)
            A = _midpoint_step(method, A, s)
            h5z[:, s + 1] = vec(A)
            h5t[s + 1] = times[s]
        end
    finally
        close(h5)
    end

    dist.coefficients .= A
    return dist
end
