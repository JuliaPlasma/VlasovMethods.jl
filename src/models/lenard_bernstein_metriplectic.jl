# The residual of the Picard solve of one step, `x - vn - Δt·dv((x + vn)/2)`: `f!` writes its
# negation, so the sign flips once, and `SimpleSolvers`' fixed-point step is `x ← x - α F(x)`.
# It reads the step's data from `params`, so that one solver serves every step.
function _mlb_residual!(y, x, params)
    f!(y, x, params.vn, params, params.Δt, params.mlb)
    y .*= -1
end

# The unaccelerated Picard step reads no Jacobian. Passing this one keeps `SimpleSolvers` from
# building a `ForwardDiff` configuration, whose type depends on the length of the iterate.
function _mlb_jacobian!(j, x, params)
    throw(ArgumentError("the Picard solve of the metriplectic Lenard–Bernstein step " *
                        "evaluates no Jacobian"))
end

function _mlb_picard_solver(x::Vector{T}, abstol, reltol, maxiters) where {T}
    problem = SimpleSolvers.NonlinearProblem(_mlb_residual!, _mlb_jacobian!, x)
    SimpleSolvers.NonlinearSolver(SimpleSolvers.Picard(), x, problem;
        f_abstol = abstol, f_reltol = reltol, max_iterations = maxiters, verbosity = 0)
end

# A Picard solver of iterates of length `n`, its state, and the options it was built for.
mutable struct MLBPicard{T, ST, SST}
    solver::ST
    state::SST
    abstol::T
    reltol::T
    maxiters::Int
end

function MLBPicard{T}(n, abstol, reltol, maxiters) where {T}
    solver = _mlb_picard_solver(zeros(T, n), abstol, reltol, maxiters)
    state = SimpleSolvers.SolverState(solver)
    MLBPicard{T, typeof(solver), typeof(state)}(solver, state, abstol, reltol, maxiters)
end

# The solver of a step of `mlb` at the element type `T` with these options. The model keeps one
# solver, for the element type of its particles, and builds it again for other options. Another
# element type gets a new solver on each call.
function _mlb_picard!(mlb, ::Type{T}, abstol, reltol, maxiters) where {T}
    p = mlb.picard
    n = size(mlb.dist.particles.v, 2)
    p isa MLBPicard{T} || return MLBPicard{T}(n, abstol, reltol, maxiters)
    if !isequal(p.abstol, T(abstol)) || !isequal(p.reltol, T(reltol)) ||
       p.maxiters != maxiters
        p.solver = _mlb_picard_solver(zeros(T, n), abstol, reltol, maxiters)
        p.state = SimpleSolvers.SolverState(p.solver)
        p.abstol, p.reltol, p.maxiters = abstol, reltol, maxiters
    end
    return p
end

struct MLBCache{T, PT <: ParticleDistribution, ST <: SplineDistribution{T}} <: Cache{T}
    pdist::PT
    sdist::ST

    v::Vector{T}

    J::Vector{T}
    dS::Vector{T}
    dS_dg::Vector{T}

    # The work arrays of the vector field: the squared velocities of the moments, the load
    # vector of the projection, the samples of `1 + log|f_s|` at the quadrature nodes, and the
    # values or derivatives of the basis functions that are non-zero at one point.
    v²::Vector{T}
    rhs::Vector{T}
    g::Vector{T}
    values::Vector{T}

    function MLBCache{T}(pdist, sdist) where {T}
        M = length(sdist)
        N = size(pdist.particles.v, 2)

        # zeros(T, N), not zeros(N): `f!` writes the midpoint (vn .+ vp) ./ 2 into this
        # buffer, and a Float64 buffer truncates it under any wider element type -- a dual
        # number from a Jacobian, or extended precision.
        v = zeros(T, N)

        J = zeros(T, M)
        dS = zeros(T, N)
        dS_dg = zeros(T, N)

        v² = zeros(T, N)
        rhs = zeros(T, M)
        g = zeros(T, length(quadrature_nodes(sdist.quadrature)))
        values = _local_buffers(sdist.basis, T)

        new{T, typeof(pdist), typeof(sdist)}(
            pdist, sdist, v, J, dS, dS_dg, v², rhs, g, values)
    end
end

function MLBCache(pdist::ParticleDistribution{T}, sdist::SplineDistribution{T}) where {T}
    MLBCache{T}(pdist, sdist)
end

@collision_cache MLBCache

# The Picard solver of a new model, for the element type and the count of its particles. It is
# built at `abstol = reltol = 0` and `maxiters = 1`, and the first step builds it again for its
# own options.
function _mlb_picard(pdist::ParticleDistribution{T}) where {T}
    MLBPicard{T}(size(pdist.particles.v, 2), zero(T), zero(T), 1)
end

@doc raw"""
    LenardBernsteinBracket

The metric bracket of [`MetriplecticLenardBernstein`](@ref): [`DoughertyBracket`](@ref) or
[`JeyakumarBracket`](@ref). Both conserve momentum ``\sum_\gamma w_\gamma v_\gamma`` and energy
``\sum_\gamma w_\gamma v_\gamma^2`` for any entropy gradient, and their right-hand sides differ
by the factor ``n_h / (\varepsilon_h - u_h^2)``, which rescales friction against diffusion.
"""
abstract type LenardBernsteinBracket end

@doc raw"""
    DoughertyBracket()

The bracket `eq:lb_bracket_discrete` of the metriplectic Lenard-Bernstein manuscript, with the
kernel ``\varepsilon_h - u_h^2 + (v - u_h)(v' - u_h)``. Its equation of motion is
`eq:semidiscrete_equations_of_motion`,

```math
\dot{v}_\gamma = \nu \left[ - (\varepsilon_h - u_h^2) \, \frac{S_\gamma}{w_\gamma}
  + \frac{\varepsilon_h - u_h v_\gamma}{n_h} \sum_\beta S_\beta
  + \frac{v_\gamma - u_h}{n_h} \sum_\beta v_\beta S_\beta \right] ,
```

with ``S_\gamma = \partial S_h / \partial v_\gamma``, the moments ``n_h = \sum_\beta w_\beta``,
``u_h = \sum_\beta w_\beta v_\beta / n_h`` and ``\varepsilon_h = \sum_\beta w_\beta v_\beta^2 / n_h``,
and the weight ``w_\gamma`` of each particle. The default bracket of
[`MetriplecticLenardBernstein`](@ref).
"""
struct DoughertyBracket <: LenardBernsteinBracket end

@doc raw"""
    JeyakumarBracket()

The bracket of [`DoughertyBracket`](@ref) times ``n_h / (\varepsilon_h - u_h^2)``. Its
equation of motion is `eq:semidiscrete_equations_of_motion` of the metriplectic
Lenard-Bernstein manuscript times the same factor,

```math
\dot{v}_\gamma = \nu \left[ - n_h \, \frac{S_\gamma}{w_\gamma}
  + \frac{\varepsilon_h - u_h v_\gamma}{\varepsilon_h - u_h^2} \sum_\beta S_\beta
  + \frac{v_\gamma - u_h}{\varepsilon_h - u_h^2} \sum_\beta v_\beta S_\beta \right] ,
```

in the notation of [`DoughertyBracket`](@ref).
"""
struct JeyakumarBracket <: LenardBernsteinBracket end

@doc raw"""
    MetriplecticLenardBernstein(dist, ent; ν = 1.0, bracket = DoughertyBracket())

The metriplectic Lenard-Bernstein collision operator on the particles of `dist`, with the
entropy `ent`, the collision frequency `ν` and the metric bracket `bracket`, a
[`LenardBernsteinBracket`](@ref).

The right-hand side of either bracket throws a `DomainError` when the velocities have no spread
that the precision resolves: when
``\sum_\gamma w_\gamma (v_\gamma - u_h)^2 / n_h \le \epsilon \, \varepsilon_h`` with the machine
epsilon ``\epsilon`` of the element type, as for one particle or all particles at one velocity,
or when the computed ``\varepsilon_h - u_h^2`` is not positive. There
``\varepsilon_h - u_h^2`` is zero, or a rounding residue of either sign.
"""
struct MetriplecticLenardBernstein{
    D, XD, VD, DT <: DistributionFunction{<:Any, XD, VD}, ET <: Entropy, CT <: CacheDict,
    PST <: MLBPicard, BT <: LenardBernsteinBracket} <: VlasovModel
    dist::DT
    entropy::ET
    ν::D

    cache::CT
    picard::PST
    bracket::BT

    function MetriplecticLenardBernstein(
            dist::DistributionFunction{<:Any, XD, VD}, ent::Entropy; ν::D = 1.0,
            bracket::LenardBernsteinBracket = DoughertyBracket()) where {
            D, XD, VD}
        cache = CacheDict(MLBCache(dist, ent.dist))
        picard = _mlb_picard(dist)
        new{D, XD, VD, typeof(dist), typeof(ent), typeof(cache), typeof(picard),
            typeof(bracket)}(dist, ent, ν, cache, picard, bracket)
    end
end

# The value at `x` of the spline `c` on the basis `b`, with the buffer `values`: the sum of
# `SimpleSplines._evaluate_block!` in the same order, so the same value to the last bit, without
# the buffer that the point evaluation allocates. On a `Dirichlet` or `Neumann` basis the block
# of an end cell reaches past the basis, and those terms are skipped.
function _spline_value!(values, b, c, x)
    j₀ = evaluate_all!(values, b, x, 0)
    s = zero(eltype(values))
    for t in eachindex(values)
        j = basis_index(b, j₀ + t - 1)
        1 ≤ j ≤ nbasis(b) || continue
        s += values[t] * c[j]
    end
    return s
end

@doc raw"""
    compute_J!(J, g, values, sdist, ::MetriplecticLenardBernstein)

The vector
``\mathbb{L}_k = \sum_i \mathbb{M}^{-1}_{ik} \int \varphi_i \, (1 + \log |f_s|) \, dv`` of
`eq:defn_Lk` of the Landau and Lenard-Bernstein metriplectic manuscripts, which is exactly the
``L^2`` projection of ``1 + \log |f_s|`` onto the basis.

This discretisation projects ``1 + \log |f_s|`` onto the spline space and differentiates the
*projection*, where `eq:velocity_ode` of the Lenard-Bernstein manuscript uses the pointwise
ratio ``f_s'(v_\alpha)/f_s(v_\alpha)``. The two agree only up to the projection error of the
logarithm. That is deliberate and is what makes this a genuine discrete-gradient system: `J`
is ``\mathbb{M}^{-1}`` times ``\partial S_h / \partial f_i`` for the entropy
``S_h = \int f_s \log |f_s| \, dv`` of `eq:disc_entropy`.

``\log |f_s| = \frac{1}{2} \log (f_s^2)`` is defined where the ``L^2`` projection of a particle
distribution undershoots below zero. A quadrature node where ``f_s`` is zero or not finite
throws a `DomainError`.

`g`, of the length of the quadrature nodes, and `values`, of the local width of the basis, are
work arrays.
"""
function compute_J!(J, g, values, sdist::SplineDistribution{T, XD, 1},
        ::MetriplecticLenardBernstein) where {T, XD}
    q = sdist.quadrature
    fs = sdist.spline
    x = quadrature_nodes(q)

    # Sampled on the quadrature grid the basis already carries, then contracted and solved
    # against the mass operator -- which is what `l2_projection!` does.
    for r in eachindex(x)
        f = _spline_value!(values, fs.basis, fs.coefficients, x[r])
        _check_log_defined(f, x[r])
        g[r] = 1 + log(abs(f))
    end

    l2_projection!(J, q, g)
    return J
end

@doc raw"""
    compute_dS!(dS, buf, J, v, sdist, ::MetriplecticLenardBernstein, pdist)

``\partial S_h / \partial v_\alpha = w_\alpha \sum_k \mathbb{L}_k \, \varphi_k'(v_\alpha)``,
the Landau manuscript's `eq:entropy_derivative`.

The weight ``w_\alpha`` is read from `pdist`, which is why the particle distribution is passed
in: it equals ``1/N`` for every initialiser in `examples/` but not for the importance-sampling
one, so `1/length(dS)` is not a substitute. `buf`, of the local width of the basis, is a work
array.
"""
function compute_dS!(
        dS, buf, J, v::AbstractArray{ST}, sdist::SplineDistribution{ST, XD, 1},
        ::MetriplecticLenardBernstein, pdist::ParticleDistribution) where {ST, XD}
    b = sdist.basis
    N = nbasis(b)
    w = pdist.particles.w

    dS .= 0
    for i in eachindex(dS)
        j₀ = evaluate_all!(buf, b, v[i], 1)
        for t in eachindex(buf)
            j = basis_index(b, j₀ + t - 1)
            1 ≤ j ≤ N || continue
            dS[i] += w[1, i] * J[j] * buf[t]
        end
    end
    return dS
end

@doc raw"""
    compute_entropy(f, mlb::MetriplecticLenardBernstein)

``S_h = \int f \log |f| \, dv``, the entropy of `eq:disc_entropy` for ``f = f_s``: the
H-functional that decreases, and the quantity the manuscript's entropy figures report.

Integrated on the Gauß-Legendre grid of the basis rather than adaptively, so that it is the
same quadrature the discretisation itself uses. A node where `f` is zero or not finite throws a
`DomainError`, because ``\log |f|`` is undefined there.
"""
function compute_entropy(f, mlb::MetriplecticLenardBernstein)
    sdist = mlb.entropy.dist
    q = sdist.quadrature
    x = quadrature_nodes(q)
    w = quadrature_weights(q)

    S = zero(eltype(w))
    for r in eachindex(x)
        fr = f(x[r])
        _check_log_defined(fr, x[r])
        S += w[r] * fr * log(abs(fr))
    end
    return S
end

function compute_dS_discrete_gradient!(dS_dg, dS_midpoint, v_new::AbstractArray{ST},
        vn::AbstractArray{ST}, mlb::MetriplecticLenardBernstein) where {ST}
    sdist = mlb.cache[ST].sdist

    projection(vn, mlb.dist, sdist)
    S_n = compute_entropy(x -> sdist.spline(x), mlb)

    projection(v_new, mlb.dist, sdist)
    S_n_plus_1 = compute_entropy(x -> sdist.spline(x), mlb)

    dS_dg .= dS_midpoint .+
             (v_new .- vn) .* (S_n_plus_1 - S_n - dot(v_new .- vn, dS_midpoint)) ./
             norm(v_new .- vn) .^ 2
end

# `v²`, of the length of `v`, is a work array for the squared velocities.
function compute_moments(v::AbstractArray{ST}, v², pdist::ParticleDistribution,
        ::MetriplecticLenardBernstein) where {ST}
    n = sum(pdist.particles.w)
    u = dot(pdist.particles.w, v) / n
    v² .= v .^ 2
    eps = dot(pdist.particles.w, v²) / n

    return n, u, eps
end

# The right-hand side of the bracket, without the collision frequency: the equation of motion of
# the bracket's docstring divided by `ν`. `w` holds the weight of each particle, and `eps` is the
# energy moment.
function bracket_rhs!(v̇, v, w, n, u, eps, dS, dS_sum, dS_v_sum, ::DoughertyBracket)
    v̇ .= .-(eps - u^2) ./ w .* dS .+ (eps .- u .* v) ./ n .* dS_sum .+
         (v .- u) ./ n .* dS_v_sum
end

function bracket_rhs!(v̇, v, w, n, u, eps, dS, dS_sum, dS_v_sum, ::JeyakumarBracket)
    v̇ .= .-n ./ w .* dS .+ (eps .- u .* v) ./ (eps - u^2) .* dS_sum .+
         (v .- u) ./ (eps - u^2) .* dS_v_sum
end

# Throws where the velocities have no spread that the precision resolves: where the two-pass
# variance `σ²` is at most `eps(T) ε`, or where the computed `ε − u²` is not positive. For
# particles at one velocity, `ε − u²` is zero or a rounding residue of either sign, which no fixed
# multiple of `eps(T) ε` separates from a true spread, while their `σ²` is of order `eps(T)² v²`.
# Just above that threshold, the residue of `ε − u²` can still make it zero or negative. The
# right-hand sides divide by or multiply with `ε − u²`, which is positive wherever this passes, so
# the Jeyakumar form stays finite and the Dougherty diagonal keeps its sign. An overflow of `v²`
# that makes `σ²` and `ε` both infinite, or an underflow that makes both zero, throws too. A `NaN`
# fails both comparisons and propagates.
function _check_velocity_spread(v, w, n, u, ε)
    σ² = zero(u)
    for i in eachindex(v, w)
        σ² += w[i] * (v[i] - u)^2
    end
    σ² /= n
    (σ² ≤ eps(typeof(σ²)) * ε || ε - u^2 ≤ 0) && throw(DomainError(ε - u^2,
        "the velocities have no spread that the precision resolves, so ε − u² of the " *
        "metriplectic Lenard–Bernstein bracket is not positive"))
    return nothing
end

# The velocity field `ν 𝕃 dS` of the model's bracket at the velocities `v`, for the entropy
# gradient `dS`, written into `v̇`.
function metriplectic_rhs!(v̇::AbstractArray{ST}, v::AbstractArray{ST}, dS,
        mlb::MetriplecticLenardBernstein) where {ST}
    cache = mlb.cache[ST]
    w = view(mlb.dist.particles.w, 1, :)

    ds_sum = sum(dS)
    ds_v_sum = dot(v, dS)

    n, u, eps = compute_moments(v, cache.v², mlb.dist, mlb)
    _check_velocity_spread(v, w, n, u, eps)

    bracket_rhs!(v̇, v, w, n, u, eps, dS, ds_sum, ds_v_sum, mlb.bracket)

    # Scale by the collision frequency, as every other collision model here does.
    v̇ .*= mlb.ν
end

function collisional_vectorfield!(v̇::AbstractArray{ST}, v::AbstractArray{ST}, params,
        mlb::MetriplecticLenardBernstein) where {ST}
    cache = mlb.cache[ST]

    sdist = cache.sdist

    _projection!(cache.rhs, cache.values, v, mlb.dist, sdist)

    compute_J!(cache.J, cache.g, cache.values, sdist, mlb)

    compute_dS!(cache.dS, cache.values, cache.J, v, sdist, mlb, mlb.dist)

    metriplectic_rhs!(v̇, v, cache.dS, mlb)
end

function f!(f::AbstractArray{T}, vn::AbstractArray{T}, vp, params,
        Δt, mlb::MetriplecticLenardBernstein) where {T}
    v_midpoint = mlb.cache[T].v
    v_midpoint .= (vn .+ vp) ./ 2

    collisional_vectorfield!(f, v_midpoint, params, mlb)

    f .*= Δt
    f .-= (vn .- vp)

    return nothing
end

# The guess of the first three steps: a midpoint extrapolation of the vector field from
# `(t - Δt, vn)` to `t`, written into `v`. It allocates. The `EquationProblem` constructor of
# `GeometricEquations` sets the first type parameter of the problem with `eval`, so the
# assertion `::GeometricEquations.ODEProblem` states it, and the `extrapolate!` call dispatches
# statically.
function _mlb_startup_guess!(v, vn, t, Δt, params, mlb)
    problemGNI = GeometricEquations.ODEProblem(
        (v̇, t, v, params) -> collisional_vectorfield!(v̇, v, params, mlb),
        (t, t+Δt), Δt, vn; parameters = params)::GeometricEquations.ODEProblem
    extrapolate!(t - Δt, vn, t, v, problemGNI, MidpointExtrapolation(5))
end

"""
    Picard_iterate_over_particles!(v, vn, vn_minus_one, dv_history, ti, t, Δt, abstol, reltol,
        mlb::MetriplecticLenardBernstein; maxiters = 1000)

Solve one implicit-midpoint step of the metriplectic Lenard–Bernstein model from the velocities
`vn` by a Picard iteration, write the solved velocities into `v` and return `v`.

`ti` counts the steps from 1, and `t` is the time of the guess. From the fourth step on, the
guess is the Hermite extrapolation from `vn_minus_one` and `vn` with the derivatives in the
columns 2 and 1 of `dv_history`. The first three steps (`ti < 4`) take the guess from a
midpoint extrapolation of the vector field instead. After the solve, column 2 of `dv_history`
takes column 1, and column 1 the vector field at `v`. `abstol`, `reltol` and `maxiters` are the
`f_abstol`, `f_reltol` and `max_iterations` of the `SimpleSolvers` solve.

The model `mlb` keeps the solver and its state for the element type of its particles, and the
cache of `mlb` for the element type of `vn` keeps the work arrays of the vector field. The first
call builds the solver, and a call with other options builds it again. A step with `ti ≥ 4`
allocates nothing when `vn` has the element type of the particles of `mlb`; the guess of the
first three steps allocates, and a `vn` of another element type builds a solver on each call.

A solve that does not converge throws an `ErrorException` with the residual and the iteration
count. An exception of the model, such as the `DomainError` of a particle that leaves the
velocity domain, reaches the caller unchanged.

Throws a `DimensionMismatch` when `v` and `vn` differ in length, when `vn` or `vn_minus_one` has
another length than the particle count of `mlb`, or when `dv_history` does not have one row per
particle and two columns. Throws an `ArgumentError` when `v` shares memory with `vn`,
`vn_minus_one` or `dv_history`, or when `vn` or `vn_minus_one` shares memory with `dv_history`.
Each check runs before the step writes to any argument.
"""
function Picard_iterate_over_particles!(v::AbstractArray{ST}, vn::AbstractArray{ST},
        vn_minus_one::AbstractArray{ST}, dv_history, ti, t, Δt, abstol, reltol,
        mlb::MetriplecticLenardBernstein; maxiters::Int = 1000) where {ST}
    length(v) == length(vn) || throw(DimensionMismatch(
        "the output `v` has $(length(v)) entries and `vn` has $(length(vn))"))
    N = size(mlb.dist.particles.v, 2)
    length(vn) == N || throw(DimensionMismatch(
        "`vn` has $(length(vn)) entries and the model has $(N) particles"))
    length(vn_minus_one) == N || throw(DimensionMismatch(
        "`vn_minus_one` has $(length(vn_minus_one)) entries and the model has $(N) particles"))
    size(dv_history) == (N, 2) || throw(DimensionMismatch(
        "`dv_history` has size $(size(dv_history)) and the model has $(N) particles"))
    Base.mightalias(v, vn) && throw(ArgumentError("the output `v` shares memory with `vn`"))
    Base.mightalias(v, vn_minus_one) &&
        throw(ArgumentError("the output `v` shares memory with `vn_minus_one`"))
    Base.mightalias(v, dv_history) &&
        throw(ArgumentError("the output `v` shares memory with `dv_history`"))
    Base.mightalias(vn, dv_history) &&
        throw(ArgumentError("`vn` shares memory with `dv_history`"))
    Base.mightalias(vn_minus_one, dv_history) &&
        throw(ArgumentError("`vn_minus_one` shares memory with `dv_history`"))

    picard = _mlb_picard!(mlb, ST, abstol, reltol, maxiters)
    params = (dist = mlb.dist, ent = mlb.entropy, vn = vn, Δt = Δt, mlb = mlb)

    # use Hermite extrapolation to get an initial guess
    if ti ≥ 4
        extrapolate!(t - 2Δt, vn_minus_one, view(dv_history, :, 2), t - Δt, vn,
            view(dv_history, :, 1), t, v, HermiteExtrapolation())
    else
        _mlb_startup_guess!(v, vn, t, Δt, params, mlb)
    end

    # The state is passed in, not built by the solve: a non-finite step throws a
    # `NonlinearSolverException` before any status exists, and the state lets the error below
    # report that step's residual and iteration count.
    status = try
        SimpleSolvers.solve_with_status!(v, picard.solver, picard.state, params)
    catch e
        e isa SimpleSolvers.NonlinearSolverException || rethrow()
        SimpleSolvers.status(picard.solver, picard.state)
    end

    SimpleSolvers.isconverged(status) || error(
        "the Picard solve of the metriplectic Lenard–Bernstein step did not converge: " *
        "residual = $(status.rfₐ), iterations = $(status.iterations)")

    @views dv_history[:, 2] .= dv_history[:, 1]
    collisional_vectorfield!(view(dv_history, :, 1), v, params, mlb)

    return v
end

"""
    Picard_iterate_over_particles(dv, vn, vn_minus_one, dv_history, ti, t, Δt, m, β, abstol,
        reltol, mlb::MetriplecticLenardBernstein; maxiters = 1000)

The step of [`Picard_iterate_over_particles!`](@ref) into a new vector, which it returns. It
ignores `dv`, `m` and `β`.
"""
function Picard_iterate_over_particles(dv::AbstractArray{ST}, vn::AbstractArray{ST},
        vn_minus_one::AbstractArray{ST}, dv_history, ti, t, Δt, m, β,
        abstol, reltol, mlb::MetriplecticLenardBernstein;
        maxiters::Int = 1000) where {ST}
    Picard_iterate_over_particles!(similar(vn), vn, vn_minus_one, dv_history, ti, t, Δt,
        abstol, reltol, mlb; maxiters = maxiters)
end
