struct MLBCache{T, PT <: ParticleDistribution, ST <: SplineDistribution{T}} <: Cache{T}
    pdist::PT
    sdist::ST

    v::Vector{T}

    J::Vector{T}
    dS::Vector{T}
    dS_dg::Vector{T}

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

        new{T, typeof(pdist), typeof(sdist)}(pdist, sdist, v, J, dS, dS_dg)
    end
end

function MLBCache(pdist::ParticleDistribution{T}, sdist::SplineDistribution{T}) where {T}
    MLBCache{T}(pdist, sdist)
end

@collision_cache MLBCache

struct MetriplecticLenardBernstein{
    D, XD, VD, DT <: DistributionFunction{<:Any, XD, VD}, ET <: Entropy, CT <: CacheDict} <:
       VlasovModel
    dist::DT
    entropy::ET
    ν::D

    cache::CT

    function MetriplecticLenardBernstein(
            dist::DistributionFunction{<:Any, XD, VD}, ent::Entropy; ν::D = 1.0) where {
            D, XD, VD}
        cache = CacheDict(MLBCache(dist, ent.dist))
        new{D, XD, VD, typeof(dist), typeof(ent), typeof(cache)}(dist, ent, ν, cache)
    end
end

@doc raw"""
    compute_J!(J, sdist, ::MetriplecticLenardBernstein)

The vector ``\mathbb{L}_k = \sum_i \mathbb{M}^{-1}_{ik} \int \varphi_i \, (1 + \log f_s) \, dv``
of the Landau manuscript's `eq:defn_Lk`, which is exactly the ``L^2`` projection of
``1 + \log f_s`` onto the basis.

This discretisation projects ``1 + \log f_s`` onto the spline space and differentiates the
*projection*, where `eq:velocity_ode` of the Lenard-Bernstein manuscript uses the pointwise
ratio ``f_s'(v_\alpha)/f_s(v_\alpha)``. The two agree only up to the projection error of the
logarithm. That is deliberate and is what makes this a genuine discrete-gradient system: `J`
is ``\partial S_h / \partial f_i`` for ``S_h = \int f_s \log f_s \, dv``.

!!! warning "The logarithm is unguarded in the manuscripts and guarded here"
    The earlier implementation wrote `0.5 * log(f_s^2)`. That is `log|f_s|` exactly, and its
    derivative is `f_s'/f_s` for either sign, so the drift it produces is finite and
    plausible-looking wherever `f_s < 0`. What it is not is the entropy: `S = ∫ f log f`
    requires `f > 0`, and where `f_s < 0` the H-theorem reverses — `dS/dt = -ν ∫ F²/f dv`
    becomes *positive*. So the construction converted a positivity violation from a visible
    `DomainError` into an invisible wrong answer, and removed the only diagnostic that would
    have caught it. Positivity is checked here instead.
"""
function compute_J!(J, sdist::SplineDistribution{T, XD, 1},
        ::MetriplecticLenardBernstein) where {T, XD}
    q = sdist.quadrature
    fs = sdist.spline
    x = quadrature_nodes(q)

    # Sampled on the quadrature grid the basis already carries, then contracted and solved
    # against the mass operator -- which is what `l2_projection!` does.
    g = similar(J, length(x))
    for r in eachindex(x)
        f = fs(x[r])
        f > 0 || throw(DomainError(f,
            "the projected distribution is non-positive at v = $(x[r]), so " *
            "log f_s and hence the discrete entropy are undefined there. Writing this as " *
            "0.5*log(f_s^2) would return log|f_s| and hide the violation; the H-theorem " *
            "reverses where f_s < 0."))
        g[r] = 1 + log(f)
    end

    l2_projection!(J, q, g)
    return J
end

@doc raw"""
    compute_dS!(dS, J, v, sdist, ::MetriplecticLenardBernstein, pdist)

``\partial S_h / \partial v_\alpha = w_\alpha \sum_k \mathbb{L}_k \, \varphi_k'(v_\alpha)``,
the Landau manuscript's `eq:entropy_derivative` up to its overall sign.

The weight ``w_\alpha`` is read from `pdist`, which is why the particle distribution is passed
in: it equals ``1/N`` for every initialiser in `examples/` but not for the importance-sampling
one, so `1/length(dS)` is not a substitute.
"""
function compute_dS!(dS, J, v::AbstractArray{ST}, sdist::SplineDistribution{ST, XD, 1},
        ::MetriplecticLenardBernstein, pdist::ParticleDistribution) where {ST, XD}
    b = sdist.basis
    N = nbasis(b)
    w = pdist.particles.w
    buf = zeros(ST, local_width(b))

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

``S_h = \int f_s \log f_s \, dv``, the quantity the manuscript's entropy figures report.

Integrated on the Gauß-Legendre grid of the basis rather than adaptively, so that it is the
same quadrature the discretisation itself uses. Returns the value alone; the earlier version
returned a `quadgk` error estimate alongside it and wrote `0.5*log(f^2)`, so it reported
``\int f \log|f|`` rather than the entropy.
"""
function compute_entropy(f, mlb::MetriplecticLenardBernstein)
    sdist = mlb.entropy.dist
    q = sdist.quadrature
    x = quadrature_nodes(q)
    w = quadrature_weights(q)

    S = zero(eltype(w))
    for r in eachindex(x)
        fr = f(x[r])
        fr > 0 || throw(DomainError(fr,
            "the distribution is non-positive at v = $(x[r]), so the entropy " *
            "∫ f log f dv is undefined there"))
        S += w[r] * fr * log(fr)
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

function compute_moments(v::AbstractArray{ST}, pdist::ParticleDistribution,
        ::MetriplecticLenardBernstein) where {ST}
    n = sum(pdist.particles.w)
    u = dot(pdist.particles.w, v) / n
    eps = dot(pdist.particles.w, v .^ 2) / n

    return n, u, eps
end

@doc raw"""
``(\varepsilon_h - u_h^2) \, \mathbb{L} \, \partial S_h / \partial v``, i.e. the bracket of
[`rhs_downstairs_factor!`](@ref) with its denominator cleared.

Useful where the common factor is supplied elsewhere; it is **not** ``\dot{v}``, which is what
`rhs_downstairs_factor!` returns.
"""
function rhs!(
        v̇::AbstractArray{ST}, v::AbstractArray{ST}, pdist::ParticleDistribution, n, u, eps,
        dS::AbstractArray{ST}, dS_sum, dS_v_sum, ::MetriplecticLenardBernstein) where {ST}
    w = view(pdist.particles.w, 1, :)
    v̇ .= .-n ./ w .* (eps - u^2) .* dS .+ (eps .- u .* v) .* dS_sum .+
         (v .- u) .* dS_v_sum
end

@doc raw"""
The metriplectic Lenard-Bernstein right-hand side, ``\dot{v} = \mathbb{L} \, \partial S_h / \partial v``
with the bracket

```math
\mathbb{L}_{\alpha\beta} = - \frac{n_h}{w_\alpha} \, \delta_{\alpha\beta}
  + \frac{(\varepsilon_h - u_h v_\alpha) + (v_\alpha - u_h) v_\beta}{\varepsilon_h - u_h^2} .
```

``\mathbb{L}`` is symmetric and satisfies ``\sum_\alpha w_\alpha \mathbb{L}_{\alpha\beta} = 0``
and ``\sum_\alpha w_\alpha v_\alpha \mathbb{L}_{\alpha\beta} = 0``, so momentum and energy are
exact Casimirs of the bracket **whatever** `dS` is — a stronger structure than the
manuscript's, and the reason these runs conserve.

The weight is per particle, read from `pdist` for each ``\alpha``: the two degeneracies above
are what fail if a single weight is used for every particle, so with non-uniform weights the
conservation the bracket is built for would be lost.
"""
function rhs_downstairs_factor!(
        v̇::AbstractArray{ST}, v::AbstractArray{ST}, pdist::ParticleDistribution, n, u, eps,
        dS::AbstractArray{ST}, dS_sum, dS_v_sum, ::MetriplecticLenardBernstein) where {ST}
    w = view(pdist.particles.w, 1, :)
    v̇ .= .-n ./ w .* dS .+ (eps .- u .* v) ./ (eps - u^2) .* dS_sum .+
         (v .- u) ./ (eps - u^2) .* dS_v_sum
end

function collisional_vectorfield!(v̇::AbstractArray{ST}, v::AbstractArray{ST}, params,
        mlb::MetriplecticLenardBernstein) where {ST}
    cache = mlb.cache[ST]

    sdist = cache.sdist

    projection(v, mlb.dist, sdist)

    compute_J!(cache.J, sdist, mlb)

    compute_dS!(cache.dS, cache.J, v, sdist, mlb, mlb.dist)

    ds_sum = sum(cache.dS)
    ds_v_sum = dot(v, cache.dS)

    n, u, eps = compute_moments(v, mlb.dist, mlb)

    rhs_downstairs_factor!(v̇, v, mlb.dist, n, u, eps, cache.dS, ds_sum, ds_v_sum, mlb)

    # Scale by the collision frequency, as every other collision model here does.
    v̇ .*= mlb.ν
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

function Picard_iterate_over_particles(dv::AbstractArray{ST}, vn::AbstractArray{ST},
        vn_minus_one::AbstractArray{ST}, dv_history, ti, t, Δt, m, β,
        abstol, reltol, mlb::MetriplecticLenardBernstein;
        maxiters::Int = 1000) where {ST}

    # set up vectors for storing intermediates
    # v_new = copy(vn)
    v_prev = copy(vn)

    dist = mlb.dist
    ent = mlb.entropy
    # cache = mlb.cache[ST]
    # sdist = cache.sdist

    # vnew_vec = zeros(1) # 1-vector for passing to SimpleSolvers as it does not support scalars
    # fvec = zero(vn)

    params = (dist = dist, ent = ent)

    # use Hermite extrapolation to get an initial guess
    if ti ≥ 4
        extrapolate!(t - 2Δt, vn_minus_one, dv_history[:, 2], t - Δt, vn,
            dv_history[:, 1], t, v_prev, HermiteExtrapolation())
    else
        problemGNI = GeometricEquations.ODEProblem(
            (v̇, t, v, params) -> collisional_vectorfield!(v̇, v, params, mlb),
            (t, t+Δt), Δt, vn; parameters = params)
        extrapolate!(
            t - Δt, vn, t, v_prev, problemGNI, MidpointExtrapolation(5))
    end

    # `SimpleSolvers` solves `F(x) = 0` by the fixed-point step `x ← x - α F(x)`. The implicit
    # midpoint residual here is `F(x) = x - vn - Δt·dv((x + vn)/2)`, and `f!` writes its negation,
    # so the sign flips once; the unaccelerated step ignores the `m` and `β` the signature carries.
    probN = SimpleSolvers.NonlinearProblem(
        (f, v, p) -> (f!(f, v, vn, params, Δt, mlb); f .*= -1), v_prev)

    # Build the solver and state directly: a non-finite step throws a `NonlinearSolverException`
    # before any status exists, and owning the state lets the error below report that step's
    # residual and iteration count.
    solver = SimpleSolvers.NonlinearSolver(SimpleSolvers.Picard(), v_prev, probN;
        f_abstol = abstol, f_reltol = reltol, max_iterations = maxiters, verbosity = 0)
    state = SimpleSolvers.SolverState(solver)

    status = try
        SimpleSolvers.solve_with_status!(v_prev, solver, state)
    catch e
        e isa SimpleSolvers.NonlinearSolverException || rethrow()
        SimpleSolvers.status(solver, state)
    end

    SimpleSolvers.isconverged(status) || error(
        "the Picard solve of the metriplectic Lenard–Bernstein step did not converge: " *
        "residual = $(status.rfₐ), iterations = $(status.iterations)")

    dv_history[:, 2] .= dv_history[:, 1]
    # dv_history[:, 1] .= dv
    collisional_vectorfield!(view(dv_history, :, 1), v_prev, params, mlb)

    return v_prev
end
