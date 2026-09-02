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

function Cache(AT, c::MLBCache{DT, PT, ST}) where {DT, PT, ST}
    MLBCache{AT}(c.pdist, similar(AT, c.sdist))
end
function CacheType(AT, c::MLBCache{DT, PT, ST}) where {DT, PT, ST}
    MLBCache{AT, PT, similar_type(AT, c.sdist)}
end

struct MetriplecticLenardBernstein{
    D, XD, VD, DT <: DistributionFunction{XD, VD}, ET <: Entropy, CT <: CacheDict} <:
       VlasovModel
    dist::DT
    entropy::ET
    ν::D

    cache::CT

    function MetriplecticLenardBernstein(
            dist::DistributionFunction{XD, VD}, ent::Entropy; ν::D = 1.0) where {D, XD, VD}
        cache = CacheDict(MLBCache(dist, ent.dist))
        new{D, XD, VD, typeof(dist), typeof(ent), typeof(cache)}(dist, ent, ν, cache)
    end
end

# `mlb.dist` and `mlb.entropy.dist`: the fields are named `dist` and `entropy`, so the earlier
# `mlb.pdist` / `mlb.sdist` would have thrown. Unreachable in practice — `CacheDict`'s parent
# is an `MLBCache`, so the `::MLBCache` methods above are the ones that fire — but wrong as
# written.
function Cache(AT, mlb::MetriplecticLenardBernstein)
    MLBCache{AT}(mlb.dist, similar(AT, mlb.entropy.dist))
end
function CacheType(AT, mlb::MetriplecticLenardBernstein)
    MLBCache{AT, typeof(mlb.dist), similar_type(AT, mlb.entropy.dist)}
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
    # against the mass operator -- which is what `l2_projection!` does. The earlier version
    # called `BSplineKit.galerkin_projection` and then `ldiv!` with the mass factorisation,
    # i.e. the same two steps through a different package.
    g = similar(J, length(x))
    for r in eachindex(x)
        f = fs(x[r])
        f > 0 || throw(ErrorException(
            "the projected distribution is non-positive, f_s = $(f) at v = $(x[r]), so " *
            "log f_s and hence the discrete entropy are undefined there. Writing this as " *
            "0.5*log(f_s^2) would return log|f_s| and hide the violation; the H-theorem " *
            "reverses where f_s < 0."))
        g[r] = 1 + log(f)
    end

    l2_projection!(J, q, g)
    return J
end

@doc raw"""
    compute_dS!(dS, J, v, sdist, ::MetriplecticLenardBernstein)

``\partial S_h / \partial v_\alpha = w_\alpha \sum_k \mathbb{L}_k \, \varphi_k'(v_\alpha)``,
the Landau manuscript's `eq:entropy_derivative` up to its overall sign.

The particle weight is ``w_\alpha``. The earlier version divided by `length(dS)` instead,
which equals ``w_\alpha`` only when every weight is ``1/N`` — true of every initialiser in
`examples/`, and false for the importance-sampling one.
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
        fr > 0 || throw(ErrorException(
            "the distribution is non-positive, f = $(fr) at v = $(x[r]), so the entropy " *
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

The weight is per particle. Both right-hand sides used `pdist.particles.w[1]`, particle one's
weight, for every particle; the commented-out single-particle version above them had it right.
The two degeneracies above are what fail when it is wrong, so with non-uniform weights the
conservation the bracket is built for is lost.
"""
function rhs!(
        v̇::AbstractArray{ST}, v::AbstractArray{ST}, pdist::ParticleDistribution, n, u, eps,
        dS::AbstractArray{ST}, dS_sum, dS_v_sum, ::MetriplecticLenardBernstein) where {ST}
    w = view(pdist.particles.w, 1, :)
    v̇ .= .-n ./ w .* (eps - u^2) .* dS .+ (eps .- u .* v) .* dS_sum .+
         (v .- u) .* dS_v_sum
end

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

    # The collision frequency was declared, stored, and then never read: `ν` appeared nowhere
    # after the constructor, so `MetriplecticLenardBernstein(dist, ent; ν = 0.1)` silently ran
    # at ν = 1. Every other collision model here multiplies by it.
    rhs_downstairs_factor!(v̇, v, mlb.dist, n, u, eps, cache.dS, ds_sum, ds_v_sum, mlb)
    v̇ .*= mlb.ν
end

# # single particle rhs for Newton solve
# function rhs!(v̇, v, ind::Int, v_all, pdist::ParticleDistribution, n, u, eps, dS::AbstractArray{ST}, dS_sum, dS_v_sum, ::MetriplecticLenardBernstein) where {ST}
#     # dS_sum = sum(dS)
#     v̇[1] = - n / pdist.particles.w[ind] * (eps - u^2) * dS[ind] + (eps - u * v[1]) * dS_sum + (v[1] - u) * dS_v_sum
#     # v̇[1] = - n / pdist.particles.w[ind] * (eps - u^2) * dS[ind] + (eps - u * v[1]) * dS_sum + (v[1] - u) * dot(v_all, dS)
# end

# function rhs(v, params)
#     v̇ = zero(v)
#     rhs!(v̇, v, params.ind, params.v_all, params.dist, params.n, params.u, params.eps, params.dS, params.ds_sum, params.ds_v_sum, params.mlb)
#     return v̇
# end

function f!(f::AbstractArray{T}, vn::AbstractArray{T}, vp, params,
        Δt, mlb::MetriplecticLenardBernstein) where {T}
    v_midpoint = mlb.cache[T].v
    v_midpoint .= (vn .+ vp) ./ 2

    collisional_vectorfield!(f, v_midpoint, params, mlb)

    f .*= Δt
    f .-= (vn .- vp)

    return nothing
end

# function IM_update!(y_new, y_new_guess, yₙ, f, Δt, params, ::MetriplecticLenardBernstein)
#     # println("define F!")
#     F!(g, x) = IM_rule!(g, x, yₙ, f, Δt, params)
#     obj = zero(yₙ)
#     # opt = Options(f_abstol = 1e-14)
#     # opt = Options(x_reltol = 1e-12, f_reltol = 1e-10, f_abstol = 1e-12)
#     # println("construct NewtonSolver")
#     nl = NewtonSolver(y_new_guess, obj)
#     # nl = NewtonSolver(y_new_guess, obj, config = opt)
#     # println("Newton solve")
#     SimpleSolvers.solve!(y_new, F!, nl)

#     return y_new
# end

function Picard_iterate_over_particles(dv::AbstractArray{ST}, vn::AbstractArray{ST},
        vn_minus_one::AbstractArray{ST}, dv_history, ti, t, Δt, m, β,
        abstol, reltol, mlb::MetriplecticLenardBernstein) where {ST}

    # err = 1
    f = 1
    # j = 0

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

    probN = NonlinearProblem{true}((f, v, p) -> f!(f, v, vn, params, Δt, mlb), v_prev)

    # NonlinearSolve.jl using Picard w/ anderson acceleration
    @time sol = NonlinearSolve.solve(probN,
        NonlinearSolve.NLsolveJL(; method = :anderson, m = m, beta = β);
        abstol = abstol, reltol = reltol, show_trace = Val(true))

    dv_history[:, 2] .= dv_history[:, 1]
    # dv_history[:, 1] .= dv
    collisional_vectorfield!(view(dv_history, :, 1), sol.u, params, mlb)

    return sol
    # return v_new, v_prev, j, err, f 

end

# function Picard_iterate_over_particles(dv::AbstractArray{ST}, vn::AbstractArray{ST}, vn_minus_one::AbstractArray{ST}, dv_history, ti, t,  Δt, max_iters, tol, ftol, mlb::MetriplecticLenardBernstein) where {ST}

#     err = 1
#     f = 1
#     j = 0

#     # set up vectors for storing intermediates
#     v_new = copy(vn)
#     v_prev = copy(vn)

#     dist = mlb.dist
#     ent = mlb.entropy
#     cache = mlb.cache[ST]
#     sdist = cache.sdist

#     vnew_vec = zeros(1) # 1-vector for passing to SimpleSolvers as it does not support scalars
#     fvec = zero(vn)

#     params = (dist = dist, ent = ent)

#     # use Hermite extrapolation to get an initial guess
#     if ti ≥ 4
#         Extrapolators.extrapolate!(t - 2Δt, vn_minus_one, dv_history[:, 2], t - Δt, vn, dv_history[:, 1], t, v_prev, Extrapolators.HermiteExtrapolation())
#     else
#         problemGNI = GeometricEquations.ODEProblem((v̇,t,v, params) -> collisional_vectorfield!(v̇,v,params,mlb), (t, t+Δt), Δt, vn; parameters = params)
#         Extrapolators.extrapolate!(t - Δt, vn, t, v_prev, problemGNI, Extrapolators.MidpointExtrapolation(5))
#     end

#     ## CHECK how good the initial guess is (compute midpoint using initial guess and then compute RHS of fixed point map using midpoint)
#     v_midpoint = (v_prev .+ vn)./2

#     projection(v_midpoint, mlb.dist, sdist)

#     compute_J!(cache.J, sdist, mlb)

#     compute_dS!(cache.dS, cache.J, v_midpoint, sdist, mlb)

#     # compute_dS_discrete_gradient!(cache.dS_dg, cache.dS, v_prev, vn, mlb)

#     n, u, eps = compute_moments(v_midpoint, mlb.dist, mlb)

#     # @show ds_sum = sum(cache.dS)
#     # @show ds_v_sum = dot(v_prev, cache.dS)
#     ds_sum = sum(cache.dS)
#     ds_v_sum = dot(v_prev, cache.dS)

#     for i in axes(dist.particles.v, 2)
#         # @show i
#         params_i = (ind = i, v_all = v_prev, dist = dist, n = n, u = u, eps = eps, dS = cache.dS, ds_sum = ds_sum, ds_v_sum = ds_v_sum, mlb)
#         # IM_update!(vnew_vec, [v_prev[i]], [vn[i]], rhs, Δt, params_i, mlb)
#         # v_new[i] = vnew_vec[1]
#         IM_rule!(view(fvec, i), [v_new[i]], [vn[i]], rhs, Δt, params_i)
#     end

#     @show f = norm(fvec) # this is the fixed point error of the initial guess 

#     # Now enter the Picard loop
#     while err > tol && f > ftol && j < max_iters
#         @show j

#         v_midpoint = (v_prev .+ vn)./2

#         projection(v_midpoint, mlb.dist, sdist)

#         compute_J!(cache.J, sdist, mlb)

#         compute_dS!(cache.dS, cache.J, v_midpoint, sdist, mlb)

#         # compute_dS_discrete_gradient!(cache.dS_dg, cache.dS, v_prev, vn, mlb)

#         n, u, eps = compute_moments(v_midpoint, mlb.dist, mlb)
#         ds_sum = sum(cache.dS)
#         ds_v_sum = dot(v_midpoint, cache.dS)

#         for i in axes(dist.particles.v, 2)
#             # @show i
#             params_i = (ind = i, v_all = v_prev, dist = dist, n = n, u = u, eps = eps, dS = cache.dS, ds_sum = ds_sum, ds_v_sum = ds_v_sum, mlb)
#             IM_update!(vnew_vec, [v_prev[i]], [vn[i]], rhs, Δt, params_i, mlb)
#             v_new[i] = vnew_vec[1]
#         end

#         v_midpoint = (v_new .+ vn)./2

#         projection(v_midpoint, mlb.dist, sdist)

#         compute_J!(cache.J, sdist, mlb)

#         compute_dS!(cache.dS, cache.J, v_midpoint, sdist, mlb)

#         # compute_dS_discrete_gradient!(cache.dS_dg, cache.dS, v_new, vn, mlb)

#         n, u, eps = compute_moments(v_midpoint, mlb.dist, mlb)

#         ds_sum = sum(cache.dS)
#         ds_v_sum = dot(v_midpoint, cache.dS)

#         # check norm of f
#         println("computing f")        
#         for i in axes(dist.particles.v, 2)

#             params_i = (ind = i, v_all = v_prev, dist = dist, n = n, u = u, eps = eps, dS = cache.dS, ds_sum = ds_sum, ds_v_sum = ds_v_sum, mlb)

#             IM_rule!(view(fvec, i), [v_new[i]], [vn[i]], rhs, Δt, params_i)

#             # dv[i] = rhs(v_midpoint[i], params_i)

#         end
#         # @show Δt * dot(vn, dv) + Δt^2 / 2 * norm(dv)^2

#         # collisional_vectorfield!(dv, (v_prev .+ vn)./2, params, mlb)
#         # v_new .= vn .+ Δt * dv

#         # @show err_sqeuc = sqeuclidean(v_new, v_prev) 
#         @show err = norm((v_new .- v_prev)./v_new)

#         # check rhs of fixed point map
#         # collisional_vectorfield!(dv, (v_new .+ vn)./2, params, mlb)
#         # @show f = norm(v_new .- vn .- Δt .* dv)
#         @show f = norm(fvec)

#         j += 1

#         if j == max_iters
#             println("WARNING: MAX Picard iterations reached")
#             @show err, f
#         end

#         v_prev .= v_new

#     end

#     dv_history[:, 2] .= dv_history[:, 1]
#     dv_history[:, 1] .= dv

#     println("ITERATIONS FINISHED")
#     # @show j

#     # v_midpoint = (v_prev .+ vn)./2

#     # projection(v_midpoint, mlb.dist, sdist)

#     # compute_J!(cache.J, sdist, mlb)

#     # compute_dS!(cache.dS, cache.J, v_midpoint, sdist, mlb)

#     # # compute_dS_discrete_gradient!(cache.dS_dg, cache.dS, v_prev, vn, mlb)

#     # n, u, eps = compute_moments(v_midpoint, mlb.dist, mlb)

#     # ds_sum = sum(cache.dS_dg)
#     # ds_v_sum = dot(v_prev, cache.dS_dg)

#     return v_new
#     # return v_new, v_prev, j, err, f 

# end
