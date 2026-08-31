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

        v = zeros(N)

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

function Cache(AT, mlb::MetriplecticLenardBernstein)
    MLBCache{AT}(mlb.pdist, similar(AT, mlb.sdist))
end
function CacheType(AT, mlb::MetriplecticLenardBernstein)
    MLBCache{AT, typeof(mlb.pdist), similar_type(AT, mlb.sdist)}
end

function compute_J!(J, sdist::SplineDistribution{T, 1, 1}, ::MetriplecticLenardBernstein) where {T}
    # knot_list = BSplineKit.knots(sdist.basis)
    # for i in eachindex(J)
    #     J[i], _ = quadgk(x -> (1 + 0.5 * log((sdist.spline(x))^2))*sdist.basis[i](x), knot_list[1], knot_list[end], atol = 1e-14)
    # end
    J .= BSplineKit.galerkin_projection(x -> 1 + 0.5 * log((sdist.spline(x))^2), sdist.basis)
    ldiv!(sdist.mass_fact, J)
end

function compute_dS!(dS, J, v::AbstractArray{ST}, sdist::SplineDistribution{ST, 1, 1},
        ::MetriplecticLenardBernstein) where {ST}
    dS .= 0
    for i in eachindex(dS)
        ilast, bs = BSplineKit.evaluate_all(sdist.basis, v[i], BSplineKit.Derivative(1))
        for (δi, bi) in pairs(bs)
            if abs(bi) > eps()
                j = ilast + 1 - δi
                dS[i] += J[j] * bi / length(dS)
                if isnan(dS[i])
                    println("NANs")
                end
            else
                continue
            end
        end
    end
end

function compute_entropy(f, mlb::MetriplecticLenardBernstein)
    k = BSplineKit.knots(mlb.entropy.dist.basis)
    S, err = quadgk(x -> f(x) * 0.5 * log((f(x))^2), k[1], k[end], atol = 1e-14)
    return S, err
end

function compute_dS_discrete_gradient!(dS_dg, dS_midpoint, v_new::AbstractArray{ST},
        vn::AbstractArray{ST}, mlb::MetriplecticLenardBernstein) where {ST}
    sdist = mlb.cache[ST].sdist

    projection(vn, mlb.dist, sdist)
    @show S_n, err_n = compute_entropy(x -> sdist.spline(x), mlb)

    projection(v_new, mlb.dist, sdist)
    @show S_n_plus_1, err_n_plus_1 = compute_entropy(x -> sdist.spline(x), mlb)

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

function rhs!(
        v̇::AbstractArray{ST}, v::AbstractArray{ST}, pdist::ParticleDistribution, n, u, eps,
        dS::AbstractArray{ST}, dS_sum, dS_v_sum, ::MetriplecticLenardBernstein) where {ST}
    # dS_sum = sum(dS)
    v̇ .= - n / pdist.particles.w[1] * (eps - u^2) .* dS .+ (eps .- u .* v) * dS_sum .+
         (v .- u) * dS_v_sum
end

function rhs_downstairs_factor!(
        v̇::AbstractArray{ST}, v::AbstractArray{ST}, pdist::ParticleDistribution, n, u, eps,
        dS::AbstractArray{ST}, dS_sum, dS_v_sum, ::MetriplecticLenardBernstein) where {ST}
    # dS_sum = sum(dS)
    # v̇ .= - one(ST) / pdist.particles.w[1] .* dS .+ (eps .- u .* v) ./ (eps - u^2) * dS_sum .+ (v .- u) ./ (eps - u^2) * dS_v_sum
    v̇ .= - n / pdist.particles.w[1] .* dS .+ (eps .- u .* v) ./ (eps - u^2) * dS_sum .+
         (v .- u) ./ (eps - u^2) * dS_v_sum
end

function collisional_vectorfield!(v̇::AbstractArray{ST}, v::AbstractArray{ST}, params,
        mlb::MetriplecticLenardBernstein) where {ST}
    cache = mlb.cache[ST]

    sdist = cache.sdist

    projection(v, mlb.dist, sdist)

    compute_J!(cache.J, sdist, mlb)

    compute_dS!(cache.dS, cache.J, v, sdist, mlb)

    ds_sum = sum(cache.dS)
    ds_v_sum = dot(v, cache.dS)

    n, u, eps = compute_moments(v, mlb.dist, mlb)

    # rhs!(v̇, v, mlb.dist, n, u, eps, cache.dS, ds_sum, ds_v_sum, mlb)
    rhs_downstairs_factor!(v̇, v, mlb.dist, n, u, eps, cache.dS, ds_sum, ds_v_sum, mlb)
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
        Extrapolators.extrapolate!(t - 2Δt, vn_minus_one, dv_history[:, 2], t - Δt, vn,
            dv_history[:, 1], t, v_prev, Extrapolators.HermiteExtrapolation())
    else
        problemGNI = GeometricEquations.ODEProblem(
            (v̇, t, v, params) -> collisional_vectorfield!(v̇, v, params, mlb),
            (t, t+Δt), Δt, vn; parameters = params)
        Extrapolators.extrapolate!(
            t - Δt, vn, t, v_prev, problemGNI, Extrapolators.MidpointExtrapolation(5))
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
