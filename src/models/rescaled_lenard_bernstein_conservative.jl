struct RCLBCache{T, PT <: ParticleDistribution, ST <: SplineDistribution{T}} <: Cache{T}
    pdist::PT
    sdist::ST

    function RCLBCache{T}(pdist, sdist) where {T}
        new{T, typeof(pdist), typeof(sdist)}(pdist, sdist)
    end
end

function RCLBCache(pdist::ParticleDistribution{T}, sdist::SplineDistribution{T}) where {T}
    RCLBCache{T}(pdist, sdist)
end

function Cache(AT, c::RCLBCache{DT, PT, ST}) where {DT, PT, ST}
    RCLBCache{AT}(c.pdist, similar(AT, c.sdist))
end
function CacheType(AT, c::RCLBCache{DT, PT, ST}) where {DT, PT, ST}
    RCLBCache{AT, PT, similar_type(AT, c.sdist)}
end

struct RescaledConservativeLenardBernstein{
    XD, VD, DT <: DistributionFunction{XD, VD}, ET <: Entropy, T, CT <: CacheDict} <:
       CollisionOperator
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 

    cache::CT

    function RescaledConservativeLenardBernstein(
            dist::DistributionFunction{XD, VD}, ent::Entropy; ν::T = 1.0) where {XD, VD, T}
        cache = CacheDict(RCLBCache(dist, ent.dist))
        new{XD, VD, typeof(dist), typeof(ent), T, typeof(cache)}(dist, ent, ν, cache)
    end
end

function Cache(AT, clb::RescaledConservativeLenardBernstein)
    RCLBCache{AT}(clb.dist, similar(AT, clb.ent.dist))
end
function CacheType(AT, clb::RescaledConservativeLenardBernstein)
    RCLBCache{AT, typeof(clb.dist), similar_type(AT, clb.ent.dist)}
end

function compute_coefficients_rclb(
        distribution::SplineDistribution, particle_dist::ParticleDistribution,
        vp::AbstractArray{VT}) where {VT}
    # function compute_coefficients(distribution::SplineDistribution, particle_dist::ParticleDistribution, vp::AbstractArray{VT}) where {VT}

    # compute discrete moments of fₚ
    n = length(vp)
    nu = sum(vp)
    neps = sum(vp .^ 2)

    # compute derivative of fₛ and required terms for A1 and A2
    dfs = Derivative(1) * distribution.spline
    B1_α = [one(VT)/distribution.spline(v_α)*dfs(v_α) for v_α in vp]
    B2_α = [v_α/distribution.spline(v_α)*dfs(v_α) for v_α in vp]

    B1 = sum(B1_α)
    B2 = sum(B2_α)

    # directly compute A1 and A2 coefficients 
    A1 = (n * neps - (nu)^2) / (nu * B1 - n * B2)
    A2 = (nu * B2 - neps * B1) / (nu * B1 - n * B2)

    # check for NaNs
    if isnan(A1) || isnan(A2)
        throw(ErrorException("NaNs in computation of A1 or A2, use a smaller timestep."))
    end

    return A1, A2
end

function RCLB_rhs_GI!(v, t, q::AbstractArray{ST}, params) where {ST}
    dist = params.model.cache[ST].sdist

    fs = projection(q, params.idist, dist)

    dfdv = BSplineKit.Derivative(1) * fs

    A = compute_coefficients(dist, params.idist, q)

    v .= -params.ν .* ((A[1] ./ fs.(q)) .* dfdv.(q) .+ (A[2] .+ q))
end

# used for plotting
function RCLB_rhs(v::AbstractVector{ST}, params, fs::BSplineKit.Spline) where {ST}
    dist = params.model.cache[ST].sdist

    dfdv = BSplineKit.Derivative(1) * fs

    A = compute_coefficients(dist, params.idist, v)
    v̇ = -params.ν .* ((A[1] ./ fs.(v)) .* dfdv.(v) .+ (A[2] .+ v))
    # v̇ = -params.ν .* ((one(ST) ./ fs.(v)) .* dfdv.(v) .+( A[1] .+ A[2] .* v))
    # v̇ = -params.ν .* (dfdv.(v) .+( A[1] .+ A[2] .* v) .* fs.(v))

    return v̇
end

function GeometricIntegrator(model::RescaledConservativeLenardBernstein, tspan::Tuple, tstep::Real) where {DT}
    # collect parameters
    # params = (ϕ = model.potential, model = model)
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # create geometric problem
    equ = GeometricEquations.ODEProblem(
        CLB_rhs_GI!,
        tspan, tstep, copy(model.dist.particles.v[1, :]);
        parameters = params)

    # create integrator
    int = Integrators.GeometricIntegrator(equ, Integrators.ImplicitMidpoint())
    # int = Integrators.Integrator(equ, Integrators.ImplicitMidpoint())
    # int = Integrators.Integrator(equ, Integrators.RK438())
    # int = Integrators.Integrator(equ, Integrators.CrankNicolson())

    # put together splitting method
    GeometricIntegrator(model, equ, int)
    # return int
end
