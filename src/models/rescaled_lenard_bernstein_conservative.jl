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

@doc raw"""
    compute_coefficients_rclb(distribution, particle_dist, vp)

The coefficients of the **rescaled** parametrisation, in which the coefficient of ``v`` is
normalised to one:

```math
\dot{v}_\alpha = - \nu \left( A_1 \frac{f_s'(v_\alpha)}{f_s(v_\alpha)} + A_2 + v_\alpha \right) .
```

This is the conservative operator with time rescaled by ``\varepsilon_h - u_h^2``, so that
``A_1^{\mathrm{resc}} = 1 / A_2^{\mathrm{cons}} = \varepsilon_h - u_h^2`` and
``A_2^{\mathrm{resc}} = A_1^{\mathrm{cons}} / A_2^{\mathrm{cons}} = -u_h``. It is *not* the
fixed-coefficient cumulant-scaling experiment of the manuscript appendix, which sets the
coefficients by hand rather than solving for them.

The particle weights appear in every sum, for the reason given at
[`compute_coefficients`](@ref).
"""
function compute_coefficients_rclb(
        distribution::SplineDistribution, particle_dist::ParticleDistribution,
        vp::AbstractArray{VT}) where {VT}
    w = particle_dist.particles.w
    fs = distribution.spline
    dfs = derivative(fs)

    n = zero(VT)
    nu = zero(VT)
    neps = zero(VT)
    B1 = zero(VT)
    B2 = zero(VT)

    for α in eachindex(vp)
        wα = w[1, α]
        vα = vp[α]
        fα = fs(vα)

        fα > 0 || throw(ErrorException(
            "the projected distribution is non-positive, f_s(v) = $(fα) at v = $(vα); see " *
            "the positivity note in `compute_coefficients`"))

        n += wα
        nu += wα * vα
        neps += wα * vα^2

        g = dfs(vα) / fα
        B1 += wα * g
        B2 += wα * vα * g
    end

    den = nu * B1 - n * B2
    A1 = (n * neps - nu^2) / den
    A2 = (nu * B2 - neps * B1) / den

    if isnan(A1) || isnan(A2)
        throw(ErrorException(
            "the rescaled coefficient system is singular: n = $(n), nu = $(nu), " *
            "neps = $(neps), denominator $(den)."))
    end

    return A1, A2
end

# `compute_coefficients_rclb`, not `compute_coefficients`. Both call sites here used the
# *conservative* coefficients, which belong to a different parametrisation: A₁ᶜᵒⁿˢ = -u/σ²
# was being used as the multiplier of f'/f and A₂ᶜᵒⁿˢ = 1/σ² as the constant drift, giving a
# vector field with no conservation property at all.
function RCLB_rhs_GI!(v, t, q::AbstractArray{ST}, params) where {ST}
    dist = params.model.cache[ST].sdist

    fs = projection(q, params.idist, dist)
    dfdv = derivative(fs)

    A = compute_coefficients_rclb(dist, params.idist, q)

    v .= -params.ν .* (A[1] .* dfdv.(q) ./ fs.(q) .+ (A[2] .+ q))
end

# used for plotting
function RCLB_rhs(v::AbstractVector{ST}, params, fs::Spline) where {ST}
    dist = params.model.cache[ST].sdist

    dfdv = derivative(fs)

    A = compute_coefficients_rclb(dist, params.idist, v)

    return -params.ν .* (A[1] .* dfdv.(v) ./ fs.(v) .+ (A[2] .+ v))
end

function GeometricIntegrator(model::RescaledConservativeLenardBernstein, tspan::Tuple, tstep::Real)
    # collect parameters
    # params = (ϕ = model.potential, model = model)
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # create geometric problem
    equ = GeometricEquations.ODEProblem(
        # `RCLB_rhs_GI!`, not `CLB_rhs_GI!`. This constructor wired the *conservative*
        # right-hand side, so every run built as a RescaledConservativeLenardBernstein
        # actually integrated the unrescaled operator and `RCLB_rhs_GI!` was dead code.
        RCLB_rhs_GI!,
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
