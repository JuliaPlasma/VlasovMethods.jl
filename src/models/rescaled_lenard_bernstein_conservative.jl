# The particles, the spline and the work arrays of a rescaled conservative Lenard–Bernstein
# right-hand side in the element type `T` of a state, as `LBCache` holds them. A cache belongs
# to one task.
struct RCLBCache{T, PT <: ParticleDistribution, ST <: SplineDistribution{T}, VT} <: Cache{T}
    pdist::PT
    sdist::ST

    rhs::Vector{T}
    values::VT

    function RCLBCache{T}(pdist, sdist) where {T}
        values = _local_buffers(sdist.basis, T)
        new{T, typeof(pdist), typeof(sdist), typeof(values)}(
            pdist, sdist, zeros(T, nbasis(sdist.basis)), values)
    end
end

function RCLBCache(pdist::ParticleDistribution{T}, sdist::SplineDistribution{T}) where {T}
    RCLBCache{T}(pdist, sdist)
end

@collision_cache RCLBCache buffers

struct RescaledConservativeLenardBernstein{
    XD, VD, DT <: DistributionFunction{<:Any, XD, VD}, ET <: Entropy, T, CT <: CacheDict} <:
       CollisionOperator
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 

    cache::CT

    function RescaledConservativeLenardBernstein(
            dist::DistributionFunction{<:Any, XD, VD}, ent::Entropy; ν::T = 1.0) where {
            XD, VD, T}
        cache = CacheDict(RCLBCache(dist, ent.dist))
        new{XD, VD, typeof(dist), typeof(ent), T, typeof(cache)}(dist, ent, ν, cache)
    end
end

@doc raw"""
    compute_coefficients_rclb(distribution, particle_dist, vp[, values])

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
[`compute_coefficients`](@ref), and `values` is the buffer of the point evaluations as there.
"""
function compute_coefficients_rclb(
        distribution::SplineDistribution, particle_dist::ParticleDistribution,
        vp::AbstractArray{VT},
        values::AbstractVector = _local_buffers(
            distribution.basis, promote_type(eltype(distribution.spline), VT))) where {VT}
    w = particle_dist.particles.w
    fs = distribution.spline

    n = zero(VT)
    nu = zero(VT)
    neps = zero(VT)
    # the sums of `f_s'/f_s` take the type of the buffer where it is wider than that of the state
    BT = promote_type(VT, eltype(values))
    B1 = zero(BT)
    B2 = zero(BT)

    # the error carried by each running sum, in the type of the sum
    εn = zero(VT)
    εnu = zero(VT)
    εneps = zero(VT)
    εB1 = zero(BT)
    εB2 = zero(BT)

    for α in eachindex(vp)
        wα = VT(w[1, α])
        vα = vp[α]
        fα = evaluate!(values, fs, vα)

        fα > 0 || throw(DomainError(fα,
            "the projected distribution is non-positive at v = $(vα); see " *
            "the positivity note in `compute_coefficients`"))

        n, εn = compensated_summation(n, wα, εn)
        nu, εnu = compensated_summation(nu, wα * vα, εnu)
        neps, εneps = compensated_summation(neps, wα * vα^2, εneps)

        g = evaluate!(values, fs, vα, 1) / fα
        B1, εB1 = compensated_summation(B1, wα * g, εB1)
        B2, εB2 = compensated_summation(B2, wα * vα * g, εB2)
    end

    den = nu * B1 - n * B2
    A1 = (n * neps - nu^2) / den
    A2 = (nu * B2 - neps * B1) / den

    if !isfinite(A1) || !isfinite(A2)
        # `isfinite` rather than `isnan`, for the reason given in `compute_coefficients`.
        throw(ArgumentError(
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
    (; sdist, rhs, values) = params.model.cache[ST]

    fs = _projection!(rhs, values, q, params.idist, sdist)

    A1, A2 = compute_coefficients_rclb(sdist, params.idist, q, values)

    for α in eachindex(v, q)
        qα = q[α]
        v[α] = -params.ν *
               (A1 * evaluate!(values, fs, qα, 1) / evaluate!(values, fs, qα) + (A2 + qα))
    end
    return v
end

# used for plotting
function RCLB_rhs(v::AbstractVector{ST}, params, fs::Spline) where {ST}
    dist = params.model.cache[ST].sdist

    dfdv = derivative(fs)

    A = compute_coefficients_rclb(dist, params.idist, v)

    return -params.ν .* (A[1] .* dfdv.(v) ./ fs.(v) .+ (A[2] .+ v))
end

# The velocity space is 1D, in any spatial dimension: a scalar velocity and one weight per particle.
function GeometricIntegrator(
        model::RescaledConservativeLenardBernstein{<:Any, 1}, tspan::Tuple, tstep::Real,
        method = Integrators.ImplicitMidpoint())
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
    int = Integrators.GeometricIntegrator(equ, method)
    # int = Integrators.Integrator(equ, Integrators.ImplicitMidpoint())
    # int = Integrators.Integrator(equ, Integrators.RK438())
    # int = Integrators.Integrator(equ, Integrators.CrankNicolson())

    GeometricIntegrator(model, equ, int)
    # return int
end
