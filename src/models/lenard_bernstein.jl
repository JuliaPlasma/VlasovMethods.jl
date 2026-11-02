# The particles and the spline of a Lenard–Bernstein right-hand side in the element type `T` of
# a state, with its work arrays: `rhs`, the load vector of the projection, and `values`, the
# basis functions that are non-zero at one point, which the projection and each point
# evaluation of the spline use in turn. A cache belongs to one task: two tasks that evaluate at
# once on one cache share its buffers.
struct LBCache{T, PT <: ParticleDistribution, ST <: SplineDistribution{T}} <: Cache{T}
    pdist::PT
    sdist::ST

    rhs::Vector{T}
    values::Vector{T}

    function LBCache{T}(pdist, sdist) where {T}
        new{T, typeof(pdist), typeof(sdist)}(
            pdist, sdist, zeros(T, nbasis(sdist.basis)), _local_buffers(sdist.basis, T))
    end
end

function LBCache(pdist::ParticleDistribution{T}, sdist::SplineDistribution{T}) where {T}
    LBCache{T}(pdist, sdist)
end

@collision_cache LBCache

struct LenardBernstein{
    XD, VD, DT <: DistributionFunction{<:Any, XD, VD}, ET <: Entropy, T, CT <: CacheDict} <:
       CollisionOperator
    dist::DT    # distribution function
    ent::ET     # entropy
    ν::T        # collision frequency

    cache::CT

    function LenardBernstein(
            dist::DistributionFunction{<:Any, XD, VD}, ent::Entropy; ν::T = 1.0) where {
            XD, VD, T}
        cache = CacheDict(LBCache(dist, ent.dist))
        new{XD, VD, typeof(dist), typeof(ent), T, typeof(cache)}(dist, ent, ν, cache)
    end
end

@doc raw"""
The plain Lenard-Bernstein right-hand side,

```math
\dot{v}_\alpha = - \nu \left( \frac{f_s'(v_\alpha)}{f_s(v_\alpha)} + v_\alpha \right) ,
```

i.e. the conservative operator of `eq:cons_LB` with the coefficients fixed at
``A_1 = 0``, ``A_2 = 1``. It conserves mass but not momentum or energy — that is the whole
point of the conservative variant.

!!! note "The division by `f_s` was missing"
    Both right-hand sides computed `-ν (f_s' + v f_s)`, which is the collisional **flux**
    ``F[f_s]``, not the advection coefficient. The two differ by a factor of ``f_s(v_\alpha)``:
    the equilibrium condition ``f_s'/f_s = -v`` is the same either way, so the fixed point was
    right, but particles in the tails — where ``f_s`` is small — were barely advected and the
    whole transient was wrong.
"""
function LB_rhs!(v̇, v::AbstractArray{ST}, params, t) where {ST}
    # `params.fdist`, not `params.model.ent.cache`: `CollisionEntropy` has a `dist` field and
    # no `cache` field, so the earlier line threw.
    dist = params.fdist

    fs = projection(v, params.idist, dist)
    dfdv = derivative(fs)

    v̇ .= -params.ν .* (dfdv.(v) ./ fs.(v) .+ v)
end

# The right-hand side of `LB_rhs!` with the buffers of the model's cache of the element type of
# `q`, which allocates nothing once that cache exists. `LB_rhs!` takes `params` without the
# model, and allocates the projection and the point evaluations.
function LB_rhs_GI!(v, t, q::AbstractArray{ST}, params) where {ST}
    (; sdist, rhs, values) = params.model.cache[ST]

    fs = _projection!(rhs, values, q, params.idist, sdist)

    for α in eachindex(v, q)
        qα = q[α]
        v[α] = -params.ν * (evaluate!(values, fs, qα, 1) / evaluate!(values, fs, qα) + qα)
    end
    return v
end

# used for plotting
function LB_rhs(v, params, fs::Spline)
    dfdv = derivative(fs)

    return -params.ν .* (dfdv.(v) ./ fs.(v) .+ v)
end

function GeometricIntegrator(model::LenardBernstein{1, 1}, tspan::Tuple, tstep::Real)
    # collect parameters
    # params = (ϕ = model.potential, model = model)
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # create geometric problem
    equ = GeometricEquations.ODEProblem(
        LB_rhs_GI!,
        tspan, tstep, copy(model.dist.particles.v[1, :]);
        parameters = params)

    # create integrator
    int = Integrators.GeometricIntegrator(equ, Integrators.RK438())
    # int = Integrators.GeometricIntegrator(equ, Integrators.CrankNicolson())

    # put together splitting method
    GeometricIntegrator(model, equ, int)
end
