struct LenardBernstein{XD, VD, DT <: DistributionFunction{XD, VD}, ET <: Entropy, T} <:
       CollisionOperator
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 

    function LenardBernstein(dist::DistributionFunction{XD, VD}, ent::Entropy; ν::T = 1.0) where {
            XD, VD, T}
        new{XD, VD, typeof(dist), typeof(ent), T}(dist, ent, ν)
    end
end

# function update_distribution!(model::LenardBernstein, v_new::VT) where {VT}
#     model.dist.particles.v .= v_new'
# end

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

function LB_rhs_GI!(v, t, q::AbstractArray{ST}, params) where {ST}
    LB_rhs!(v, q, params, t)
end

# used for plotting
function LB_rhs(v, params, fs::Spline)
    dfdv = derivative(fs)

    return -params.ν .* (dfdv.(v) ./ fs.(v) .+ v)
end

function DiffEqIntegrator(model::LenardBernstein{1, 1}, tspan::Tuple, tstep::Real)
    # parameters for computing vector field
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # u0 = copy(model.dist.particles.v[1,:])
    # construct DifferentialEquations ODEProblem
    equ = DifferentialEquations.ODEProblem(
        LB_rhs!,
        copy(model.dist.particles.v[1, :]),
        tspan,
        params
    )

    # choose integrator
    int = DifferentialEquations.TRBDF2()
    # int = DifferentialEquations.Trapezoid()

    DiffEqIntegrator(model, equ, int, tstep)
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
