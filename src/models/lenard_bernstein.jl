struct LenardBernstein{XD, VD, DT <: DistributionFunction{XD,VD}, ET <: Entropy, T} <: CollisionOperator
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 
    
    function LenardBernstein(dist::DistributionFunction{XD,VD}, ent::Entropy; ν::T=1.) where {XD, VD, T}
        new{XD, VD, typeof(dist), typeof(ent), T}(dist, ent, ν)
    end
end

# function update_distribution!(model::LenardBernstein, v_new::VT) where {VT}
#     model.dist.particles.v .= v_new'
# end

function update_entropy!(model::LenardBernstein)
    projection!(model.dist, model.ent.dist)
end

# RHS function for solving collisions using DifferentialEquations.jl
 function LB_rhs!(v̇, v::AbstractArray{ST}, params, t) where {ST}

    dist = params.model.ent.cache[ST]

    fs = projection(v, params.idist, dist)

    dfdv = Derivative(1) * fs # TODO: does this belong here?

    v̇ .= -params.ν .* (dfdv.(v) .+ v .* fs.(v))

 end

 function LB_rhs_GI!(v, t, q::AbstractArray{ST}, params) where {ST}
    LB_rhs!(v, q, params, t)
 end

# used for plotting
 function LB_rhs(v, params, fs::Spline)

    dfdv = Derivative(1) * fs # TODO: does this belong here?

    v̇ = -params.ν .* (dfdv.(v) .+ v .* fs.(v))

    return v̇
 end

 

 function DiffEqIntegrator(model::LenardBernstein{1,1}, tspan::Tuple, tstep::Real)
    # parameters for computing vector field
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # u0 = copy(model.dist.particles.v[1,:])
    # construct DifferentialEquations ODEProblem
    equ = DifferentialEquations.ODEProblem(
        LB_rhs!,
        copy(model.dist.particles.v[1,:]),
        tspan,
        params
    )

    # choose integrator
    int = DifferentialEquations.TRBDF2()
    # int = DifferentialEquations.Trapezoid()

    DiffEqIntegrator(model, equ, int, tstep)
 end


function GeometricIntegrator(model::LenardBernstein{1,1}, tspan::Tuple, tstep::Real)
    # collect parameters
    # params = (ϕ = model.potential, model = model)
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # create geometric problem
    equ = GeometricEquations.ODEProblem(
            LB_rhs_GI!,
            tspan, tstep, copy(model.dist.particles.v[1,:]);
            parameters = params)

    # create integrator
    int = GeometricIntegrators.GeometricIntegrator(equ, GeometricIntegrators.RK438())
    # int = GeometricIntegrators.GeometricIntegrator(equ, GeometricIntegrators.CrankNicolson())

    # put together splitting method
    GeometricIntegrator(model, equ, int)
end
