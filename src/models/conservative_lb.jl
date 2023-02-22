struct ConservativeLenardBernstein{XD, VD, DT <: DistributionFunction{XD,VD}, ET <: Entropy, T} <: VlasovModel
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 
    
    function ConservativeLenardBernstein(dist::DistributionFunction{XD,VD}, ent::Entropy; ν::T=1.) where {XD, VD, T}
        new{XD, VD, typeof(dist), typeof(ent), T}(dist, ent, ν)
    end
end

function compute_coefficients(distribution::SplineDistribution, particle_dist::ParticleDistribution, vp::AbstractArray{VT}) where {VT}

    n, nu, neps = compute_f_densities(distribution, vp)
    B1, B2 = compute_df_densities(distribution, vp)
    B1 *= -1 
    B2 *= -1 
    A1 = (neps * B1 - nu * B2) / (n * neps - (nu)^2)
    A2 = -  (nu * B1 - n * B2) / (n * neps - (nu)^2)

    return A1, A2
end

# RHS function for solving collisions using DifferentialEquations.jl
 function CLB_rhs!(v̇, v::AbstractVector{ST}, params, t) where {ST}

    dist = params.model.ent.cache[ST]

    fs = projection(v, params.idist, dist)

    dfdv = Derivative(1) * fs 

    A = compute_coefficients(dist, params.idist, v)

    v̇ .= -params.ν .* (dfdv.(v) .+ ( A[1] .+ A[2] .* v) .* fs.(v))

 end

 function CLB_rhs_GI!(v, t, q::AbstractArray{ST}, params) where {ST}

    dist = params.model.ent.cache[ST]

    fs = projection(q, params.idist, dist)

    dfdv = Derivative(1) * fs 

    A = compute_coefficients(dist, params.idist, q)

    v .= -params.ν .* (dfdv.(q) .+ ( A[1] .+ A[2] .* q) .* fs.(q))

 end

# used for plotting
 function CLB_rhs(v::AbstractVector{ST}, params, fs::Spline) where {ST}

    dist = params.model.ent.cache[ST]

    dfdv = Derivative(1) * fs 

    A = compute_coefficients(dist, params.idist, v)

    v̇ = -params.ν .* (dfdv.(v) .+( A[1] .+ A[2] * v) .* fs.(v))

    return v̇
 end

 

 function DiffEqIntegrator(model::ConservativeLenardBernstein{1,1}, tspan::Tuple, tstep::Real)
    # parameters for computing vector field
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # u0 = copy(model.dist.particles.v[1,:])
    # construct DifferentialEquations ODEProblem
    equ = DifferentialEquations.ODEProblem(
        CLB_rhs!,
        copy(model.dist.particles.v[1,:]),
        tspan,
        params
    )

    # choose integrator
    int = DifferentialEquations.TRBDF2()
    # int = DifferentialEquations.Trapezoid()

    DiffEqIntegrator(model, equ, int, tstep)
 end


function GeometricIntegrator(model::ConservativeLenardBernstein{1,1}, tspan::Tuple, tstep::Real)
    # collect parameters
    # params = (ϕ = model.potential, model = model)
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # create geometric problem
    equ = GeometricEquations.ODEProblem(
            CLB_rhs_GI!,
            tspan, tstep, copy(model.dist.particles.v[1,:]);
            parameters = params)

    # create integrator
    int = Integrators.Integrator(equ, Integrators.RK438())
    # int = Integrators.Integrator(equ, Integrators.CrankNicolson())

    # put together splitting method
    GeometricIntegrator(model, equ, int)
end
