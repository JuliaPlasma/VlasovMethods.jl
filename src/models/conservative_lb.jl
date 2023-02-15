struct ConservativeLenardBernstein{XD, VD, DT <: DistributionFunction{XD,VD}, ET <: Entropy, T} <: VlasovModel
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 
    
    function ConservativeLenardBernstein(dist::DistributionFunction{XD,VD}, ent::Entropy; ν::T=1.) where {XD, VD, T}
        new{XD, VD, typeof(dist), typeof(ent), T}(dist, ent, ν)
    end
end

# helper function to compute coefficients A₁ and A₂ in the expansion of the conservative operator, C[f](v) = ∂ᵥ(∂ᵥf + A₁f + A₂vf).
function compute_coefficients(distribution::SplineDistribution, particle_dist::ParticleDistribution)

    D = zeros(eltype(distribution), 2, 2)
    D[1,1], D[1,2], D[2,2] = compute_v_densities(distribution, particle_dist)
    D[2,1] = D[1,2]

    B = zeros(eltype(distribution), 2)
    B .= compute_v_densities(distribution, particle_dist; isDerivative=true)
    B .*= -1
    
    return D\B
end

# RHS function for solving collisions using DifferentialEquations.jl
 function CLB_rhs!(v̇, v::AbstractVector{ST}, params, t) where {ST}

    dist = params.model.ent.cache[ST]

    fs = projection(v, params.idist, dist)

    dfdv = Derivative(1) * fs 

    A = compute_coefficients(dist, params.idist)

    v̇ .= -params.ν .* (dfdv.(v) .+ ( A[1] .+ A[2] * v) .* fs.(v))

 end

# used for plotting
 function CLB_rhs(v::AbstractVector{ST}, params, fs::Spline) where {ST}

    dist = params.model.ent.cache[ST]

    dfdv = Derivative(1) * fs 

    A = compute_coefficients(dist, params.idist)

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
    # int = DifferentialEquations.TRBDF2()
    int = DifferentialEquations.Trapezoid()

    DiffEqIntegrator(model, equ, int, tstep)
 end