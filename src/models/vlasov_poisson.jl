
struct VlasovPoisson{XD, VD, DT <: DistributionFunction{XD,VD}, PT <: Potential} <: VlasovModel
    distribution::DT
    potential::PT

    function VlasovPoisson(dist::DistributionFunction{XD,VD}, potential) where {XD,VD}
        new{XD, VD, typeof(dist), typeof(potential)}(dist, potential)
    end
end


function update_potential!(model::VlasovPoisson)
    projection!(model.potential, model.distribution)
    PoissonSolvers.update!(model.potential)
end


####################################################
# Define Splitting Method for Vlasov-Poisson Model #
####################################################

# vector field
function lorentz_force!(ż, t, z, params)
    update_potential!(params.model)
    for i in axes(ż, 2)
        ż[1,i] = z[2,i]
        ż[2,i] = - params.ϕ(z[1,i], Derivative(1))
    end
end

###########################################################
# Vlasov-Poisson 1D1V splitting fields for particles      #
###########################################################

# Vector field for advection
function v_advection!(ż, t, z, params)
    for i in axes(ż, 2)
        ż[1,i] = z[2,i]
        ż[2,i] = 0
    end
end

# Vector field for acceleration
function v_acceleration!(ż, t, z, params)
    update_potential!(params.model)
    for i in axes(ż, 2)
        ż[1,i] = 0
        ż[2,i] = - params.ϕ(z[1,i], Derivative(1))
    end
end

# Solution for advection
function s_advection!(z, t, z̄, t̄, params)
    for i in axes(z, 2)
        z[1,i] = z̄[1,i] + (t-t̄) * z̄[2,i]
        z[2,i] = z̄[2,i]
    end
end

# Solution for Lorentz force
function s_acceleration!(z, t, z̄, t̄, params)
    update_potential!(params.model)
    for i in axes(z, 2)
        z[1,i] = z̄[1,i]
        z[2,i] = z̄[2,i] - (t-t̄) * params.ϕ(z̄[1,i], Derivative(1))
    end
end

# Constructor for a splitting method from GeometricIntegrators
# The problem is setup such that one solution step pushes all particles.
# While this allows for a simple implementation, it is not well-suited
# for parallelisation.
function SplittingMethod(model::VlasovPoisson{1, 1, <: ParticleDistribution}, tspan::Tuple, tstep::Real)
    # collect parameters
    params = (ϕ = model.potential, model = model)

    # create geometric problem
    equ = GeometricEquations.SODEProblem(
            (v_advection!, v_acceleration!),
            (s_advection!, s_acceleration!),
            tspan, tstep, copy(model.distribution.particles.z);
            parameters = params)

    # create integrator
    int = GeometricIntegrators.GeometricIntegrator(equ, GeometricIntegrators.Strang())

    # put together splitting method
    SplittingMethod(model, equ, int)
end
