
struct VlasovPoisson{XD, VD, DT <: DistributionFunction{XD,VD}, PT <: Potential, RT <: AbstractVector} <: VlasovModel
    distribution::DT
    potential::PT
    rhs::RT

    function VlasovPoisson(dist::DistributionFunction{XD,VD}, potential) where {XD,VD}
        rhs = zero(potential.coefficients)
        new{XD, VD, typeof(dist), typeof(potential), typeof(rhs)}(dist, potential, rhs)
    end
end



function update_potential!(model::VlasovPoisson)
    projection!(model.rhs, model.potential, model.distribution)
    PoissonSolvers.update!(model.potential, model.rhs)
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

# splitting fields
function v_advection!(ż, t, z, params)
    for i in axes(ż, 2)
        ż[1,i] = z[2,i]
        ż[2,i] = 0
    end
end

function v_lorentz_force!(ż, t, z, params)
    update_potential!(params.model)
    for i in axes(ż, 2)
        ż[1,i] = 0
        ż[2,i] = - params.ϕ(z[1,i], Derivative(1))
    end
end

function s_advection!(z, t, z̄, t̄, params)
    for i in axes(z, 2)
        z[1,i] = z̄[1,i] + (t-t̄) * z̄[2,i]
        z[2,i] = z̄[2,i]
    end
end

function s_lorentz_force!(z, t, z̄, t̄, params)
    update_potential!(params.model)
    for i in axes(z, 2)
        z[1,i] = z̄[1,i]
        z[2,i] = z̄[2,i] - (t-t̄) * params.ϕ(z̄[1,i], Derivative(1))
    end
end


function SplittingMethod(model::VlasovPoisson{1,1}, tspan::Tuple, tstep::Real)
    # collect parameters
    params = (ϕ = model.potential, model = model)

    # create geometric problem
    equ = GeometricEquations.SODEProblem(
            (v_advection!, v_lorentz_force!),
            (s_advection!, s_lorentz_force!),
            tspan, tstep, copy(model.distribution.particles.z);
            parameters = params)

    # create integrator
    int = Integrators.Integrator(equ, Integrators.Strang())

    # put together splitting method
    SplittingMethod(model, equ, int)
end
