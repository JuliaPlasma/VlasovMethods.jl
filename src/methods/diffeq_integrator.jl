struct DiffEqIntegrator{MT, ET, IT, TT} <: ParticleMethod
    model::MT
    equation::ET
    integrator::IT
    timestep::TT

    function DiffEqIntegrator(model::MT, equation::ET, integrator::IT, timestep::TT) where {MT, ET, IT, TT}
        new{MT,ET,IT,TT}(model, equation, integrator, timestep)
    end
end

function run!(method::DiffEqIntegrator, sol)
    sol .= DifferentialEquations.solve(method.equation, method.integrator, dt = method.timestep)
end


function run(method::DiffEqIntegrator)
    sol = DifferentialEquations.solve(
        method.equation, 
        method.integrator, 
        dt = method.timestep, 
        adaptive = false, 
        progress = true,
        progress_steps = 1
        )
    return sol
end
