struct GeometricIntegrator{MT, ET, IT} <: ParticleMethod
    model::MT
    equation::ET
    integrator::IT

    function GeometricIntegrator(model::MT, equation::ET, integrator::IT) where {MT, ET <: GeometricEquations.GeometricProblem, IT}
        new{MT,ET,IT}(model, equation, integrator)
    end
end


function run!(int::GeometricIntegrator, h5file)
    # initial conditions
    t₀ = int.equation.tspan[begin]
    z₀ = int.equation.ics.q

    # dimensions and number of particles
    np = size(z₀,1)

    # create HDF5 file and copy initial conditions
    h5  = h5open(h5file, "w")
    h5z = create_dataset(h5, "z", eltype(z₀), ((np, ntime(int.equation)+1), (np, -1)), chunk=(np,1))
    h5t = create_dataset(h5, "t", eltype(t₀), ((ntime(int.equation)+1,), (-1,)), chunk=(1,))
    h5z[:,1] = z₀
    h5t[1] = t₀

    solstep = Integrators.SolutionStep(Integrators.problem(int.integrator), Integrators.method(int.integrator))

    Integrators.initialize!(solstep, Integrators.problem(int.integrator))

    # loop over time steps showing progress bar
    try
        @showprogress 5 for n in 1:ntime(int.equation)
            Integrators.integrate!(solstep, int.integrator)
            h5z[:,n+1] = solstep.q
            h5t[n+1] = solstep.t
        end
    finally
        # close HDF5 file
        close(h5)
    end

    int.model.dist.particles.v[1,:] .= solstep.q

    return int.model.dist
end
