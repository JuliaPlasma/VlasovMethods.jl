struct GeometricIntegrator{MT, ET, IT} <: ParticleMethod
    model::MT
    equation::ET
    integrator::IT

    function GeometricIntegrator(model::MT, equation::ET, integrator::IT) where {MT, ET <: GeometricEquations.GeometricProblem, IT}
        new{MT,ET,IT}(model, equation, integrator)
    end
end


function run!(method::GeometricIntegrator, h5file)
    # initial conditions
    t₀ = method.equation.tspan[begin]
    z₀ = method.equation.ics.q

    # dimensions and number of particles
    np = size(z₀,1)

    # create HDF5 file and copy initial conditions
    h5  = h5open(h5file, "w")
    h5z = create_dataset(h5, "z", eltype(z₀), ((np, ntime(method.equation)+1), (np, -1)), chunk=(np,1))
    h5t = create_dataset(h5, "t", eltype(t₀), ((ntime(method.equation)+1,), (-1,)), chunk=(1,))
    h5z[:,1] = z₀
    h5t[1] = t₀

    Integrators.initialize!(method.integrator)

    # loop over time steps showing progress bar
    try
        @showprogress 5 for n in 1:ntime(method.equation)
            Integrators.integrate!(method.integrator)
            h5z[:,n+1] = method.integrator.solstep.q
            h5t[n+1] = method.integrator.solstep.t
        end
    finally
        # close HDF5 file
        close(h5)
    end

    method.model.dist.particles.v[1,:] .= method.integrator.solstep.q

    return method.model.dist
end
