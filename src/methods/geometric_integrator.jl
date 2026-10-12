struct GeometricIntegrator{MT, ET, IT} <: VlasovMethod
    model::MT
    equation::ET
    integrator::IT

    function GeometricIntegrator(model::MT, equation::ET,
            integrator::IT) where {MT, ET <: GeometricEquations.GeometricProblem, IT}
        new{MT, ET, IT}(model, equation, integrator)
    end
end

function run!(int::GeometricIntegrator, h5file)
    # initial conditions
    prob = int.equation
    t₀ = GeometricEquations.initialtime(prob)
    z₀ = prob.ics.q

    # the length of the state
    np = size(z₀, 1)

    # create HDF5 file and copy initial conditions
    h5 = h5open(h5file, "w")
    h5z = create_dataset(
        h5, "z", eltype(z₀), ((np, ntime(prob)+1), (np, -1)), chunk = (np, 1))
    h5t = create_dataset(
        h5, "t", eltype(t₀), ((ntime(prob)+1,), (-1,)), chunk = (1,))
    h5z[:, 1] = z₀
    h5t[1] = t₀

    # One solution step is advanced through the run and written out after each step, as in
    # `run!(::SplittingMethod, h5file)`.
    solstep = GeometricIntegrators.GeometricIntegratorsBase.solutionstep(
        int.integrator, GeometricEquations.initialstate(prob))
    times = t₀:GeometricEquations.timestep(prob):GeometricEquations.finaltime(prob)

    # loop over time steps showing progress bar
    try
        @showprogress 5 for n in 1:ntime(prob)
            GeometricIntegrators.reset!(solstep, times[n + 1])
            GeometricIntegrators.integrate!(solstep, int.integrator)
            h5z[:, n + 1] = solstep.q
            h5t[n + 1] = times[n + 1]
        end
    finally
        # close HDF5 file
        close(h5)
    end

    return _write_back!(int.model, solstep.q)
end

# The particle models step the first row of `dist.particles.v`; a grid model keeps the whole state in `distribution.values` and has its own method, in
# `src/gridbased/vlasov_poisson_1d1v.jl`. Both return the distribution.
function _write_back!(model, q)
    model.dist.particles.v[1, :] .= q

    return model.dist
end
