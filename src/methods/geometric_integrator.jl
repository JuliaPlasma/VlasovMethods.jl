"""
    GeometricIntegrator(model, tspan, tstep[, method])

The method that advances `model` over the time span `tspan` with the step `tstep` as one
`GeometricEquations.ODEProblem`, solved by the `GeometricIntegrators` method `method`.

The fields `model`, `equation` and `integrator` hold the model, the problem and the integrator.
Each model defines its own constructor, with its own state and default `method`; the docstrings
of the constructors below list them. [`run!`](@ref) advances the integrator.
"""
struct GeometricIntegrator{MT, ET, IT} <: VlasovMethod
    model::MT
    equation::ET
    integrator::IT

    function GeometricIntegrator(model::MT, equation::ET,
            integrator::IT) where {MT, ET <: GeometricEquations.GeometricProblem, IT}
        new{MT, ET, IT}(model, equation, integrator)
    end
end

"""
    run!(int::GeometricIntegrator, h5file)

Advance the method `int` over its time span, write every state into the HDF5 file `h5file`, copy
the final state into the distribution of the model, and return that distribution.

The file holds the dataset `z` of the size `(n, ntime)`, one column per state, and the dataset
`t` of the `ntime` times of `t₀:Δt:t₁`, the initial time included. The file is created, or
overwritten, and closed on return, also when a step throws. A progress bar is shown for a run of
more than five seconds.
"""
function run!(int::GeometricIntegrator, h5file)
    # initial conditions
    prob = int.equation
    t₀ = GeometricEquations.initialtime(prob)
    z₀ = prob.ics.q

    # the length of the state
    np = size(z₀, 1)

    # The time grid of a `Solution` of the same problem, as in `run!(::SplittingMethod, h5file)`:
    # one point per state, the initial one included, which can be one point fewer than
    # `ntime(prob) + 1` where `ntime` counts a partial last step that the range drops.
    times = t₀:GeometricEquations.timestep(prob):GeometricEquations.finaltime(prob)

    # create HDF5 file and copy initial conditions
    h5 = h5open(h5file, "w")
    h5z = create_dataset(
        h5, "z", eltype(z₀), ((np, length(times)), (np, -1)), chunk = (np, 1))
    h5t = create_dataset(
        h5, "t", eltype(t₀), ((length(times),), (-1,)), chunk = (1,))
    h5z[:, 1] = z₀
    h5t[1] = t₀

    # One solution step is advanced through the run and written out after each step, as in
    # `run!(::SplittingMethod, h5file)`.
    solstep = GeometricIntegrators.GeometricIntegratorsBase.solutionstep(
        int.integrator, GeometricEquations.initialstate(prob))

    # loop over time steps showing progress bar
    try
        @showprogress 5 for n in 1:(length(times) - 1)
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
