
"""
    SplittingMethod(model, tspan, tstep)

The splitting method of the particle model `model` over the time span `tspan` with the step
`tstep`: the `GeometricEquations.SODEProblem` of the parts of the model's vector field and
the step of each part, and the `GeometricIntegrators` integrator of their Strang composition. The
state is a copy of `model.distribution.particles.z`.

A 1D1V [`VlasovPoisson`](@ref) model of a `ParticleDistribution` splits into the advection
``\\dot{x} = v`` and the acceleration ``\\dot{v} = E(x)``, each solved exactly. A
[`VlasovWithExternalFields`](@ref) model splits into the drift ``\\dot{x} = v``, solved exactly,
and the Lorentz force, whose step is a Boris push, an approximation of its exact solution.

The fields `model`, `equation` and `integrator` hold the three parts. [`run!`](@ref) advances
the integrator, and `GeometricIntegrators.integrate(method.integrator)` returns the solution in
memory.
"""
struct SplittingMethod{MT, ET, IT} <: ParticleMethod
    model::MT
    equation::ET
    integrator::IT

    function SplittingMethod(model::MT,
            equation::ET,
            integrator::IT) where {
            MT <: AbstractVlasovModel, ET <: GeometricEquations.GeometricProblem, IT}
        new{MT, ET, IT}(model, equation, integrator)
    end
end

# The splitting method of a model from GeometricIntegrators, over a copy of the model's particle
# state: one solution step pushes all particles. This keeps the implementation simple, but does
# not suit parallelisation. `fields` and `solutions` are the vector fields and their exact
# solutions, one pair per component, and `splitting` is the composition.
function _splitting_method(model, fields::Tuple, solutions::Tuple, params, tspan::Tuple,
        tstep::Real, splitting)
    equ = GeometricEquations.SODEProblem(fields, solutions, tspan, tstep,
        copy(model.distribution.particles.z); parameters = params)
    int = GeometricIntegrators.GeometricIntegrator(equ, splitting)
    return SplittingMethod(model, equ, int)
end

# solution storage
function copy_to_hdf5(h5z, z, n)
    h5z[:, :, n + 1] = z
end

"""
    run!(method::SplittingMethod, h5file)

Advance the splitting method `method` over its time span, write every state into the HDF5 file
`h5file`, copy the final state into `method.model.distribution.particles.z`, and return that
distribution.

The file holds the dataset `z` of the size `(nd, npart, ntime)`: `z[:, :, s]` is the state of
the particles at the `s`-th time of the solution that `GeometricIntegrators.integrate` returns for
the same problem, the initial state included. The file is created,
or overwritten, and closed on return, also when a step throws.
"""
function run!(method::SplittingMethod, h5file)
    # initial conditions
    z₀ = method.equation.ics.q

    # dimensions and number of particles
    nd = size(z₀, 1)
    np = size(z₀, 2)

    # The time grid of a `Solution` of the same problem, built by the constructor that
    # `GeometricSolutions` builds it with, so a run holds one state per time of
    # `GeometricIntegrators.integrate`, the initial one included. `times[n]` is the time after `n`
    # steps. GeometricSolutions up to 0.7.0 builds the range `t₀:Δt:t₁`, and 0.7.1 and later
    # `ntime(equation) + 1` times, so the two differ by one time on a span whose last step ends a
    # few ulps past `t₁`.
    prob = method.equation
    t₀ = GeometricEquations.initialtime(prob)
    t₁ = GeometricEquations.finaltime(prob)
    times = GeometricIntegrators.TimeSeries(t₀, t₁, GeometricEquations.timestep(prob))

    # create HDF5 file and copy initial conditions. A chunk dimension of zero is invalid, and a
    # zero-particle run has one: the dataset itself may be empty, its chunk may not.
    h5 = h5open(h5file, "w")
    h5z = create_dataset(
        h5, "z", eltype(z₀), ((nd, np, length(times)), (nd, np, -1)),
        chunk = (nd, max(np, 1), 1))
    copy_to_hdf5(h5z, z₀, 0)

    # One solution step is advanced through the run and written out after each step, so only one
    # state is held. The step times are those of a `Solution` of the same problem.
    solstep = GeometricIntegrators.GeometricIntegratorsBase.solutionstep(
        method.integrator, GeometricEquations.initialstate(prob))
    try
        for n in 1:(length(times) - 1)
            GeometricIntegrators.reset!(solstep, times[n])
            GeometricIntegrators.integrate!(solstep, method.integrator)
            copy_to_hdf5(h5z, solstep.q, n)
        end
    finally
        # close HDF5 file
        close(h5)
    end

    copy!(method.model.distribution.particles.z, solstep.q)

    return method.model.distribution
end
