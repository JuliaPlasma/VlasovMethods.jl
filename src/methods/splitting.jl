
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

function run!(method::SplittingMethod, h5file)
    # initial conditions
    z₀ = method.equation.ics.q

    # dimensions and number of particles
    nd = size(z₀, 1)
    np = size(z₀, 2)

    # The time grid of a `Solution` of the same problem, the range `t₀:Δt:t₁`. Its length is the
    # number of states a run holds, the initial one included. It is one point shorter than
    # `ntime(equation) + 1` wherever the span is not an exact multiple of the step, because
    # `ntime` rounds the step count up (`div(t₁ - t₀, Δt, RoundUp)`).
    prob = method.equation
    t₀ = GeometricEquations.initialtime(prob)
    t₁ = GeometricEquations.finaltime(prob)
    times = t₀:GeometricEquations.timestep(prob):t₁

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
            GeometricIntegrators.reset!(solstep, times[n + 1])
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
