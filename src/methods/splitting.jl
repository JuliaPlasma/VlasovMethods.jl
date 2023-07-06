
struct SplittingMethod{MT, ET, IT} <: ParticleMethod
    model::MT
    equation::ET
    integrator::IT

    function SplittingMethod(model::MT, equation::ET, integrator::IT) where {MT <: VlasovModel, ET <: GeometricEquations.GeometricProblem, IT}
        new{MT,ET,IT}(model, equation, integrator)
    end
end





# solution storage
function copy_to_hdf5(h5z, z, n)
    h5z[:,:,n+1] = z
end



function run!(method::SplittingMethod, h5file)
    # initial conditions
    z₀ = method.equation.ics.q

    # dimensions and number of particles
    nd = size(z₀,1)
    np = size(z₀,2)

    # create HDF5 file and copy initial conditions
    h5  = h5open(h5file, "w")
    h5z = create_dataset(h5, "z", eltype(z₀), ((nd, np, ntime(method.equation)+1), (nd, np, -1)), chunk=(nd,np,1))
    copy_to_hdf5(h5z, z₀, 0)

    Integrators.initialize!(method.integrator)

    # loop over time steps showing progress bar
    try
        @showprogress 5 for n in 1:ntime(method.equation)
            Integrators.integrate!(method.integrator)
            copy_to_hdf5(h5z, method.integrator.solstep.q, n)
        end
    finally
        # close HDF5 file
        close(h5)
    end

    copy!(method.model.distribution.particles.z, method.integrator.solstep.q)

    return method.model.distribution
end
