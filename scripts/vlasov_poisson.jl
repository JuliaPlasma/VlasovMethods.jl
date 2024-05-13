# import libraries
using PoissonSolvers
using VlasovMethods

# parameters
npart = 10000  # number of particles
nknot = 16     # number of grid points
order = 3      # spline order
tstep = 0.1    # time step size
tspan = (0.0, 20.0)    # integration time interval
domain = (0.0, 1.0)

# output file
h5file = "vlasov_poisson.hdf5"

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 1, npart), NormalDistribution())
# dist = initialize!(ParticleDistribution(1, 1, npart), BumpOnTail())

# create electrostatic potential
potential = Potential(PeriodicBasisBSplineKit(domain, order, nknot))

# create Vlasov-Poisson model
model = VlasovPoisson(dist, potential)

# create integrator
integrator = SplittingMethod(model, tspan, tstep)

# integrate
run!(integrator, h5file)


# load HDF5 and Plots packages
using HDF5
using Plots

# read array from HDF5 file
z = h5read(h5file, "z")

# compute plot ranges
vmax = ceil(maximum(abs.(z[2,:,begin])))
xlim = (0, 1)
vlim = (-vmax, +vmax)

# plot initial condition
scatter(mod.(z[1,:,begin], 1), z[2,:,begin],
        marker = 3,
        xlim = xlim,
        ylim = vlim,
        title = "Vlasov-Poisson",
        legend = false,
        size = (800, 600)
)

# save figure to file
savefig("vlasov_poisson_z₀.png")


# create animation
anim = @animate for n in axes(z,3)
    scatter(mod.(z[1,:,n], 1), z[2,:,n],
        marker = 3,
        xlim = xlim,
        ylim = vlim,
        title = "Vlasov-Poisson",
        legend = false,
        size = (800, 600)
    )
end

# save animation to file
gif(anim, "vlasov_poisson.gif", fps=10)
