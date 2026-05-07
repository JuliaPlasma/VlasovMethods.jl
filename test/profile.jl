# run with julia --project --track-allocation=user test/profile.jl

using BSplineKit
using Profile
using VlasovMethods

# output file
h5file = "profile.hdf5"

# parameters
npart = 1000   # number of particles
nknot = 41     # number of grid points
order = 4      # spline order
tstep = 1e-3   # time step size
tspan = (0.0, 1e-2)    # integration time interval
domainv = (-10., 10.)

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 1, npart), NormalDistribution())

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 1, nknot, order, domainv, :Dirichlet)
entropy = CollisionEntropy(sdist)

# create LenardBernstein model
# model = LenardBernstein(dist, entropy)
model = ConservativeLenardBernstein(dist, entropy)

# create integrator
integrator = GeometricIntegrator(model, tspan, tstep)

# clear profile cache
Profile.clear()
Profile.clear_malloc_data()

# run integrator 
Profile.Allocs.@profile VlasovMethods.run!(integrator, h5file)
