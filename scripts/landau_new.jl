# import stuff
using BSplineKit
using VlasovMethods
using GeometricIntegrators.Extrapolators
using JLD2

# parameters
npart = 10000  # number of particles
nknot = 5     # number of grid points in one direction
order = 4      # spline order
tstep = 2e-6    # time step size
tspan = (0.0, 2e-2)     # integration time interval
domainv = (-3., 3.)
# domainv = (-7., 7.)

filename = "tstep=1e-5_T=5e-2_nls_beta=0.5_test"

# create and initialize particle distribution function
# dist = initialize!(ParticleDistribution(1, 2, npart), NormalDistribution())
dist = initialize!(ParticleDistribution(1, 2, npart), UniformDistribution())
# dist = initialize!(ParticleDistribution(1, 2, npart), DoubleMaxwellian(shift = 2.0))

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 2, nknot, order, domainv, :Dirichlet)
# second spline dist for diagnostics 
sdist2 = SplineDistribution(1, 2, nknot, order, domainv, :Dirichlet)
# entropy = CollisionEntropy(sdist)

trange = (tspan[1] - tstep):tstep:tspan[2]

v_full = zeros(2, npart, length(trange))
v_full[:,:,2] .= dist.particles.v

rhs_prev = zeros(2, npart, 2)

tol = 5e-4 # Picard iteration tolerance for |x|_2
ftol = 5e-3 # Picard iteration tolerance for |f(x)|_∞  
max_iters = 15 # max number of Picard iterations
β = 0.5 #damping parameter for the Picard iterations
m = 2 # depth for anderson acceleration
@time for (i,t) in pairs(trange[3:end])
    println("i=",i, " t =",t)
    # v_full[:,:,i+2] = VlasovMethods.Picard_iterate_Landau!(dist, sdist, tol, β, tstep, i+2, t, v_full[:,:,i], rhs_prev, sdist2, max_iters, m )
    v_full[:,:,i+2] .= VlasovMethods.Picard_iterate_Landau_nls!(dist, sdist, tol, ftol, β, tstep, i+2, t, v_full[:,:,i+1], v_full[:,:,i], rhs_prev, sdist2, max_iters, m )
end

# using JLD2
# jldsave("data_" * filename;v_full)

# # Uncomment below to produce animation of f in time

# using GLMakie

# stepsize = 0.1
# x = domainv[1]:stepsize:domainv[2]
# y = domainv[1]:stepsize:domainv[2]

# F = Figure()
# Ax = Axis3(F[1,1])


# record(F, "f_surface_anim_" * filename * ".gif", 2:length(trange)) do i 
#     # record(F, filename * ".gif", 2:length(trange)) do i 
#     _S = projection(v_full[:,:,i], dist, sdist)
#     z = [_S([xa, ya]) for xa in x, ya in y]
#     empty!(F)
#     Ax = Axis3(F[1,1], title = "i = $i")
#     surf = surface!(Ax,x,y,z)
# end 