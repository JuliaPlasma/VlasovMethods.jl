# import stuff
using VlasovMethods
using SciMLBase

# using JLD2

# using BenchmarkTools
# using ProfileView
# using Cthulhu
using Profile

# parameters
# npart = 200   # number of particles
npart = 2000  # number of particles
nknot = 4     # number of grid points in one direction (in the interior)
order = 2      # spline order
tstep = 1e-3    # time step size
tspan = (0.0, 1.0)     # integration time interval
domainv = (-1.75, 1.75) # size of interior domain (i.e. excluding outer "big" cell on either side)
length_big_cell = 6.0 # set this to 0 to not construct a grid with a large cell on either side
ν = 1.0  # collision frequency

filename = "tstep=1e-4_bc=4_particlebounds_tests_normal"

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 2, npart), NormalDistribution())
# dist = initialize!(ParticleDistribution(1, 2, npart), UniformDistribution())
# dist = initialize!(ParticleDistribution(1, 2, npart), DoubleMaxwellian(shift = 1.))

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 2, nknot, order, domainv, length_big_cell, :Periodic)

# second spline dist for diagnostics 
sdist2 = SplineDistribution(1, 2, nknot, order, domainv, length_big_cell, :Periodic)

# construct entropy 
entropy = CollisionEntropy(sdist)

trange = (tspan[1] - tstep):tstep:tspan[2]

S = projection(dist.particles.v, dist, sdist)

# construct Landau operator
landau = Landau(dist, entropy; ν = ν)

# closure for vector field
landau_rhs!(v̇, v, params) = VlasovMethods.collisional_vectorfield!(v̇, v, params, landau)

params = (sdist2 = sdist2, n = 2)
rhs = zero(dist.particles.v)

# J = VlasovMethods.compute_J_gl(sdist, 2)
# sdist2.coefficients .= J

v_full = zeros(2, npart, length(trange))
v_full[:, :, 2] .= dist.particles.v

rhs_full = zeros(2, npart, length(trange))
landau_rhs!(rhs_full[:, :, 2], dist.particles.v, params)

rhs_prev = zeros(2, npart, 2)

tol = 5e-4 # Picard iteration tolerance for |x|_2
ftol = 5e-3 # Picard iteration tolerance for |f(x)|_∞  
max_iters = 15 # max number of Picard iterations
β = 1.0 #damping parameter for the Picard iterations
m = 2 # depth for anderson acceleration
n = 2
chunksize = 100

### Run profiler

# i = 3
# t = trange[i]

# VlasovMethods.Picard_iterate_Landau_nls!(landau, tol, ftol, β, tstep, i+2, t, v_full[:,:,i+1], v_full[:,:,i], rhs_prev, m, n, chunksize)

# Profile.clear()
# Profile.clear_malloc_data()

# Profile.Allocs.@profile VlasovMethods.Picard_iterate_Landau_nls!(landau, tol, ftol, β, tstep, i+2, t, v_full[:,:,i+1], v_full[:,:,i], rhs_prev, m, n, chunksize)

### Run actual code

@time for (i, t) in pairs(trange[3:end])
    println("i=", i, " t =", t)
    # v_full[:,:,i+2] = VlasovMethods.Picard_iterate_Landau!(dist, sdist, tol, β, tstep, i+2, t, v_full[:,:,i], rhs_prev, sdist2, max_iters, m )
    sol = VlasovMethods.Picard_iterate_Landau_nls!(
        landau, tol, ftol, β, tstep, i+2, t, v_full[:, :, i + 1],
        v_full[:, :, i], rhs_prev, m, n, chunksize)
    v_full[:, :, i + 2] .= dist.particles.v
    rhs_full[:, :, i + 2] .= rhs_prev[:, :, 1]
    if !SciMLBase.successful_retcode(sol.retcode)
        break
    end
end

# # # mom = [mapreduce(p -> p, +, v_full[:,:,n]) for n in axes(v_full,3)]./npart
# # # enr = [mapreduce(p -> p[1].^2, +, v_full[:,:,n]) for n in axes(v_full,3)]./npart

# using JLD2
# jldsave("data_" * filename;v_full)

# # # # # # Uncomment below to produce animations of f, v̇_alpha, {v_alpha} in time

# using GLMakie

# stepsize = 0.1
# x = domainv[1]-length_big_cell:stepsize:domainv[2]+length_big_cell
# y = domainv[1]-length_big_cell:stepsize:domainv[2]+length_big_cell

# F = Figure()
# Ax = Axis3(F[1,1])

# # plot initial vector field (first component)
# GLMakie.scatter!(Ax, v_full[1,:,2], v_full[2,:,2], rhs_full[1,:,2])

# plot initial projected f
# S = projection(v_full[:,:,2], dist, sdist)
# z = [S(xa, ya) for xa in x, ya in y]
# GLMakie.surface!(Ax, x, y, z)

# animation parameters
# anim_step = 1
# final_ind = 24

# F = Figure()
# record(F, "f_v_anim_" * filename * ".gif", 2:anim_step:final_ind) do i 
#     empty!(F)
#     Ax = Axis(F[1,1], title = "i = $i")
#     plt = GLMakie.scatter!(Ax, v_full[1,:,i], v_full[2,:,i])
# end 

# record(F, "f_surface_anim_" * filename * ".gif", 2:anim_step:final_ind) do i 
#     _S = projection(v_full[:,:,i], dist, sdist)
#     z = [_S([xa, ya]) for xa in x, ya in y]
#     empty!(F)
#     Ax = Axis3(F[1,1], title = "i = $i")
#     surf = GLMakie.surface!(Ax,x,y,z)
# end 

# F = Figure()
# record(F, "f_v1dot_anim_" * filename * ".gif", 2:final_ind) do i 
#     empty!(F)
#     Ax = Axis3(F[1,1], title = "i = $i")
#     GLMakie.scatter!(Ax, v_full[1,:,i], v_full[2,:,i], rhs_full[1,:,i])
# end 

# F = Figure()

# record(F, "f_v2dot_anim_" * filename * ".gif", 2:final_ind) do i 
#     empty!(F)
#     Ax = Axis3(F[1,1], title = "i = $i")
#     GLMakie.scatter!(Ax, v_full[1,:,i], v_full[2,:,i], rhs_full[2,:,i])
# end 
