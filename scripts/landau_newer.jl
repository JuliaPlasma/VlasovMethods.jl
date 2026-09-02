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
tspan = (0.0, 1e-1)     # integration time interval
# tspan = (0.0, 1.0)     # integration time interval
domainv = (-1.75, +1.75) # size of interior domain (i.e. excluding outer "big" cell on either side)
boundary_layer = 6.0 # set this to 0 to not construct a grid with a large cell on either side
ν = 1.0  # collision frequency

tol = 5e-4 # Picard iteration tolerance for |x|_2
ftol = 5e-3 # Picard iteration tolerance for |f(x)|_∞  
niter = 5 # fixed number of Picard iterations
max_iters = 15 # max number of Picard iterations
n = 1 # number of quadrature nodes
β = 1.0 # damping parameter for the Picard iterations
m = 2 # depth for anderson acceleration
chunksize = 100

filename = "nknot=$(nknot)_tstep=1e-3_tend=$(tspan[end])_"

trange = (tspan[1] - tstep):tstep:tspan[2]

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 2, npart), NormalDistribution())
# dist = initialize!(ParticleDistribution(1, 2, npart), UniformDistribution())
# dist = initialize!(ParticleDistribution(1, 2, npart), DoubleMaxwellian(shift = 1.))

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 2, nknot, order, domainv, boundary_layer, :Periodic)

# second spline dist for diagnostics 
sdist2 = SplineDistribution(1, 2, nknot, order, domainv, boundary_layer, :Periodic)

# construct entropy 
entropy = CollisionEntropy(sdist)
s = zero(trange)

S = projection(dist.particles.v, dist, sdist)
s[2] = entropy()

# construct Landau operator
landau = Landau(dist, entropy; ν = ν)

# closure for vector field
landau_rhs!(v̇, v, params) = VlasovMethods.collisional_vectorfield!(v̇, v, params, landau)

params = (sdist2 = sdist2, n = n)
rhs = zero(dist.particles.v)

# J = VlasovMethods.compute_J_gl(sdist, 2)
# sdist2.coefficients .= J

v_iter = zeros(2, npart, length(trange), niter+1)
v̇_iter = zeros(2, npart, length(trange), niter+1)

v_full = zeros(2, npart, length(trange))
v_full[:, :, 2] .= dist.particles.v

rhs_full = zeros(2, npart, length(trange))
landau_rhs!(view(rhs_full, :, :, 2), dist.particles.v, params)

v_iter[:, :, 2, 1] .= dist.particles.v
v̇_iter[:, :, 2, 1] .= rhs_full[:, :, 2]

rhs_prev = zeros(2, npart, 2)

### Run actual code

@time for (i, t) in pairs(trange[3:end])
    println("i = ", i, ", t = ", t)
    # v_full[:,:,i+2] = VlasovMethods.Picard_iterate_Landau!(dist, sdist, tol, β, tstep, i+2, t, v_full[:,:,i], rhs_prev, sdist2, max_iters, m )
    sol = VlasovMethods.Picard_iterate_Landau_nls!(
        landau, tol, ftol, β, tstep, i+2, t, v_full[:, :, i + 1],
        v_full[:, :, i], rhs_prev, m, n, chunksize)
    v_full[:, :, i + 2] .= dist.particles.v
    rhs_full[:, :, i + 2] .= rhs_prev[:, :, 1]

    for q in eachindex(sol.v, sol.v̇)
        v_iter[:, :, i + 2, q] .= sol.v[q]
        v̇_iter[:, :, i + 2, q] .= sol.v̇[q]
    end

    # compute entropy
    projection(dist.particles.v, dist, sdist)
    s[2 + i] = entropy()
    println("   Entropy = ", s[2 + i])
    println()

    # SciMLBase.successful_retcode(sol.retcode) || break
end

# mom = [mapreduce(p -> p, +, v_full[:,:,n]) for n in axes(v_full,3)]./npart
# enr = [mapreduce(p -> p[1].^2, +, v_full[:,:,n]) for n in axes(v_full,3)]./npart

# using JLD2
# jldsave("data_" * filename;v_full)

# Uncomment below to produce animations of f, v̇_alpha, {v_alpha} in time

using GLMakie

stepsize = 0.1
x = (domainv[1] - boundary_layer):stepsize:(domainv[2] + boundary_layer)
y = (domainv[1] - boundary_layer):stepsize:(domainv[2] + boundary_layer)

smin = minimum(s[2:end])
smax = maximum(s[2:end])
Δs = abs(smax - smin)

# F = Figure()
# Ax = Axis3(F[1,1])

# plot initial vector field (first component)
# GLMakie.scatter!(Ax, v_full[1,:,2], v_full[2,:,2], rhs_full[1,:,2])

# plot initial projected f
# S = projection(v_full[:,:,2], dist, sdist)
# z = [S(xa, ya) for xa in x, ya in y]
# GLMakie.surface!(Ax, x, y, z)

# animation parameters
# anim_step = 1
# final_ind = length(trange)

# F = Figure(size = (1600, 900))
# record(F, "Landau_anim_" * filename * ".gif", 2:anim_step:final_ind) do i 
#     empty!(F)
#     ax11 = Axis(F[1,1], xlabel="v₁", ylabel="v₂")
#     plt = GLMakie.scatter!(ax11, v_full[1,:,i], v_full[2,:,i])
#     xlims!(ax11, -5, +5)
#     ylims!(ax11, -5, +5)

#     ax12 = Axis(F[1,2], xlabel="t", ylabel="S(t) - S(0)")#, title = "ΔS"
#     GLMakie.lines!(ax12, trange[2:i], s[2:i] .- s[2])
#     xlims!(ax12, tspan...)
#     ylims!(ax12, 0, Δs)

#     ax21 = Axis(F[2,1], xlabel="v₁", ylabel="v̇₁")#, title = "v̇₁"
#     GLMakie.scatter!(ax21, v_full[1,:,i], rhs_full[1,:,i])
#     xlims!(ax21, -5, +5)
#     ylims!(ax21, -9, +9)

#     ax22 = Axis(F[2,2], xlabel="v₂", ylabel="v̇₂")#, title = "v̇₂"
#     GLMakie.scatter!(ax22, v_full[2,:,i], rhs_full[2,:,i])
#     xlims!(ax22, -5, +5)
#     ylims!(ax22, -9, +9)

#     _S = projection(v_full[:,:,i], dist, sdist)
#     z = [_S([xa, ya]) for xa in x, ya in y]
#     axf = Axis3(F[1:2,3:4], xlabel="v₁", ylabel="v₂", zlabel="f", title = "(i = $(i-1))")
#     surf = GLMakie.surface!(axf,x,y,z)
# end

# F = Figure(size = (1600, 900))
# record(F, "Landau_vdot_" * filename * ".gif", 2:anim_step:final_ind) do i 
#     empty!(F)
#     ax1 = Axis3(F[1,1], xlabel="v₁", ylabel="v₂", title="v̇₁")
#     surf = GLMakie.scatter!(ax1, v_full[1,:,i], v_full[2,:,i], rhs_full[1,:,i])
#     xlims!(ax1, -3, +3)
#     ylims!(ax1, -3, +3)
#     zlims!(ax1, -9, +9)

#     ax2 = Axis3(F[1,2], xlabel="v₁", ylabel="v₂", title="v̇₂")
#     surf = GLMakie.scatter!(ax2, v_full[1,:,i], v_full[2,:,i], rhs_full[2,:,i])
#     xlims!(ax2, -3, +3)
#     ylims!(ax2, -3, +3)
#     zlims!(ax2, -9, +9)
# end

### Plot evolution of iterations of one time step ###

# F = Figure(size = (1600, 900))
# ind = 56
# record(F, "Landau_vdot_$(ind)_" * filename * ".gif", 1:6) do i 
#     empty!(F)
#     ax1 = Axis3(F[1,1], xlabel="v₁", ylabel="v₂", title="v̇₁")
#     surf = GLMakie.scatter!(ax1, v_full[1,:,ind], v_full[2,:,ind], v̇_iter[1,:,ind,i])
#     xlims!(ax1, -3, +3)
#     ylims!(ax1, -3, +3)
#     zlims!(ax1, -9, +9)

#     ax2 = Axis3(F[1,2], xlabel="v₁", ylabel="v₂", title="v̇₂")
#     surf = GLMakie.scatter!(ax2, v_full[1,:,ind], v_full[2,:,ind], v̇_iter[2,:,ind,i])
#     xlims!(ax2, -3, +3)
#     ylims!(ax2, -3, +3)
#     zlims!(ax2, -9, +9)
# end
