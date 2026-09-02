using VlasovMethods
using Profile

# parameters
npart = 2000  # number of particles
nknot = 4     # number of grid points in one direction (in the interior)
order = 2      # spline order
tstep = 1e-3    # time step size
tspan = (0.0, 1e-3)     # integration time interval
domainv = (-1.75, +1.75) # size of interior domain (i.e. excluding outer "big" cell on either side)
length_big_cell = 6.0 # set this to 0 to not construct a grid with a large cell on either side
ν = 1.0  # collision frequency

trange = tspan[begin]:tstep:tspan[end]

# create and initialize particle distribution function
pdist = initialize!(ParticleDistribution(1, 2, npart), NormalDistribution())

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 2, nknot, order, domainv, length_big_cell, :Periodic)

# construct entropy 
entropy = CollisionEntropy(sdist)

# construct Landau operator
landau = Landau(dist, entropy; ν = ν)

# closure for vector field
landau_rhs!(v̇, v, params) = VlasovMethods.collisional_vectorfield!(v̇, v, params, landau)

# initial projection
S = projection(dist.particles.v, dist, sdist)

params = (sdist2 = sdist2, n = 2)
rhs = zero(dist.particles.v)

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

i = 3
t = trange[i]

VlasovMethods.Picard_iterate_Landau_nls!(
    landau, tol, ftol, β, tstep, i+2, t, v_full[:, :, i + 1],
    v_full[:, :, i], rhs_prev, m, n, chunksize)

Profile.clear()
Profile.clear_malloc_data()

Profile.Allocs.@profile VlasovMethods.Picard_iterate_Landau_nls!(
    landau, tol, ftol, β, tstep, i+2, t, v_full[:, :, i + 1],
    v_full[:, :, i], rhs_prev, m, n, chunksize)
