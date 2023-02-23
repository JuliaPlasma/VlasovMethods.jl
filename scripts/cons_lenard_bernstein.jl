# import stuff
using BSplineKit
using VlasovMethods

# output file
h5file = "lenard_bernstein.hdf5"

# params
# parameters
npart = 1000  # number of particles
nknot = 41     # number of grid points
order = 4      # spline order
tstep = 1e-2    # time step size
tspan = (0.0, 5e2)    # integration time interval
domainv = (-10., 10.)

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 1, npart), DoubleMaxwellian(domainv, 2.0))
# dist = initialize!(ParticleDistribution(1, 1, npart), NormalDistribution())

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 1, nknot, order, domainv, :Dirichlet)
entropy = CollisionEntropy(sdist)

# create LenardBernstein model
model = ConservativeLenardBernstein(dist, entropy)
# model = LenardBernstein(dist, entropy)

# create integrator
integrator = GeometricIntegrator(model, tspan, tstep)
# integrator = DiffEqIntegrator(model, tspan, tstep)


# run integrator 
# sol = ... #TODO need to construct this somehow...
# run!(integrator, sol)

println("Running integrator")
VlasovMethods.run!(integrator, h5file)

# load HDF5 and Plots packages
using HDF5
using Plots

# read array from HDF5 file
z = h5read(h5file, "z")
t = h5read(h5file, "t")

mom = [mapreduce(p -> p[1], +, z[:,n]) for n in axes(z,2)]
enr = [mapreduce(p -> p[1].^2, +, z[:,n]) for n in axes(z,2)]

xgrid = -8.25:0.125:+8.25
vgrid = -8:0.01:+8
params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)

# compute plot ranges
vmax = ceil(maximum(abs.(z[:,begin])))
xlim = (0, 1)
vlim = (-vmax, +vmax)
nplot = 100
step = size(z,2) ÷ nplot
# create animation
anim = @animate for n in 1:step:size(z,2)
    println("i= $n, mom = $((mom[n] - mom[1])/mom[1]), enr=$((enr[n] - enr[1])/enr[1]), min(v) = $(minimum(z[:,n])), max(v) = $(maximum(z[:,n]))")

    # compute quantities for plotting
    f = projection(z[:,n], dist, sdist)
    df = Derivative(1) * f
    # v = LB_rhs(collect(vgrid), params, f)
    v = CLB_rhs(collect(vgrid), params, f)

    plot(xlims = [-8, +8], ylims = [-0.5, +0.5], size=(1200,800))

    histogram!(z[:,n], bins=xgrid, normalize=:probability, label=string(t[n]))
    plot!(vgrid, v, color=:red, label="v")
    plot!(xgrid, f.(xgrid), label="f")
    plot!(xgrid, df.(xgrid), label="df")

end

# save animation to file
gif(anim, "lenard_bernstein_anim.gif", fps=10)
