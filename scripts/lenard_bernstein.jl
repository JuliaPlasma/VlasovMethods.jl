# import stuff
using BSplineKit
using VlasovMethods
using Plots

# params
# parameters
npart = 1000  # number of particles
nknot = 31     # number of grid points
order = 4      # spline order
tstep = 0.5    # time step size
tspan = (0.0, 1e3)    # integration time interval
domainv = (-6., 6.)

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 1, npart), NormalDistribution())

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 1, nknot, order, domainv, :Dirichlet)
entropy = CollisionEntropy(sdist)

# create LenardBernstein model
model = LenardBernstein(dist, entropy)

# create integrator
integrator = DiffEqIntegrator(model, tspan, tstep)

# run integrator 
# sol = ... #TODO need to construct this somehow...
# run!(integrator, sol)

println("Running integrator")
sol = VlasovMethods.run(integrator)

# create animation
println("producing animation")
step = 1
xgrid = -8.25:0.125:+8.25
vgrid = -8:0.01:+8
params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)


f = projection(sol[:,1], dist, sdist)
v = LB_rhs(collect(vgrid), params, f)

anim = @animate for i in 1:step:length(sol)
    println("i=",i)
    
    # compute quantities for plotting
    f = projection(sol[:,i], dist, sdist)
    df = Derivative(1) * f
    v = LB_rhs(collect(vgrid), params, f)

    plot(xlims = [-8, +8], ylims = [-0.5, +0.5], size=(1200,800))

    histogram!(sol[:,i], bins=xgrid, normalize=:probability, label=string(tstep*(i-1)))
    plot!(vgrid, v, color=:red, label="v")
    plot!(xgrid, f.(xgrid), label="f")
    plot!(xgrid, df.(xgrid), label="df")
end
println("saving animation")
gif(anim,"lenard_bernstein.gif",fps=2)

