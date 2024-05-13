# import stuff
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using BSplineKit
using VlasovMethods
using Plots

# params
# parameters
npart = 1000  # number of particles
nknot = 41     # number of grid points
order = 4      # spline order
tstep = 0.1    # time step size
tspan = (0.0, 1e1)    # integration time interval
domainv = (-10., 10.)

# create and initialize particle distribution function
dist = initialize!(ParticleDistribution(1, 1, npart), NormalDistribution())

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 1, nknot, order, domainv, :Dirichlet)
entropy = CollisionEntropy(sdist)

# create LenardBernstein model
model = ConservativeLenardBernstein(dist, entropy)
# model = LenardBernstein(dist, entropy)

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
v = VlasovMethods.LB_rhs(collect(vgrid), params, f)
mom = [mapreduce(p -> p[1], +, sol[:,i]) for i in 1:step:length(sol)]
enr = [mapreduce(p -> p[1].^2, +, sol[:,i]) for i in 1:step:length(sol)]

anim = @animate for i in 1:step:length(sol)
    println("i= $i, mom = $((mom[i] - mom[1])/mom[1]), enr=$((enr[i] - enr[1])/enr[1])")

    # compute quantities for plotting
    f = projection(sol[:,i], dist, sdist)
    df = Derivative(1) * f
    v = VlasovMethods.LB_rhs(collect(vgrid), params, f)
    # v = VlasovMethods.CLB_rhs(collect(vgrid), params, f)

    plot(xlims = [-8, +8], ylims = [-0.5, +0.5], size=(1200,800))

    histogram!(sol[:,i], bins=xgrid, normalize=:probability, label=string(sol.t[i]))
    plot!(vgrid, v, color=:red, label="v")
    plot!(xgrid, f.(xgrid), label="f")
    plot!(xgrid, df.(xgrid), label="df")
end
println("saving animation")
gif(anim,"lenard_bernstein.gif",fps=2)
