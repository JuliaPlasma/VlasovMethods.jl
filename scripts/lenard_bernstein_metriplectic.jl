# using Logging: global_logger
# using TerminalLoggers: TerminalLogger
# global_logger(TerminalLogger())

# import stuff
using ProgressMeter
using BSplineKit
using VlasovMethods
using QuadGK
using SciMLBase

# output file
# run_name = "fixed_coefficient_N=1000_double_maxwellian_T=2_tstep=5e-3"
run_name = "downstairs_factor_shifted_uniform_N=1000_dt=1e-3_T=2/"
path = "results/" * "lenard_bernstein_metriplectic/" * run_name
if isdir(path) #check if results directory exists
    println("Run directory already exists, files will be overwritten")
else
    println("Creating run directory")
    mkpath(path)
    # mkdir(path)
end

# h5file = "lenard_bernstein_metriplectic_" * run_name * ".hdf5"
h5file = path * "run_data.hdf5"

# params
# parameters
npart = 2000  # number of particles
nknot = 41     # number of grid points
order = 4      # spline order
tstep = 1e-3    # time step size
tspan = (0.0, 2.0)    # integration time interval
length_big_cell = 0.0
domainv = (-10.0, 10.0)
domainv_p = (-2.0, 2.0)

# create and initialize particle distribution function
# dist = initialize!(ParticleDistribution(1, 1, npart), NormalDistribution())
# dist = initialize!(ParticleDistribution(1, 1, npart), ShiftedNormalV())
# dist = initialize!(ParticleDistribution(1, 1, npart), SumMaxwellian(var_two = 16.))
# dist = initialize!(ParticleDistribution(1, 1, npart), DoubleMaxwellian(shift = 2.))
# dist = initialize!(ParticleDistribution(1, 1, npart), Bump(a = 2., b = 1.))
dist = initialize!(ParticleDistribution(1, 1, npart), ShiftedUniformDistribution())

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 1, nknot, order, domainv, length_big_cell, :nothing)
entropy = CollisionEntropy(sdist)

# create LenardBernstein model
model = MetriplecticLenardBernstein(dist, entropy)

v = zeros(npart)
v .= transpose(dist.particles.v)
# v̇ = zero(v)
# max_iters = 100000
# tol = 1e-14
# ftol = 1e-14

dv = zeros(npart)

# trange = tspan[1]:tstep:0.04652
trange = tspan[1]:tstep:tspan[2]

sol = zeros(npart, length(trange))
sol[:, 1] .= transpose(model.dist.particles.v)
dv_history = zeros(npart, 2)

β = 0.5 #damping parameter for the Picard iterations
m = 3 # depth for anderson acceleration

abstol = 1e-15
reltol = 1e-50 # make this 1e-50 and check results 
# check initial and final residuals

# for i in 1:1
@showprogress for (i, t) in pairs(trange[1:(end - 1)])
    @show i, t
    if i > 1
        sol_object = Picard_iterate_over_particles(
            dv, sol[:, i], sol[:, i - 1], dv_history,
            i, t, tstep, m, β, abstol, reltol, model)
        sol[:, i + 1] .= sol_object.u
        # @show sol_object.resid

        if !SciMLBase.successful_retcode(sol_object)
            @show sol_object.retcode
            return
        end
    else
        sol_object = Picard_iterate_over_particles(
            dv, sol[:, i], sol[:, i], dv_history, i, t, tstep, m, β, abstol, reltol, model)
        sol[:, i + 1] .= sol_object.u
    end

    # @show j
end

mom = [mapreduce(p -> p[1], +, sol[:, n]) for n in axes(sol, 2)] ./ npart
enr = [mapreduce(p -> p[1] .^ 2, +, sol[:, n]) for n in axes(sol, 2)] ./ npart

using Plots

xgrid = -10.0:0.5:+10.0
vgrid = -10:0.01:+10
nplot = 100
@show step = length(trange) ÷ nplot
ent = zeros(length(trange))

# create animations
# anim = @animate for n in 1:length(trange)
# anim = @animate for (i, n) in zip(1:nplot, 1:step:length(trange))
#     println(n)

#     f = projection(sol[:,n], dist, sdist)

#     ent[n], err = VlasovMethods.compute_entropy(x -> f(x), model)

#     plot(xlabel = "v", xlims = [-8, +8], ylims = [-0.5, 0.5], size=(1200,800),legendfontsize = 12, title = "T = " * string(trange[n]))

#     histogram!(sol[:,n], bins=xgrid, normalize=:pdf, label="particle distribution")

#     plot!(vgrid, f.(vgrid), label="spline-projected distribution")
# end

# # # # # save initial condition to file
# gif(anim, "metriplectic_lenard_bernstein_anim_" * run_name * ".gif", fps=5)

plot(trange, (mom .- mom[1])/mom[1], label = "momentum",
    xlabel = "t", ylabel = "relative error")
plot!(trange, (enr .- enr[1])/enr[1], label = "energy")
savefig(path * "conservation.pdf")

# # p = plot(t[1:step:step * nplot], ent./ent[1], xlabel="t", ylabel="entropy (normalised)", legend=false)
# # savefig(p, "entropy_log_" * run_name * ".pdf")

# # # p = plot(t[1:step:step * nplot], sum(dSdt, dims=2), label="dS/dt")
# # # savefig(p, "dSdt_" * run_name * ".pdf")

function plot_distributions(
        t_ind, z, dist::ParticleDistribution, sdist::SplineDistribution, xgrid, vgrid)
    f = projection(z[:, t_ind], dist, sdist)
    # df = Derivative(1) * f
    scalefontsizes()
    p = plot(xlabel = "v", xlims = [-9.0, +9.0], ylims = [0.0, +0.6],
        size = [1200, 800], legendfontsize = 14)
    histogram!(
        p, z[:, t_ind], bins = xgrid, normalize = :pdf, label = "particle distribution")
    # scatter!(z[:, t_ind], v̇, color=:red, label=L"\partial_t v")
    plot!(p, vgrid, f.(vgrid), lw = 3, label = "spline-projected distribution")
    return p
end

# plot initial condition and final result
p1 = plot_distributions(1, sol, dist, sdist, xgrid, vgrid)
scalefontsizes(1.5)
savefig(p1, path * "initial_distribution.pdf")
# savefig(p1, "initial_distribution_" * run_name * ".pdf")

n_end = length(trange)
p2 = plot_distributions(n_end, sol, dist, sdist, xgrid, vgrid)
savefig(p2, path * "final_distribution.pdf")
# savefig(p2, "final_distribution_" * run_name * ".pdf")

plot(p1, p2, layout = (1, 2))
savefig(path * "initial_final_distribution.pdf")
# savefig("initial_final_distribution_" * run_name * ".pdf")
