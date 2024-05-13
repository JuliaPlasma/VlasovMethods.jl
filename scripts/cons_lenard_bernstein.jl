# using Logging: global_logger
# using TerminalLoggers: TerminalLogger
# global_logger(TerminalLogger())

# import stuff
using BSplineKit
using VlasovMethods
using QuadGK

# output file
# run_name = "fixed_coefficient_N=1000_double_maxwellian_T=2_tstep=5e-3"
# run_name = "shifted_normal_N=1000_T=1000_new"
run_name = "uniform_N=200_T=1000_new"
h5file = "lenard_bernstein_" * run_name * ".hdf5"

# params
# parameters
npart = 1000  # number of particles
nknot = 41     # number of grid points
order = 4      # spline order
tstep = 8e-4    # time step size
tspan = (0.0, 1)    # integration time interval
domainv = (-10., 10.)
domainv_p = (-2., 2.)

# create and initialize particle distribution function
# dist = initialize!(ParticleDistribution(1, 1, npart), NormalDistribution())
# dist = initialize!(ParticleDistribution(1, 1, npart), ShiftedNormalV())
dist = initialize!(ParticleDistribution(1, 1, npart), UniformDistribution((0.,1.0), domainv_p))
# dist = initialize!(ParticleDistribution(1, 1, npart), DoubleMaxwellian(domainv, 2.0))
# dist = initialize!(ParticleDistribution(1, 1, npart), Bump(a = 2., b = 1.))
# dist = initialize!(ParticleDistribution(1, 1, npart), ShiftedUniformDistribution())

# create spline distribution function and entropy 
sdist = SplineDistribution(1, 1, nknot, order, domainv, :Dirichlet)
entropy = CollisionEntropy(sdist)

# create LenardBernstein model
model = ConservativeLenardBernstein(dist, entropy)
# model = LenardBernstein(dist, entropy)

# # create integrator
# integrator = GeometricIntegrator(model, tspan, tstep)
# # integrator = DiffEqIntegrator(model, tspan, tstep)

# # run integrator 
# println("Running integrator")
# VlasovMethods.run!(integrator, h5file)
# # sol = VlasovMethods.run(integrator)


# load HDF5 and Plots packages
using HDF5
using Plots
using LaTeXStrings

# read array from HDF5 file
z = h5read(h5file, "z")
t = h5read(h5file, "t")

t ./= npart

# z = z[:, 1:201]
# t = t[1:201]
# z = sol
# t = sol.t

mom = [mapreduce(p -> p[1], +, z[:,n]) for n in axes(z,2)]./npart
enr = [mapreduce(p -> p[1].^2, +, z[:,n]) for n in axes(z,2)]./npart
# third_mom = [mapreduce(p -> (p[1].^3)./npart + 2*mom[n]^3 - 3*mom[n]*enr[n], +, z[:,n]) for n in axes(z,2)]
# third_mom = [mapreduce(p -> (p[1] - mom[n]).^3, +, z[:,n]) for n in axes(z,2)]./npart
# # fourth_mom = [mapreduce(p -> (p[1] - mom[1]).^4, +, z[:,n]) for n in axes(z,2)]./npart
# fourth_cum = ([mapreduce(p -> (((p[1] - mom[n]).^4)./npart - 3*(enr[n] - mom[n]^2)^2), +, z[:,n]) for n in axes(z,2)].^2)
# fifth_cum = ([mapreduce(p -> (((p[1] - mom[n]).^5)./npart - 10*(enr[n] - mom[n]^2)*third_mom[n]), +, z[:,n]) for n in axes(z,2)].^2)

# κ₃ = copy(mom)
# κ₃_s = copy(mom)
# κ₄_s = copy(mom)
# κ₄ = copy(mom)
# κ₅ = copy(mom)
# for i in axes(z,2)
#     # κ₃[i] = sum(z[:,i].^3)./npart + 2*mom[i]^3 - 3*mom[i]*enr[i] 

#     f = projection(z[:,i], dist, sdist)
#     κ₃[i] = sum((z[:,i] .- mom[i]).^3)./npart  
#     κ₃_s[i], err = quadgk(v -> (v .- mom[i]).^3 * f(v), domainv[1], domainv[2])
#     κ₄[i] = sum((z[:,i] .- mom[i]).^4)./npart - 3*(enr[i] - mom[i]^2)^2
#     κ₄_s[i], err = quadgk(v -> ((v .- mom[i]).^4 .- 3*(enr[i] - mom[i]^2)^2) * f(v), domainv[1], domainv[2])
#     κ₅[i] = sum((z[:,i] .- mom[i]).^5)./npart - 10*(enr[i] - mom[i]^2)*κ₃[i]
# end


scalefontsizes()
# # plot(t, abs.(κ₃./κ₃[1]), yaxis=:log10, label=L"\kappa_3 / \kappa_3(t=0)")
# plot(t, abs.(κ₄./κ₄[1]), label=L"\kappa_4_s / \kappa_4(t=0)", yaxis=:log10, xlabel = "t", legendfontsize = 12)
# plot!(t, abs.(κ₅./κ₅[1]), label=L"\kappa_5 / \kappa_5(t=0)", yaxis=:log10)
# # plot!(t, exp.(-3 .* t), linestyle=:dash, label=L"y = \exp{(-3t)}")
# plot!(t, exp.(-4 .* t), linestyle=:dash, label=L"y = \exp{(-4t)}")
# plot!(t, exp.(-5 .* t), linestyle=:dash, label=L"y = \exp{(-5t)}")
# scalefontsizes(1.25)
# savefig("cumulant_decay_" * run_name * ".pdf")


xgrid = -8.:0.25:+8.
vgrid = -8:0.01:+8
params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)

# compute plot ranges
vmax = ceil(maximum(abs.(z[:,begin])))
xlim = (0, 1)
vlim = (-vmax, +vmax)
nplot = 100
knot_list = knots(sdist.basis)
@show step_size = size(z,2) ÷ nplot
act_len = length(1:step_size:size(z,2))
ent = zeros(act_len)
dSdt = zeros(act_len, size(z,1))
v̇ = zeros(act_len, size(z,1))
# create animations
anim = @animate for (i, n) in zip(1:act_len, 1:step_size:size(z,2))
    println(i)
    # anim = @animate for n in 1:step:size(z,2)
    println("i= $n, mom = $((mom[n] - mom[1])/mom[1]), enr=$((enr[n] - enr[1])/enr[1]), min(v) = $(minimum(z[:,n])), max(v) = $(maximum(z[:,n]))")

    # compute quantities for plotting
    f = projection(z[:,n], dist, sdist)
    ent[i], err = quadgk(x -> 0.5 * f(x) * log(f(x)^2), domainv[1], domainv[2])
    # df = Derivative(1) * f
    # # v = LB_rhs(collect(vgrid), params, f)
    # v = CLB_rhs(collect(vgrid), params, f)
    # v̇[i, :] = CLB_rhs(z[:,n], params, f)
    # dSdt[i, :] .= df.(z[:,n]) .* v̇[i, :]


    # plot(xlabel = "v", xlims = [-8, +8], ylims = [-0.5, 0.5], size=(1200,800),legendfontsize = 12)

    # histogram!(z[:,n], bins=xgrid, normalize=:pdf, label="particle distribution (t ="*string(t[n])*")")
    # # scatter!(knot_list, zeros(length(knot_list)))
    # plot!(vgrid, v, color=:red, label=L"\partial_t v")
    # plot!(vgrid, f.(vgrid), label="spline-projected distribution")
    # plot!(xgrid, df.(xgrid), label="df")
end
# gif(anim, "lenard_bernstein_anim_" * run_name * ".gif", fps=5) 

# entropy evolution
scalefontsizes()
p = plot(t[1:step_size:size(z,2)], ent./abs(ent[1]), xlabel="t", ylabel=L"S/|S(t=0)|", legend=false)
scalefontsizes(1.5)
savefig(p, "final_entropy_log_" * run_name * ".pdf")


# p = plot(t[1:step:step * nplot], sum(dSdt, dims=2), label="dS/dt")
# savefig(p, "dSdt_" * run_name * ".pdf")


# # # save initial condition to file
function plot_distributions(t_ind, z,  dist::ParticleDistribution, sdist::SplineDistribution, xgrid, vgrid)
    f = projection(z[:,t_ind], dist, sdist)
    df = Derivative(1) * f
    # v = CLB_rhs(collect(vgrid), params, f)
    v̇ = CLB_rhs(z[:, t_ind], params, f)
    scalefontsizes()
    p = plot(xlabel = "v", xlims = [-7., +7], ylims = [0, 0.4], size = [1200,800], legendfontsize=14)
    histogram!(p, z[:,t_ind], bins=xgrid, normalize=:pdf, label="particle distribution")
    # scatter!(z[:, t_ind], v̇, color=:red, label=L"\partial_t v")
    plot!(p, vgrid, f.(vgrid), lw = 3, label="spline-projected distribution")
    # plot!(p, vgrid, df.(vgrid), lw = 3, label="df")
    
    return p
end

# tsteps = Int.([1, floor(length(t)/2), length(t)])
# nbins = Int.(floor.([nknot./2, nknot, 2*nknot]))
# domain_length = domainv[2] - domainv[1]
# folder = "histogram_scans/"
# for ts in tsteps
#     for nb in nbins
#         @show spacing = domain_length./nb
#         xgrid = -8.:spacing:+8.
#         p = plot_distributions(ts, z, dist, sdist, xgrid, vgrid)
#         savefig(p, folder * "shifted_normal/" * "N=" * string(npart) * "_M=" * string(nknot) * "_nbins=" * string(nb) * "_tstep=" * string(ts) * ".pdf")
#     end
# end

# # plot initial condition and final result
# p1 = plot_distributions(1, z, dist, sdist, xgrid, vgrid)
# scalefontsizes(1.5)
# savefig(p1, "initial_distribution_" * run_name * ".pdf")

# # savefig(p1, "shifted_maxwellian_" * "N=" * string(npart) * "_M=" * string(nknot) * "_tstep=1" * ".pdf")

# # n_end = length(t)
# # p2 = plot_distributions(n_end, z, dist, sdist, xgrid, vgrid)
# # savefig(p2, "shifted_maxwellian_" *  "N=" * string(npart) * "_M=" * string(nknot) * "_tstep=" * string(n_end) * ".pdf")


# n_end = length(t)
# p2 = plot_distributions(n_end, z, dist, sdist, xgrid, vgrid)
# savefig(p2, "final_distribution_" * run_name * ".pdf")


# plot(p1, p2, layout=(1,2))
# savefig("initial_final_distribution_" * run_name * ".pdf")

# plot(t, (mom .- mom[1])/mom[1], label = "momentum", xlabel="t", ylabel = "relative error")
# plot!(t, (enr .- enr[1])/enr[1], label = "energy")


# f_eval = zeros(npart, length(t))
# for i in 1:length(t)
#     local f = projection(z[:,i], dist, sdist)
#     f_eval[:,i] .= f.(z[:,i])
# end
# mom_s = [mapreduce(p -> p[1], +, f_eval[:,n]) for n in axes(z,2)]./npart
# enr_s = [mapreduce(p -> p[1].^2, +, f_eval[:,n]) for n in axes(z,2)]./npart
# third_s_mom = [mapreduce(p -> (p[1] - mom_s[1]).^3, +, f_eval[:,n]) for n in axes(z,2)]./npart
# fourth_s_cum = ([mapreduce(p -> ((p[1] - mom_s[1]).^4 - 3*enr_s[1]^2), +, f_eval[:,n]) for n in axes(z,2)].^2)./npart

# # plot(third_mom, label = "particle third moment")
# # plot(third_s_mom, yaxis=:log10, label = "spline third moment")
# # savefig("third_moment" * run_name * ".pdf")

# # plot(fourth_cum, label = "particle fourth cumulant")
# # plot!(fourth_s_cum, label = "spline fourth cumulant")
# # savefig("fourth_cumulant" * run_name * ".pdf"