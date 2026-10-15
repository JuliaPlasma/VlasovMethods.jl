# The Lenard–Bernstein driver: a normal-distribution particle cloud on a cubic spline basis,
# advanced by the model's own `GeometricIntegrator`, with the right-hand side at the initial
# condition and the distribution at both ends of the run drawn as figures.
#
#   julia --project=scripts scripts/lenard_bernstein.jl [--model lb|clb|rclb] [--npart 1000] \
#       [--nknot 41] [--order 4] [--tstep 0.01] [--final-time 1.0] \
#       [--runs-dir dir] [--results-dir dir]
#
# The HDF5 trajectory goes to `--runs-dir`, the figure to `--results-dir`.
#
# `--model` picks the collision operator: `lb` is `LenardBernstein`, `clb` the conservative
# variant that solves its two coefficients from the particle momentum and energy, and `rclb` the
# rescaled conservative variant. `clb` is the default.
#
# The basis is unconstrained at the ends (`Free()`), the one condition under which the moments of
# the spline representation equal the particles'; on a `Dirichlet()` basis the conservative
# models warn that they do not.
#
# **The defaults are dt = 0.01 over t ∈ [0, 1].** A dt = 0.1 over t ∈ [0, 10] fails for every one
# of the three models, because the particles drift out of the velocity domain and the projection
# then throws on the ones that left.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using CairoMakie
using HDF5
using VlasovMethods

const MODELS = Dict("lb" => LenardBernstein, "clb" => ConservativeLenardBernstein,
    "rclb" => RescaledConservativeLenardBernstein)

const OPT = parse_arguments(ARGS; model = "clb", npart = 1000, nknot = 41, order = 4,
    tstep = 1e-2, final_time = 1.0)
const MODEL = choose(MODELS, "--model", OPT.model)
const TSPAN = (0.0, OPT.final_time)
const DOMAIN = (-10.0, 10.0)

dist = initialize!(ParticleDistribution(1, 1, OPT.npart), NormalDistribution())
sdist = SplineDistribution(1, 1, OPT.nknot, OPT.order, DOMAIN, 0, :nothing)
model = MODEL(dist, CollisionEntropy(sdist))
integrator = VlasovMethods.GeometricIntegrator(model, TSPAN, OPT.tstep)

# The right-hand side at the initial condition, through the vector field and the parameters the
# integrator steps, evaluated before the run writes the final state back into `dist.particles.v`.
prob = integrator.equation
v₀ = dist.particles.v[1, :]
v̇₀ = zero(v₀)
prob.equation.v(v̇₀, TSPAN[1], v₀, prob.parameters)

h5file = runs_path("lenard_bernstein_" * OPT.model * ".hdf5")
run!(integrator, h5file)

@printf("lenard_bernstein: model %s, %d particles, %d knots, order %d, dt = %g, t ∈ %s\n",
    OPT.model, OPT.npart, OPT.nknot, OPT.order, OPT.tstep, TSPAN)

# `z` is `(npart, ntime + 1)`: the velocity vector the ODE steps, one column per time in `t`.
z = h5read(h5file, "z")
t = h5read(h5file, "t")
@printf("  wrote %s (%d time slices)\n", h5file, length(t))

xgrid = range(DOMAIN[1], DOMAIN[2]; length = 201)
edges = range(DOMAIN[1], DOMAIN[2]; length = 81)

fig = Figure(size = (1600, 450))

ax1 = Axis(fig[1, 1], xlabel = "v", ylabel = "v̇", title = "right-hand side at t = 0")
scatter!(ax1, v₀, v̇₀, markersize = 2)

for (panel, slice) in enumerate((1, length(t)))
    f = projection(z[:, slice], dist, sdist)
    df = derivative(f)
    ax = Axis(fig[1, panel + 1], xlabel = "v", ylabel = "f",
        title = @sprintf("t = %.2f", t[slice]))
    hist!(ax, z[:, slice], bins = edges, normalization = :pdf, color = (:steelblue, 0.5),
        label = "particles")
    lines!(ax, xgrid, f.(xgrid), color = :firebrick, label = "spline projection")
    lines!(ax, xgrid, df.(xgrid), color = :black, linestyle = :dash, label = "derivative")
    panel == 1 && axislegend(ax; position = :lt)
end

figure = results_path("lenard_bernstein_" * OPT.model * ".png")
save(figure, fig)
@printf("  wrote %s\n", figure)
