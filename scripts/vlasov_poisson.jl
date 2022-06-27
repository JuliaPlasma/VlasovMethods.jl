
# import libraries
using GeometricIntegrators
using ProgressMeter

# import ParticleMethods package
using ParticleMethods

# import HDF5 package for storing solution
using HDF5

# parameters
Δt = 0.1    # time step size
nt = 200    # number of time steps
np = 10000  # number of particles
nx = 16     # number of grid points

# output file
h5file = "vlasov_poisson.hdf5"

# random initial conditions
x₀ = randn(np)
v₀ = randn(np)

# shift x₀ to the interval [0,1]
xmax = ceil(maximum(abs.(x₀)))
x₀ .+= xmax
x₀ ./= 2*xmax

# concatenate x₀ and v₀
z₀ = collect(hcat(x₀, v₀)')

# vector field
function lorentz_force!(t, z, ż, p::PoissonSolver)
    Particles.solve!(p, z[1,:])
    for i in axes(ż, 2)
        ż[1,i] = z[2,i]
        ż[2,i] = eval_field(p, z[1,i])
    end
end

# splitting fields
function v_advection!(t, z, ż)
    for i in axes(ż, 2)
        ż[1,i] = z[2,i]
        ż[2,i] = 0
    end
end

function s_advection!(t, z, s, h)
    for i in axes(s, 2)
        s[1,i] = z[1,i] + h * z[2,i]
        s[2,i] = z[2,i]
    end
end

function v_lorentz_force!(t, z, ż, p::PoissonSolver)
    Particles.solve!(p, z[1,:])
    for i in axes(ż, 2)
        ż[1,i] = 0
        ż[2,i] = eval_field(p, z[1,i])
    end
end

function s_lorentz_force!(t, z, s, h, p::PoissonSolver)
    Particles.solve!(p, z[1,:])
    for i in axes(s, 2)
        s[1,i] = z[1,i]
        s[2,i] = z[2,i] + h * eval_field(p, z[1,i])
    end
end

# solution storage
function copy_to_hdf5(h5z, z, n)
    h5z[:,:,n+1] = z
end

# create Poisson solver
p = PoissonSolverFFT{eltype(z₀)}(nx)
Particles.solve!(p, x₀)

# create an ODE instance
# ode = ODE((t, z, ż) -> lorentz_force!(t, z, ż, p), z₀)

# create a splitting ODE instance
sode = SODE((v_advection!, (t, z, ż) -> v_lorentz_force!(t, z, ż, p)),
            (s_advection!, (t, z, s, h) -> s_lorentz_force!(t, z, s, h, p)),
            z₀)

# create HDF5 file and copy initial conditions
h5  = h5open(h5file, "w")
h5z = create_dataset(h5, "z", eltype(z₀), ((2, np, nt+1), (2, np, -1)), chunk=(2,np,1))
copy_to_hdf5(h5z, z₀, 0)

# create integrator
# int = IntegratorExplicitEuler(ode, Δt)
int = Integrator(sode, TableauStrang(), Δt)

# create atomic solution
# asol = AtomicSolution(ode)
asol = AtomicSolution(sode)

# copy initial conditons to atomic solution
# set_initial_conditions!(asol, ode)
set_initial_conditions!(asol, sode)

# initilize integrator
initialize!(int, asol)

# loop over time steps showing progress bar
try
    @showprogress 5 for n in 1:nt
        integrate_step!(int, asol)
        copy_to_hdf5(h5z, asol.q, n)
    end
finally
    # close HDF5 file
    close(h5)
end


# load Plots package
using Plots

# read array from HDF5 file
z = h5read(h5file, "z")

# compute plot ranges
vmax = ceil(maximum(abs.(v₀)))
xlim = (0, 1)
vlim = (-vmax, +vmax)

# plot initial condition
scatter(mod.(z[1,:,begin], 1), z[2,:,begin],
        marker = 3,
        xlim = xlim,
        ylim = vlim,
        title = "Vlasov-Poisson",
        legend = false,
        size = (800, 600)
)

# save figure to file
savefig("vlasov_poisson_z₀.png")


# create animation
anim = @animate for n in axes(z,3)
    scatter(mod.(z[1,:,n], 1), z[2,:,n],
        marker = 3,
        xlim = xlim,
        ylim = vlim,
        title = "Vlasov-Poisson",
        legend = false,
        size = (800, 600)
    )
end

# save animation to file
gif(anim, "vlasov_poisson_anim.gif", fps=10)
