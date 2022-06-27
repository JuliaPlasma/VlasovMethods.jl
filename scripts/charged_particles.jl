
# import libraries
using GeometricIntegrators
using ProgressMeter

# import ParticleMethods package
using ParticleMethods

# import HDF5 package for storing solution
using HDF5

# parameters
Δt = 0.1    # time step size
nt = 100    # number of time steps
ns = 5      # number of samples

# output file
h5file = "charged_particles.hdf5"

# random initial conditions
z₀ = collect(hcat(rand(ns) .* 2 .- 1, rand(ns) .* 2 .- 1, zeros(ns), rand(ns) .- 0.5, rand(ns) .- 0.5, ones(ns))')

# electric field
E(x::Vector{DT}) where {DT} = DT[0, 0, cos(2π*x[3])]

# magnetic field
B(x::Vector{DT}) where {DT} = DT[0, 0, 1]

# vector field
function lorentz_force!(t, z, ż)
    for i in axes(z,2)
        x = z[1:3,i]
        v = z[4:6,i]

        e = E(x)
        b = B(x)

        ż[1,i] = v[1]
        ż[2,i] = v[2]
        ż[3,i] = v[3]

        ż[4,i] = e[1] + v[2] * b[3] - v[3] * b[2]
        ż[5,i] = e[2] + v[3] * b[1] - v[1] * b[3]
        ż[6,i] = e[3] + v[1] * b[2] - v[2] * b[1]
    end
end

# solution storage
function copy_to_hdf5(h5z, z, n)
    h5z[:,:,n+1] = z
end

# create an Equation instance
equ = ODE(lorentz_force!, z₀)

# create HDF5 file and copy initial conditions
h5  = h5open(h5file, "w")
h5z = create_dataset(h5, "z", eltype(z₀), ((6, ns, nt+1), (6, ns, -1)), chunk=(6,ns,1))
copy_to_hdf5(h5z, z₀, 0)

# create integrator
int = IntegratorExplicitEuler(equ, Δt)

# create atomic solution
asol = AtomicSolution(equ)

# copy initial conditons to atomic solution
set_initial_conditions!(asol, equ)

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

# set plot ranges
xlim = (-2, +2)
ylim = (-2, +2)
zlim = ( 0, 20)

# create empty 3d plot with ns empty series
plt = plot3d(
    ns,
    xlim = xlim,
    ylim = ylim,
    zlim = zlim,
    title = "Charged Particles",
    legend = false,
    marker = 5,
    size = (800, 600)
)

# plot solutions
for j in axes(z,3)
    for i in axes(z,2)
        push!(plt, i, z[1,i,j], z[2,i,j], z[3,i,j])
    end
end

# save figure to file
savefig("charged_particles.png")
