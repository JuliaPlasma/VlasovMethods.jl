```@meta
CurrentModule = VlasovMethods
```

# Getting started

This page runs two complete simulations: a collisionless Vlasov–Poisson run on particles, and a
collision run of the conservative Lenard–Bernstein operator. Both follow the same three steps of
[Running a model](methods.md): build a distribution, build a model, and advance it with a method.

## A Vlasov–Poisson run

The first run advances 10 000 particles in their own electrostatic field on the periodic interval
``[0, 1)``.

```@example start-vp
using HDF5
using PoissonSolvers
using Random
using VlasovMethods

Random.seed!(42)

# the distribution: 10 000 particles of one position and one velocity component
dist = initialize!(ParticleDistribution(1, 1, 10_000), NormalDistribution())

# the model: the particles and a quadratic periodic spline potential (order 3) on 16 cells
potential = Potential(PeriodicBasisSpline((0.0, 1.0), 3, 16))
model = VlasovPoisson(dist, potential)

# the method: a Strang splitting over t ∈ [0, 4] with the step 0.1
method = SplittingMethod(model, (0.0, 4.0), 0.1)

h5file = tempname() * ".hdf5"
run!(method, h5file);
```

[`run!`](@ref) writes the trajectory into `h5file` and returns the distribution with the final
state. The dataset `z` holds the positions and the velocities of all particles at all 41 times:

```@example start-vp
z = h5read(h5file, "z")
size(z)
```

The plotting functions draw the phase space of the trajectory:

```@example start-vp
using CairoMakie
CairoMakie.activate!(type = "png")

plot_phase_space(z, [1, 21, 41]; xdomain = (0.0, 1.0),
    titles = ["t = 0", "t = 2", "t = 4"])
```

## A collision run

The second run relaxes a shifted normal distribution of 200 particles under the conservative
Lenard–Bernstein operator. The operator acts on the velocities only. It needs a spline
distribution beside the particles, onto which it projects the particles at every evaluation.

```@example start-clb
using HDF5
using Random
using VlasovMethods

Random.seed!(1)

# the particles: 200 velocities of the mean 2 and the variance 1, of equal weight
pdist = ParticleDistribution(1, 1, 200)
pdist.particles.v .= 2 .+ randn(1, 200)
pdist.particles.w .= 1 / 200

# the spline: a cubic basis on 21 breakpoints of [-10, 10]
sdist = SplineDistribution(1, 1, 21, 4, (-10.0, 10.0))

model = ConservativeLenardBernstein(pdist, CollisionEntropy(sdist); ν = 1.0)
method = VlasovMethods.GeometricIntegrator(model, (0.0, 0.5), 0.02)

h5file = tempname() * ".hdf5"
run!(method, h5file)

v = h5read(h5file, "z")
t = h5read(h5file, "t")
size(v)
```

The default time integrator of the conservative operator is the implicit midpoint rule. The
operator conserves the momentum and the energy of the particles, and the implicit midpoint rule
conserves both quadratic invariants, so they stay constant to rounding:

```@example start-clb
w = pdist.particles.w[1, :]
momentum = [sum(w .* v[:, s]) for s in axes(v, 2)]
energy = [sum(w .* v[:, s] .^ 2) for s in axes(v, 2)]
(maximum(abs, momentum .- momentum[1]) / momentum[1],
 maximum(abs, energy .- energy[1]) / energy[1])
```

The entropy ``S_h = \int f_s \log |f_s| \, dv`` of the projected distribution decreases:

```@example start-clb
entropy = CollisionEntropy(sdist)
S = map(axes(v, 2)) do s
    projection(v[:, s], pdist, sdist)
    entropy()
end

using CairoMakie
CairoMakie.activate!(type = "png")
plot_entropy(t, S)
```

## Where to go next

  - [Distribution functions](distributions.md) and [Projections](projections.md) explain the
    particle, spline and grid representations that these runs use.
  - [Vlasov–Poisson](vlasov_poisson.md) adds the grid discretisation and the energy of the run.
  - [Lenard–Bernstein collisions](lenard_bernstein.md) and [Landau collisions](landau.md) cover
    all collision operators and their time integrators.
  - [Reproducing the papers](../papers/index.md) runs the experiments of the papers.
