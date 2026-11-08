```@meta
CurrentModule = VlasovMethods
```

# Plotting

Every plotting function is declared in `src/plotting.jl` and has its methods in the package
extension `VlasovMethodsMakieExt`, which loading `CairoMakie` activates. Each returns a `Figure`
and writes no file: a driver in `scripts/` loads `CairoMakie`, calls these functions, and saves
each `Figure` under a name of its own.

## Distributions and series

```@docs
plot_distribution
plot_particles
plot_series
plot_groups
plot_damping
plot_phase_space
plot_velocity_projections
```

## The collision figures

```@docs
plot_collision_distribution
plot_drift
plot_entropy
plot_error
plot_convergence
plot_cost
```
