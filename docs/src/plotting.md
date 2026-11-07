```@meta
CurrentModule = VlasovMethods
```

# Plotting

The plotting functions are defined in a package extension: load `CairoMakie` to make them
available. Each returns a Makie `Figure`.

## Distributions and series

```@docs
plot_distribution
plot_particles
plot_series
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
