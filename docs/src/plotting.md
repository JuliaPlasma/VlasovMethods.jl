```@meta
CurrentModule = VlasovMethods
```

# Plotting

Every plotting function is declared in `src/plotting.jl` and has its methods in the package
extension `VlasovMethodsMakieExt`, which loading `CairoMakie` activates. Each returns a `Figure`
and writes no file: a driver in `scripts/` loads `CairoMakie`, calls these functions, and saves
each `Figure` under a name of its own. The functions of the collision figures are on the page of
[The figures of the collision manuscripts](collisions.md).

```@docs
plot_distribution
plot_particles
plot_series
plot_groups
plot_damping
plot_phase_space
```
