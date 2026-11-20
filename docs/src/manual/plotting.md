```@meta
CurrentModule = VlasovMethods
```

# Plotting

The plotting functions are declared in `src/plotting.jl` and have their methods in the package
extension `VlasovMethodsMakieExt`, which loading `CairoMakie` activates:

```julia
using VlasovMethods, CairoMakie
```

Each function returns a `Figure` and writes no file. A driver in `scripts/` calls these functions
and saves each `Figure` under a name of its own. The figures of the papers are drawn by them, and
the examples of this manual use them too.

## Distributions and series

| function | draws |
|:--|:--|
| [`plot_distribution`](@ref) | a grid distribution as a heatmap, or the weighted histogram of 1D1V particles |
| [`plot_particles`](@ref) | the particles of a distribution or of a saved state of [`Snapshots`](@ref) as a scatter |
| [`plot_series`](@ref) | several series against time on one axis |
| [`plot_groups`](@ref) | groups of series, one colour per group |
| [`plot_damping`](@ref) | the peaks of ``\log |E_1|`` of damping runs with their fitted rates |
| [`plot_phase_space`](@ref) | the particles of a trajectory at several time slices, one panel each |
| [`plot_velocity_projections`](@ref) | the right-hand side of a collision model and the spline projections of the particles |

## The collision figures

| function | draws |
|:--|:--|
| [`plot_collision_distribution`](@ref) | the initial and the final distribution, with the Maxwellian of the initial moments |
| [`plot_drift`](@ref) | the relative drift of the conserved quantities, on a logarithmic axis |
| [`plot_entropy`](@ref) | the entropy ``S_h`` against time |
| [`plot_error`](@ref) | the error of a run against a closed form, and the closed form |
| [`plot_convergence`](@ref) | an error against the cell width, with the fitted power law |
| [`plot_cost`](@ref) | the error of each run against its wall-clock cost |

## Conventions of the logarithmic axes

A quantity that is conserved exactly has a relative drift of zero, which a logarithmic axis cannot
show. [`plot_drift`](@ref), [`plot_error`](@ref) and [`plot_cost`](@ref) draw such a line one
decade below the smallest positive value on the axis, and drop the points where a time or a value
is not finite. An axis on which every value is the same spans a decade on each side of it. Each
function throws an `ArgumentError` where nothing finite is left to draw.
