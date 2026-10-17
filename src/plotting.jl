# Plotting is a package extension. This file declares the three functions and documents them;
# `ext/VlasovMethodsMakieExt.jl` holds every method, so the package loads with no Makie
# installed, and each name is a function with no methods until `CairoMakie` is loaded.

@doc raw"""
    plot_distribution(dist) -> Figure
    plot_distribution(dist; xdomain, vdomain, nx = 100, nv = 100) -> Figure

A phase-space portrait of `dist` on one axis, `x` horizontally and `v` vertically.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_distribution(dist)
```

The one-argument method takes a [`GridDistribution`](@ref) and draws `dist.values` as a single
heatmap on the grid's own nodes, untransposed: the columns of the matrix are the `v` nodes, so
the figure and the matrix have the same orientation.

The method with the two keywords takes a `ParticleDistribution` of one position and one
velocity dimension and draws the weighted particle histogram on `nx × nv` cells of `xdomain`
and `vdomain`. Cell `(i, j)` holds

```math
F_{ij} = \frac{1}{h_x h_v} \sum_k w_k ,
```

summed over the particles that fall in it, with ``h_x = (x_b - x_a) / n_x`` and
``h_v = (v_b - v_a) / n_v``, so the histogram integrates to the total weight of the particles
it keeps. A cell is half-open, `[x_a + (i-1) h_x, x_a + i h_x)` and likewise in `v`.

`xdomain` is periodic and `x` wraps into it. `vdomain` is bounded, and a particle outside it —
`v < v_a` or `v >= v_b` — is dropped from the histogram rather than clamped into an end cell.
A particle with a non-finite `x`, which the wrap cannot place, is dropped as well.
A distribution with no particles throws an `ArgumentError`.

There are no other keywords, and the two domains are required. Use the returned `Figure` to set
a title, a size or an axis label.
"""
function plot_distribution end

@doc raw"""
    plot_particles(dist; xdomain, vdomain) -> Figure
    plot_particles(SS, ts, p; xdomain, vdomain) -> Figure

The particles of a `ParticleDistribution` of one position and one velocity dimension, as
a single scatter of `(x, v)`.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_particles(dist; xdomain = (0.0, 2π), vdomain = (-4.0, 4.0))
```

`xdomain` is periodic and `x` wraps into it. Every particle's velocity is plotted as it is, and
the axis limits are set to `xdomain` horizontally and `vdomain` vertically, so a particle at
`v` outside `vdomain` is left at that `v` and lies outside the axes. A distribution with no
particles throws an `ArgumentError`.

The second method takes a [`Snapshots`](@ref), and draws its saved state: the positions and the
velocities of the one component, `SS.X[1, :, ts, p]` against `SS.V[1, :, ts, p]`, as the scatter
above, with `x` wrapped into `xdomain`. It plots the state that a march saved, so a frame of
`ts = saved_indices(nₛ)[1]` shows the final state of the sweep. `ts` and `p` must index the saved
times and the parameters, else an `ArgumentError` names the range it left.

There are no other keywords, and the two domains are required. Use the returned `Figure` to set
a title, a size or an axis label.
"""
function plot_particles end

@doc raw"""
    plot_series(t, series; labels = nothing, xlabel = "t", ylabel = "") -> Figure

One line per element of `series` against the abscissae `t`, on one axis.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_series(ts, (x2, xmax); labels = ["x2", "xmax"], ylabel = "error")
fig = plot_series(ts, (W_hfm, W_rom); labels = ["HFM", "ROM"], ylabel = "field energy")
fig = plot_series(μ, (γ_hfm, γ_rom); labels = ["HFM", "ROM"], xlabel = "parameter", ylabel = "γ")
```

`series` is any collection of vectors of `length(t)` values: the entries of the `NamedTuple` of
[`trajectory_errors`](@ref) for one parameter, or a column of each entry for `(nₚ, nₛ, nparam)`
input; the `potential_error` of each parameter; the field energies `W` of a full-order and a
reduced march, the energy `H` of a conservative one, or its drift `ΔH/H`, against the saved times;
or the [`growth_rate`](@ref) values of a sweep against the equation parameter. `t` holds the
abscissa of each value, so one figure holds the quantities that share it.

`labels` names the lines. It is `nothing`, which draws no legend, or a collection with one label
per element of `series`, in the same order, which is drawn as a legend of the axis. `xlabel` and
`ylabel` label the axes. The axis is linear, so set a log scale on the axis of the returned
`Figure` to read a decaying error, or values that span orders of magnitude.

# Throws

- `ArgumentError` when `series` holds no line, or a line holds no value.
- `DimensionMismatch` when a series is not as long as `t`, and when `labels` does not hold one
  label per series.
"""
function plot_series end
