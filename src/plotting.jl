# Plotting is a package extension. This file declares the two functions and documents them;
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
A distribution with no particles throws an `ArgumentError`.

There are no other keywords, and the two domains are required. Use the returned `Figure` to set
a title, a size or an axis label.
"""
function plot_distribution end

@doc raw"""
    plot_particles(dist; xdomain, vdomain) -> Figure

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

There are no other keywords, and the two domains are required. Use the returned `Figure` to set
a title, a size or an axis label.
"""
function plot_particles end
