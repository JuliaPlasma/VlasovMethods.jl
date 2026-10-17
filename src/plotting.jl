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
    plot_errors(t, series) -> Figure

The error norms of a sweep against its saved times `t`: one line per element of `series`.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_errors(t, (x2, xmax))
```

`series` is any collection of series of `length(t)` values — the entries of the `NamedTuple` of
[`trajectory_errors`](@ref), the `potential_error` of each parameter, or the two together. Each
element is one line, and the axis is labelled `t` and `error`; the axis is linear, so set a log
scale on the returned `Figure`'s axis to read a decaying error.

# Throws

- `ArgumentError` when `series` holds no line.
- `DimensionMismatch` when a series is not as long as the times.
"""
function plot_errors end

@doc raw"""
    plot_energy(t, series) -> Figure

The field energy of a sweep against its saved times `t`: one line per element of `series`.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_energy(t, (W_hfm, W_rom))
```

`series` is any collection of series of `length(t)` values: the field energies `W` of a full-order
and a reduced march, the energy `H` of a conservative one, or its drift `ΔH/H`. Each element is
one line, and the axis is labelled `t` and `energy`, so a caller plots the two quantities of one
item on one figure. The axis is linear, and the caller sets a log scale on the returned `Figure`
where the values span orders of magnitude.

# Throws

- `ArgumentError` when `series` holds no line.
- `DimensionMismatch` when a series is not as long as the times.
"""
function plot_energy end

@doc raw"""
    plot_growth_rate(t, series) -> Figure

The growth rates of a sweep against the equation parameter `t`: one line per element of `series`.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_growth_rate(μ, (γ_rom, γ_hfm))
```

`t` holds the parameter of each sample of the sweep and `series` any collection of
`length(t)`-long series of [`growth_rate`](@ref) values. Each element is one line, and the axis is
labelled `parameter` and `γ`, so the full-order rates of a sweep against the reduced ones show the
model's error across the parameter range.

# Throws

- `ArgumentError` when `series` holds no line.
- `DimensionMismatch` when a series is not as long as the parameter values.
"""
function plot_growth_rate end
