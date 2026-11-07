# Plotting is a package extension. This file declares the plotting functions and documents them;
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

@doc raw"""
    plot_collision_distribution(f, g; xdomain, vdomain, nx = 100, nv = 100) -> Figure

The initial distribution `f` and the final distribution `g` of a collision run, in two panels of
one `Figure`, titled `initial` and `final`, with the Maxwellian of `f` drawn on each as a contour.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_collision_distribution(f, g; xdomain = (0.0, 1.0), vdomain = (-10.0, 10.0))
```

`f` and `g` are `ParticleDistribution`s of one position and one velocity dimension. Each panel
is the weighted particle histogram of [`plot_distribution`](@ref) on `nx × nv` cells of
`xdomain` and `vdomain`, with the same cells, the same normalisation and the same particles
dropped, so each integrates to the total weight of the particles it keeps.

The contour is the Maxwellian that carries the mass, the momentum and the energy of `f`,
evaluated on the cell centres,

```math
M(x, v) = \frac{n}{\sqrt{2 \pi T}} \exp\left( - \frac{(v - u)^2}{2 T} \right) , \qquad
n = \frac{W}{x_b - x_a} , \quad
u = \frac{1}{W} \sum_k w_k v_k , \quad
T = \frac{1}{W} \sum_k w_k (v_k - u)^2 ,
```

with ``W = \sum_k w_k``, the sums running over the particles of `f` whose position and velocity
are finite. A collision operator that conserves mass, momentum and energy relaxes a spatially
uniform `f` towards this Maxwellian, so the final panel shows how close `g` has come to it.

# Throws

- `ArgumentError` when `nx` or `nv` is not positive, when a domain does not run from its lower
  to its upper end, when `f` or `g` holds no particle, and when `f` has no positive mass or no
  positive temperature, which leaves no Maxwellian to draw.
"""
function plot_collision_distribution end

@doc raw"""
    plot_drift(t, (mass, momentum, energy); labels = ("mass", "momentum", "energy")) -> Figure

The relative drift of each conserved quantity against `t`, one line per quantity, on a
logarithmic `y` axis.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_drift(t, (mass, momentum, energy))
```

Each element of the second argument is the series of one quantity ``Q`` at the times `t`, and
its line is the drift relative to the initial value, ``|Q(t) - Q(t_1)| / |Q(t_1)|``.
`labels` names the lines in a legend, one label per quantity, in the same order.

The drift at `t[1]` is zero, which a logarithmic axis cannot show, so a zero drift is clipped to
the smallest positive value of the element type, `floatmin`, and the line starts at the axis'
floor. A pair whose time or drift is `NaN` or `Inf` is dropped, and the rest of the line is drawn.
A series of one value has the single drift zero, and is refused.

# Throws

- `ArgumentError` when the second argument holds no quantity, when `t` or a quantity holds no
  value, when the initial value of a quantity is zero or not finite, so that its relative drift is
  not defined, when no pair of time and drift is finite, and when a drift is the single value
  zero.
- `DimensionMismatch` when a quantity is not as long as `t`, and when `labels` does not hold one
  label per quantity; the message names the shorter one.
"""
function plot_drift end

@doc raw"""
    plot_entropy(t, S) -> Figure

The entropy ``S_h`` against `t`, on a linear axis.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_entropy(t, S)
```

A pair whose time or entropy is `NaN` or `Inf` is dropped, and the rest of the line is drawn. A
series of one value is drawn as a marker.

# Throws

- `ArgumentError` when `t` or `S` holds no value, and when no pair of `t` and `S` is finite.
- `DimensionMismatch` when `t` and `S` differ in length; the message names the shorter one.
"""
function plot_entropy end

@doc raw"""
    plot_error(t, err, closed) -> Figure

The error `err` of a run against its closed form, and the closed form `closed` itself, against
`t`: one line each, labelled `error` and `closed form`, on a logarithmic `y` axis.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_error(t, abs.(T_run .- T_closed), T_closed)
```

A pair whose time or value is `NaN` or `Inf` is dropped, and the rest of the line is drawn. A zero
in a series of more than one value is clipped to the smallest positive value of the element type,
`floatmin`, the axis' floor. A series of one value is drawn as a marker.

# Throws

- `ArgumentError` when `t`, `err` or `closed` holds no value, when no pair of `t` and `err`, or
  of `t` and `closed`, is finite, when `err` or `closed` holds a negative value, and when either
  is the single value zero: `log10` has no value there.
- `DimensionMismatch` when the three differ in length; the message names the shorter one.
"""
function plot_error end

@doc raw"""
    plot_convergence(h, err) -> Figure

The error `err` against the cell width `h`, one marker per pair, on logarithmic axes, with the
fitted line ``C h^p`` and the exponent `p` in the `y` axis' label.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_convergence([0.4, 0.2, 0.1, 0.05], errors)
```

`p` and ``\ln C`` are the slope and the intercept of the least-squares line through the pairs
``(\ln h_i, \ln e_i)`` whose `h` and `err` are positive and finite, and the label shows `p` to
three significant digits.

A pair whose `h` or `err` is `NaN` or `Inf` is dropped. A zero is drawn at the smallest positive
value of the element type, `floatmin`, the axis' floor, and takes no part in the fit, because it
has no logarithm.

# Throws

- `ArgumentError` when `h` or `err` holds no value, when no pair of them is finite, when fewer
  than two pairs are positive and finite, when every `h` of those pairs is the same, so that no
  slope can be fitted, and when either holds a negative value.
- `DimensionMismatch` when `h` and `err` differ in length; the message names the shorter one.
"""
function plot_convergence end

@doc raw"""
    plot_cost(cost, err) -> Figure

The error `err` of each run against its wall-clock cost `cost`, one marker per run, on
logarithmic axes.

Defined in a package extension: load `CairoMakie` to make it available.

```julia
using VlasovMethods, CairoMakie
fig = plot_cost(seconds, errors)
```

A run whose cost or error is `NaN` or `Inf` is dropped. A zero in a series of more than one value
is clipped to the smallest positive value of the element type, `floatmin`, the axis' floor.

# Throws

- `ArgumentError` when `cost` or `err` holds no value, when no pair of them is finite, when
  either holds a negative value, and when either is the single value zero: `log10` has no value
  there.
- `DimensionMismatch` when `cost` and `err` differ in length; the message names the shorter one.
"""
function plot_cost end
