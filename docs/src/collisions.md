```@meta
CurrentModule = VlasovMethods
```

# The figures of the collision manuscripts

Each figure on this page has one driver in `scripts/` and one command that makes it.
A driver reads its flags through `scripts/options.jl`, defaults to the manuscript's parameters,
writes its figures into `--results-dir` and one HDF5 file per run into `--runs-dir`, and prints
its wall-clock time, the package loads included. Each HDF5 file carries the attribute `commit`,
the `git rev-parse HEAD` of the checkout or `"unknown"` outside one, the attribute `command`, the
driver's arguments joined by spaces, the attribute `script`, the driver's file name, and one
attribute per option. A directory that does not exist is created.

A driver loads `CairoMakie` and draws each figure with one of the plotting functions below; it
builds no figure itself.

The flag `--quick` replaces some of the manuscript's parameters with smaller ones.
`scripts/gate.jl` runs every driver with it, in a process of its own and into a temporary
directory, and checks that each exits 0 within 240 s, writes its figures and its HDF5 files, and
records `commit`, `command` and the parameters of `--quick`:

```
julia --startup-file=no --project=scripts scripts/gate.jl
```

The measured times below are of one run each on an Apple M4 Max with Julia 1.13.3, one driver at
a time.

## The plotting functions

```@docs
plot_collision_distribution
plot_drift
plot_entropy
plot_error
plot_convergence
plot_cost
```

## Particle discretisation of the Lenard–Bernstein operator

The drivers `scripts/pd_lb_*.jl` make the eleven figures of the Lenard–Bernstein manuscript
(PD-LB). They share `scripts/pd_lb_common.jl`. Every spline is cubic, on the velocity domain
`[−10, 10]`, and `M`, the flag `--nbasis`, is the number of its basis functions: a cubic basis on
`ncell` cells has `ncell + 3` of them, so the manuscript's "41 elements" are 41 basis functions on
38 cells. Where a projection is not positive at a particle, or a particle leaves the velocity
domain, a run stops with the model's `DomainError`: the driver writes the HDF5 file of that run
with the attribute `stopped_at`, the time of the step that threw, and the times `t` and the
velocities `z` of the steps reached, writes no figure for that run, and exits non-zero. The
sweeps of figures 1 and 2 are the exception: a point of a sweep where the operator stops is left
out of its figure, and the driver exits 0. A run to the time `t` with the step `h` takes `t / h`
steps, rounded up, so the quick run of figures 3 to 6 takes 63 steps of `8e-4` and ends at
`t = 0.0504`. K94 of `KNOWN_ISSUES.md` records each stop at the manuscript's parameters.

### Figures 1 and 2

```
julia --project=scripts scripts/pd_lb_fig1_2.jl
```

Figure 1 is the norm `‖v̇‖₂` of the right-hand side at `10⁵` normal velocities, which are never
projected, with the Maxwellian projected onto `M = 24, …, 78` basis functions, against the cell
width. Figure 2 is the sample variance `Σ v̇² / N` of `N` projected normal velocities on 41 basis
functions, against `1 / N`, for `N` from `10²` to `10⁵` in steps of a quarter decade. Files:
`pd_lb_fig1.pdf`, `pd_lb_fig2.pdf`, `pd_lb_fig1.hdf5`, `pd_lb_fig2.hdf5`.

Measured: 24.0 s, and the driver exits 0. Where the projection is not positive at a particle, the
`DomainError` of `compute_coefficients`, the operator stops at that point of the sweep: the
driver prints the stop, leaves the point out of the figure and goes on, and the HDF5 file of the
figure records the first such `M` or `N` as the attribute `stopped_at` beside the points it
computed. With the defaults, figure 1 stops at `M = 25`, where the projected Maxwellian is
`−2.4e-5` at two of the `10⁵` velocities, and draws the other 54 values of `M`; figure 2 stops at
`N = 5623` and `N = 31623`, and draws the other 11 values of `N`. A sweep that computes fewer
than two points leaves `plot_convergence` nothing to fit, and its `ArgumentError` stops the
driver.

`--quick` takes 2000 velocities, `M ∈ {24, 41, 78}` and `N` from 100 to 1778.

### Figures 3 to 6

```
julia --project=scripts scripts/pd_lb_fig3_6.jl
```

The relaxation of the normal distribution of mean 2 under `ConservativeLenardBernstein` and its
implicit midpoint rule: 1000 particles, 41 basis functions, `h = 8e-4`, `ν = 1`, to `t = 1`.
Figure 3 is the initial and the final distribution, figure 4 the relative drift of the mass, the
momentum and the energy, figure 5 the entropy `S_h`, and figure 6 the initial and the final
distribution of the same run with 200 particles. Files: `pd_lb_fig3.pdf` to `pd_lb_fig6.pdf`,
`pd_lb_fig3_5.hdf5`, `pd_lb_fig6.hdf5`.

Measured: 381.4 s, and the driver exits 0.

`--quick` takes 200 particles and `t = 0.05`.

### Figures 7 and 8

```
julia --project=scripts scripts/pd_lb_fig7_8.jl
```

The relaxation of the bi-Maxwellian of the two normal distributions of mean `±2` under
`ConservativeLenardBernstein`: 1000 particles, `h = 8e-4`, to `t = 5`, the random number generator
seeded with `--seed 1`. Figure 7 is the initial and the final distribution, and figure 8 the
entropy. Files: `pd_lb_fig7.pdf`, `pd_lb_fig8.pdf`, `pd_lb_fig7_8.hdf5`.

The run with the manuscript's parameters, 6250 steps, stops after 289 s, in the step 1768 to
`t = 1.4144`: the projection of the particles is `-1.58e-7` at `v = -6.42`, so `f_s'/f_s` of the
conservative operator is undefined and the run throws the projection's `DomainError`. As the rule
above says for a single run, the driver writes `pd_lb_fig7_8.hdf5` with `stopped_at = 1.4144`
and the 1767 steps reached, writes no figure and exits non-zero; K94 records the stop.

`--quick` takes 200 particles and `t = 0.2`.

### Figures 9 and 10

```
julia --project=scripts scripts/pd_lb_fig9_10.jl
```

The relaxation of the uniform distribution on `[−2, 2]` under `ConservativeLenardBernstein`: 200
particles, `h = 1e-4`, to `t = 1`, seeded with `--seed 1`. Figure 9 is the initial and the final
distribution, and figure 10 the entropy. Files: `pd_lb_fig9.pdf`, `pd_lb_fig10.pdf`,
`pd_lb_fig9_10.hdf5`.

Measured: 379.7 s, while other Julia processes ran on the machine, and the driver exits 0.

`--quick` takes `h = 1e-3` and `t = 0.05`.

### Figure 11

```
julia --project=scripts scripts/pd_lb_fig11.jl
```

The cumulant appendix fixes the coefficients at `A₁ = 0` and `A₂ = 1`, which is the right-hand
side of `LenardBernstein`, and runs it with the implicit midpoint scheme of the time
discretisation, which the driver uses in place of the model's own `RK438`. The bi-Maxwellian of
figure 7 with 1000 particles, 41 basis functions, `h = 5e-3`, to `t = 1`, seeded with `--seed 1`.
The closed form of the appendix is `κ₄(t) = κ₄(0) exp(−4 ν t)` at the equilibrium variance
`1 / A₂ = 1`. The driver prints the relative error of the particles' fourth cumulant against it at
the first three steps, and exits non-zero where one exceeds `--tolerance`, 5 %. Figure 11 is the
error and the closed form, both relative to `κ₄(0)`. Files: `pd_lb_fig11.pdf`, `pd_lb_fig11.hdf5`.
With `--seed 1234`, the implicit midpoint step sends a particle out of the velocity domain at
`t = 0.035`, and the run stops (K94 of `KNOWN_ISSUES.md`).

Measured: 79.0 s, and the driver exits 0. The three relative errors are `8.4e-4`, `1.5e-3` and
`1.8e-3`.

`--quick` takes 200 particles, `h = 1e-2` and `t = 0.1`.
