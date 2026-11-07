```@meta
CurrentModule = VlasovMethods
```

# The figures of the collision manuscripts

Each figure of a collision manuscript has one driver in `scripts/` and one command that makes it.
A driver reads its flags through `scripts/options.jl`, defaults to the manuscript's parameters,
writes its figures into `--results-dir` and one HDF5 file per run into `--runs-dir`, and prints
its wall-clock time, the package loads included. Each HDF5 file carries the attribute `commit`,
the `git rev-parse HEAD` of the checkout or `"unknown"` outside one, the attribute `command`, the
driver's arguments joined by spaces, and one attribute per option. A directory that does not exist
is created.

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
domain, the run stops with the model's `DomainError`, the driver writes no figure for that run and
exits non-zero. A run to the time `t` with the step `h` takes `t / h` steps, rounded up, so the
quick run of figures 3 to 6 takes 63 steps of `8e-4` and ends at `t = 0.0504`.

### Figures 1 and 2

```
julia --project=scripts scripts/pd_lb_fig1_2.jl
```

Figure 1 is the norm `‖v̇‖₂` of the right-hand side at `10⁵` normal velocities, which are never
projected, with the Maxwellian projected onto `M = 24, …, 78` basis functions, against the cell
width. Figure 2 is the sample variance `Σ v̇² / N` of `N` projected normal velocities on 41 basis
functions, against `1 / N`, for `N` from `10²` to `10⁵` in steps of a quarter decade. Files:
`pd_lb_fig1.pdf`, `pd_lb_fig2.pdf`, `pd_lb_fig1.hdf5`, `pd_lb_fig2.hdf5`.

Measured: 22.1 s, and the driver exits non-zero. Both sweeps stop on the `DomainError` of a
projection that is not positive at a particle: figure 1 at `M = 25`, where the projected
Maxwellian is `−2.4e-5` at two of the `10⁵` velocities, and figure 2 at `N = 5623`. Neither
figure is written, and each HDF5 file holds the values before the stop and the `M` or the `N`
that stopped as the attribute `stopped_at`. `--nbasis` and `--sweep-npart` take other lists, such
as `--nbasis 24,26,28` or `--sweep-npart 100,178,316,562,1000,1778,3162`.

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

Not measured yet: the run takes 6250 steps, and 100 steps of the same configuration take 20 s.

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
