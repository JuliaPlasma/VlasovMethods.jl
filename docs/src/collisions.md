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

## Metriplectic particle methods for the Lenard–Bernstein operator

The driver `scripts/lenard_bernstein_metriplectic.jl` runs the three tests of the metriplectic
Lenard–Bernstein manuscript (MP-LB) under `MetriplecticLenardBernstein`, with the parameters of
the PD-LB figures of the same distribution, since MP-LB states none: 41 cubic basis functions on
`[−10, 10]` and `ν = 1` throughout, and

- `--case shifted-normal`, the default: 1000 velocities from the normal distribution of mean 2
  and variance 1, `h = 8e-4`, to `t = 1`, as figures 3 to 5;
- `--case double-maxwellian`: 1000 velocities, half from the normal distribution of mean 2 and
  half from that of mean `−2`, `h = 8e-4`, to `t = 5`, seeded with `--seed 1`, as figures 7 and 8;
- `--case uniform`: 200 velocities from the uniform distribution on `[−2, 2]`, `h = 1e-4`, to
  `t = 1`, seeded with `--seed 1`, as figures 9 and 10.

`--bracket` takes `dougherty`, the default, for [`DoughertyBracket`](@ref), the bracket of
MP-LB, or `jeyakumar` for [`JeyakumarBracket`](@ref). A case or a bracket that is not one of these
stops the driver with an `ArgumentError` that names the flag and the value, before any run. Each
step is the model's implicit midpoint step, solved by `Picard_iterate_over_particles!`.

A run writes `lenard_bernstein_metriplectic_<case>_<bracket>.hdf5` and three figures with the same
stem: `_initial_final_distribution.pdf`, `_conservation.pdf`, the relative drift of the mass, the
momentum and the energy, and `_entropy_log.pdf`, the relative change `|S_h(t) − S_h(0)| / |S_h(0)|`
of the entropy on a logarithmic axis. A step that throws, the `DomainError` of the projection or
the error of a Picard solve that does not converge, stops the run as a single run of PD-LB stops:
the HDF5 file carries `stopped_at`, the steps reached and, with the shifted normal, the times and
temperatures of the rate check, no figure is written, and the driver exits non-zero.

### The six commands

```
julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case shifted-normal --bracket dougherty
julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case shifted-normal --bracket jeyakumar
julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case double-maxwellian --bracket dougherty
julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case double-maxwellian --bracket jeyakumar
julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case uniform --bracket dougherty
julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case uniform --bracket jeyakumar
```

At these parameters every run stops where the Picard solve of a step does not converge, writes
its HDF5 file with `stopped_at`, writes no figure and exits non-zero; K95 of `KNOWN_ISSUES.md`
records the stops:

| case | bracket | stops at | steps reached | measured |
|:--|:--|--:|--:|--:|
| shifted normal | Dougherty | `t = 0.0048` | 5 of 1250 | 28 s |
| shifted normal | Jeyakumar | `t = 0.2208` | 275 of 1250 | 27 s |
| double Maxwellian | Dougherty | `t = 0.1248` | 155 of 6250 | 21 s |
| double Maxwellian | Jeyakumar | `t = 0.2968` | 370 of 6250 | 21 s |
| uniform | Dougherty | `t = 0.0399` | 398 of 10000 | 21 s |
| uniform | Jeyakumar | `t = 0.8046` | 8045 of 10000 | 22 s |

`--quick` takes 200 particles and `t = 0.05` for the shifted normal, `t = 0.05` for the double
Maxwellian, which keeps its 1000 particles because with 200 the solve of the Jeyakumar form does
not converge in the step to `t = 0.016`, and `h = 1e-3` and `t = 0.05` for the uniform case, as
figures 9 and 10 take. Each of the six quick runs exits 0.

### The rate check

With the shifted normal, the driver first runs the rate check of its bracket. The bracket's
density `n₀` and mean `u₀` are frozen at those of the initial particles and its temperature at
`T₀ = 2`, the coefficients of the Ornstein–Uhlenbeck limit. The particles are scaled about their
mean to the temperature `T₀ / 2 = 1`, from which they relax to `T₀`. In the continuum limit,
where `S_γ / w_γ` is `∂ᵥ log f`, the Dougherty form is
`v̇ = −ν T₀ ∂ᵥ log f − ν (v − u₀)`, and the Jeyakumar form is the same times `n₀ / T₀`, so the
temperature relaxes to `T₀` as `exp(−λ t)` with

```math
\lambda_\mathrm{Dougherty} = 2 \nu , \qquad
\lambda_\mathrm{Jeyakumar} = \frac{2 \nu n_0}{T_0} = \nu n_0 .
```

The Jeyakumar rate carries the temperature and the Dougherty rate does not. The weights of the
particles sum to `n₀ = 1`, so at `T₀ = 2` the two closed forms differ by two; at a temperature
near `n₀`, such as the unit variance of the sample, they would be equal up to its sampling
error. The driver fits `λ` to `ln(T₀ − T(t))` over the run, stops with an `ArgumentError` where
a temperature reaches or crosses `T₀`, and exits non-zero where the fit differs from the closed
form of its bracket by more than `--tolerance`, 5 %, or from the control, the closed form of the
other bracket, by less. Only the shifted normal takes `--tolerance`. With the parameters above:

| bracket | fitted `λ` | closed form | relative difference | control | control's relative difference |
|:--|--:|--:|--:|--:|--:|
| Dougherty | 1.989115 | `2ν` = 2.000000 | `5.4e-3` | `ν n₀` = 1.000000 | `9.9e-1` |
| Jeyakumar | 0.993993 | `ν n₀` = 1.000000 | `6.0e-3` | `2ν` = 2.000000 | `5.0e-1` |

Each step of the check solves the implicit midpoint rule by a fixed-point iteration damped by
`1/2`. Where a particle in a sparse tail sits where `f_s` is small, the step is stiff and the
iteration stalls; after 1000 iterations the step takes the last iterate, and the driver prints
how many steps did so: 1 of the 1250 for the Dougherty form and 2 for the Jeyakumar form, at a
largest residual of `3.1e-3` and `2.1e-4`. The quick run, 200 particles to `t = 0.05`, fits
`1.969792` against `2` and `0.998385` against `1`, within `1.5 %` and `0.16 %`, with controls
`97 %` and `50 %` away, and no step at the cap. With `--quick --npart 2000` the fits are
`1.999774` and `1.000624`, with controls `100 %` and `50 %` away.

The check runs on the shifted normal only, the default case, because the closed form is not met
at the stated tolerance on every case. With the same scaling, the same frozen `T₀ = 2` and the
rate check's solver, the fits differ from their closed forms by:

| case | parameters | Dougherty | Jeyakumar |
|:--|:--|--:|--:|
| double Maxwellian | 1000 particles, `--seed 1`, `h = 8e-4`, `t = 0.05` | `4.2e-3` | `3.3e-3` |
| uniform | 200 particles, `--seed 1`, `h = 1e-3`, `t = 0.05` | `8.5e-2` | `9.6e-3` |
| uniform | 200 particles, `--seed 1`, `h = 1e-3`, `t = 1` | `3.7e-2` | `6.9e-2` |

The uniform case misses the 5 % tolerance at one bracket of each set, and the double Maxwellian
meets it at both; the controls lie between `49 %` and `99 %` away in every row.
