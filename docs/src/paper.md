# The paper's examples

The reduced models of this package are tested on two sets of examples, the analytic quadratic
model and the self-consistent bump-on-tail instability. Four scripts in `scripts/` run them and
draw the paper's quantities, each under the number that the tables below and the comments of the
scripts use: the errors of the positions and the velocities (24), the error of the potential
(25), the field energy (26), the energy and its drift (27), the growth rates (28), phase-space
frames (29) and the spectra of the bases (11). The match with the paper's figures is qualitative:
each figure plots the same quantity at the paper's parameters, and no value is pinned against the
paper.

This page builds no figure and runs no script, because the runs take minutes.

## Running the scripts

The scripts have their own environment, which develops this package from the checkout:

```
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
```

Every script takes `--runs-dir` and `--results-dir`. The numbers go to the first, the figures to
the second, and they default to `runs/` and `results/` at the repository root, which git ignores.
An unknown flag is an error.

The run times below are those of each script with its defaults, each in a new process, on an
Apple silicon Mac with Julia 1.13.2.

## The bump-on-tail set

The three scripts run in this order and share the run name `--runid`, `bump_on_tail` by default.
Their defaults are the paper's parameters: `10⁴` particles, a cubic periodic spline of `16`
cells, `dt = 0.1` up to `T = 25`, a training sweep of ten wave numbers evenly spaced over
`κ ∈ [0.1, 0.5]`, a test sweep of ten wave numbers drawn uniform in the same range, the seed
`1234`, and the tolerances `1e-4` of the particle basis and `1e-2` of the field basis. Each
sample is the problem at `κ = 0.3` with the scale `χ = κ / 0.3`, and it marches its own
`ScaledField` at that `χ`.

| script | what it does | run time | files |
|:--|:--|:--|:--|
| `bump_on_tail_training.jl` | the full-order march of the training sweep, and the accelerations the field applied | 25 s | `<runid>_training.h5`; `<runid>_training_field_energy.png` (26), `<runid>_training_growth_rate.png` (28) |
| `bump_on_tail_projections.jl` | the cotangent-lift basis of the particles, the POD basis of the accelerations and its DEIM indices, and the reduced data `X` and `A` of the training sweep | 29 s | `<runid>_projections.h5`, `<runid>_vectorfield.h5`; `<runid>_spectrum_particle.png`, `<runid>_spectrum_field.png` (11) |
| `bump_on_tail_testing.jl` | the full-order and the reduced march, with the DEIM field, of the test sweep, and the diagnostics of the one against the other | 29 s | `<runid>_testing.h5`; `<runid>_errors_x_time.png`, `<runid>_errors_v_time.png`, `<runid>_errors_x_parameter.png`, `<runid>_errors_v_parameter.png` (24), `<runid>_potential_error.png` (25), `<runid>_field_energy.png` (26), `<runid>_energy.png`, `<runid>_energy_drift.png` (27), `<runid>_growth_rate.png` (28), `<runid>_particles_initial.png`, `<runid>_particles_final_hfm.png`, `<runid>_particles_final_rom.png` (29) |

The commands with the paper's parameters are the defaults:

```
julia --project=scripts scripts/bump_on_tail_training.jl
julia --project=scripts scripts/bump_on_tail_projections.jl
julia --project=scripts scripts/bump_on_tail_testing.jl
```

Each basis keeps the leading vectors that carry all but the tolerance of the energy of its
spectrum. The default pair, which ReducedBasisMethods' projection script runs, gives
`kₚ = 10` and `kₑ = 55` on the default training run. That script names a second pair, `1e-8` and
`1e-4`, which gives `kₚ = 120` and `kₑ = 345`:

```
julia --project=scripts scripts/bump_on_tail_projections.jl --particle-tol 1e-8 --field-tol 1e-4
```

`bump_on_tail_testing.jl` fits each growth rate to the first three local maxima of the field
energy, because the reduced run of the smallest test `κ` has no more at the default pair.

`<runid>_projections.h5` holds the particle basis, the field basis and the DEIM indices, and
`<runid>_vectorfield.h5` holds the reduced coordinates `X` and accelerations `A`, each of the size
`(k, nₛ, nparam)`. Both carry the HDF5 attributes `commit`, `command` and `seed`, which name the
run that wrote them. They are the input files of the HNN of the reduced field.

## The analytic quadratic set

`analytic_quadratic.jl` runs the `n = 1000` particles of the paper in the field `E(x) = β²x`,
whose closed form is the full model. It trains a POD basis of `k = 5` and a cotangent-lift basis
of `k = 20` on six `β ∈ [5.95, 6.05]`, and tests them off-sample at `β = 6.0` with the reduced ODE
form and RK4 and the reduced PODE form and Lobatto IIIA–IIIB(2), at the reduced step `10⁻⁴`.

| script | run time | files |
|:--|:--|:--|
| `analytic_quadratic.jl` | 121 s | `analytic_quadratic.h5`; `analytic_quadratic_errors_x.png`, `analytic_quadratic_errors_v.png` (24), `analytic_quadratic_energy_drift.png` (27), `analytic_quadratic_spectrum_pod.png`, `analytic_quadratic_spectrum_cotangent_lift.png` (11) |

Its default trains on every 1000th step of `dt = 0.001` and runs to `T = 100`. The paper trains
on every 10th step up to `T = 1000`, which is a snapshot of 9.6 GB, and its PODE run alone is
`10⁷` reduced steps. The command with the paper's parameters is

```
julia --project=scripts scripts/analytic_quadratic.jl --train-step 10 --final-time 1000.0
```
