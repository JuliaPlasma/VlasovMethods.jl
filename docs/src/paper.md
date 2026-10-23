# Reproducing the paper's examples

The reduced-order paper has two examples: the analytic quadratic model, whose solution has a closed
form, and the self-consistent bump-on-tail problem. Four scripts in `scripts/` compute both. This
page says how to run them and what each one writes. It runs no script, and it shows no figure: the
figures are files that the scripts write.

## Running a script

The scripts have their own environment, `scripts/Project.toml`, which takes VlasovMethods from the
checkout it is in. Instantiate it once, from the root of the repository:

```
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
```

Every script takes `--runs-dir` and `--results-dir`. A script writes its data, the HDF5 files, to
the first, and its figures and tables to the second. The defaults are `runs/` and `results/` at the
root of the repository, which git ignores. A flag that a script does not know stops it with an
`ArgumentError`, which names the flags it takes.

```
julia --project=scripts scripts/bump_on_tail_training.jl --runs-dir /tmp/runs --results-dir /tmp/results
```

Each figure is a PNG named `item<NN>_<what>.png`, after the item of the paper that it draws. A
figure plots the quantity of its item at the script's parameters. It does not reproduce the
paper's numbers: the paper's draws ran on Julia 1.3.1, and no value is pinned against them. Every
HDF5 file carries the root attributes `commit`, the commit of VlasovMethods that wrote it, and
`command`, the script and its arguments.

## The bump-on-tail problem

Three scripts run in this order, each with the same options, because the options set the name of
the files that the next script reads. Each takes the options of the paper's training run: `--dt`
(`0.1`), `--final-time` (`25.0`), `--particles` (`10000`), `--cells` of the spline potential
(`16`), `--degree` of the spline (`3`) and `--seed` (`1234`). The particles are drawn once, at the
reference wave number `κ₀ = 0.3`, and a sample at the wave number `κ` marches them on the same
domain with the scale `χ = κ/κ₀`. Every default is the paper's value, so each script runs the
paper's parameters without options.

### `bump_on_tail_training.jl`

The training sweep: the full-order march of the ten wave numbers `κ` evenly spaced in `[0.1, 0.5]`.

- **Paper command:** `julia --project=scripts scripts/bump_on_tail_training.jl`
- **Measured run time:** 15 s.
- **Writes:** `<run id>.h5` in the runs directory, the snapshots of the sweep with the
  accelerations of every saved state, and the integrator parameters in its group `integrator`,
  with the root attribute `seed` beside `commit` and `command`.

### `bump_on_tail_projections.jl`

The reduced bases of the training sweep: the cotangent-lift basis of the particles, and the POD
basis of the field. Each keeps the leading vectors that carry all but its tolerance of the energy
of its spectrum: `--particle-tolerance` (`1e-9`) and `--field-tolerance` (`1e-5`).

- **Paper command:**
  `julia --project=scripts scripts/bump_on_tail_projections.jl --particle-tolerance 1e-9 --field-tolerance 1e-5`
- **Measured run time:** 34 s, with the bases of `k = 205` particle vectors and `kₑ = 570` field
  vectors.
- **Writes:** in the runs directory, `<run id>_projections.h5`, the particle basis with the field
  basis in its group `field`, and `<run id>_vectorfield.h5`, the reduced positions `X` and
  accelerations `A` of every saved state, each of size `(k, nₛ, nparam)`. Both carry the root
  attributes `commit`, `command` and `seed`. In the results directory, the two spectra:
  `item11_particle_spectrum.png` and `item11_field_spectrum.png`.

### `bump_on_tail_testing.jl`

The test sweep: ten wave numbers drawn uniform in `[0.1, 0.5]` after the seed, each marched by the
full-order model (HFM) and by the reduced model of the two bases with the DEIM field (ROM).
`--samples` (`10`) sets their number.

- **Paper command:** `julia --project=scripts scripts/bump_on_tail_testing.jl`
- **Measured run time:** 37 s.
- **Writes:** in the results directory, `growth_rates.csv`, the growth rate of the field energy of
  each test sample, HFM and ROM, under the header `kappa,gamma_hfm,gamma_rom`, and the figures:

  | item | files | what they plot |
  |:--|:--|:--|
  | 24 | `item24_<error>_time.png`, `item24_<error>_kappa.png` | the four trajectory errors `x2`, `xmax`, `v2` and `vmax` against time for every sample, and against `κ` at `T`, `0.75T` and `0.5T` |
  | 25 | `item25_potential_error.png` | the absolute error of the potential against time |
  | 26 | `item26_field_energy.png`, `item26_field_energy_best.png`, `item26_field_energy_worst.png` | the field energy, HFM against ROM, of every sample, and of the samples whose growth rates agree best and worst |
  | 27 | `item27_energy.png`, `item27_energy_drift.png`, `item27_energy_drift_kappa.png` | the energy and its relative drift, HFM against ROM, and the size of the drift against `κ` at `T`, `0.75T` and `0.5T` |
  | 28 | `item28_growth_rate.png` | the growth rate against `κ`, HFM and ROM |
  | 29 | `item29_particles_t0.png`, `item29_particles_T_hfm.png`, `item29_particles_T_rom.png` | the particles of the sample nearest `κ₀`, at `t = 0` and, HFM and ROM, at `T` |

  The growth rate is half the slope of the least-squares line of `log W` through the first four
  maxima of the field energy `W`, or through all of them when there are fewer, with a window of 5
  saved states on each side of a maximum. A sample with fewer than two maxima has no rate, `NaN`,
  and item 26 ranks only the samples that have one.
  The errors and the field energy are drawn on log axes, as the reference draws them; the errors
  from the first save on, because at `t = 0` they are those of the projection alone.

## The analytic quadratic model

### `analytic_quadratic.jl`

The model `H(q, p) = ½|p|² + ½β²|q|²` of 1000 particles, trained on the closed form of six `β` in
`[5.95, 6.05]` and tested at `β = 6`. The full model is the closed form on the step `0.001`. The two
reduced forms march the projected initial condition on the step `0.0001`: the ODE form in the POD
basis of `k = 5` with RK4, and the PODE form in the cotangent-lift basis of `k = 20` with Lobatto
IIIA–IIIB(2). `--final-time` and `--saves` set the final time and the number of saved states after
`t = 0`, and `--particles` (`1000`) the number of particles.

- **Paper command:**
  `julia --project=scripts scripts/analytic_quadratic.jl --final-time 1000 --saves 10000`.
  The default is `--final-time 100 --saves 1000`, a tenth of the paper's steps.
- **Measured run time:** 135 s, with the defaults.
- **Writes:** in the runs directory, `analytic_full.h5`, `analytic_ode_k5.h5` and
  `analytic_pode_k20.h5`, each with the saved times `t`, the positions `q`, the velocities `p` and
  the Hamiltonian `H` of every saved state, and the root attributes `dt` and `nsave`; the two
  reduced files also hold the reduced states `z`. In the results directory, the two spectra of the
  training set, `item11_analytic_pod_spectrum.png` and `item11_analytic_cotangent_spectrum.png`;
  the errors of both reduced forms against the closed form, `item24_analytic_x.png` and
  `item24_analytic_v.png`; and their energy drift, `item27_analytic_energy_drift.png`. A reduced
  state that is not finite stops the run with an error.

## How the run times were measured

Each time is one run of the script with its defaults, in a new Julia process, after
`Pkg.instantiate()` and `Pkg.precompile()` of the scripts environment: Julia 1.13.2 on an Apple
silicon Mac, with other jobs running beside it. The time includes the load of the packages and the
compilation of the script's own calls.
