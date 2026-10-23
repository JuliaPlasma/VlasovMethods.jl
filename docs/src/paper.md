# The paper's examples

The two example sets of the reduced-order paper are four scripts in `scripts/`, run from the
repository root with the scripts environment:

```bash
julia --project=scripts scripts/bump_on_tail_training.jl --runs-dir runs --results-dir results
```

Each script takes `--runs-dir` and `--results-dir`, which default to `runs/` and `results/` at the
repository root and are ignored by git. Every file it writes goes through one of the two, so a run
leaves the checkout itself unchanged: `git status --porcelain --ignored` is the same before and
after it. Each script prints what it wrote, and an unknown flag is refused with the list of the
flags it takes.

The tables below give, for each script, the command at the paper's parameters, the files it writes
and its run time on its defaults. The times are wall-clock seconds of one cold process, measured on
an Apple M-series machine with Julia 1.13.2, after `Pkg.instantiate()` of `scripts/` has
precompiled the environment; a first run that precompiles is longer.

## The bump-on-tail example

The model is the bump-on-tail of `κ = 0.3` on the periodic domain `2π/κ`, and the parameter of its
sweep is the time scale `χ = κₛ/κ` of ten wave numbers `κₛ ∈ [0.1, 0.5]`: one physical run, the
same particles and the same potential, marched at ten time scales. The three scripts share
`scripts/bump_on_tail_common.jl`.

| script | writes | time |
|:--|:--|--:|
| `scripts/bump_on_tail_training.jl` | `<runid>.h5`: the ten training samples | 13 s |
| `scripts/bump_on_tail_projections.jl` | `<runid>_projections.h5`, `<runid>_vectorfield.h5`, `item11_spectra.png` | 33 s |
| `scripts/bump_on_tail_testing.jl` | `growth_rates.csv`, `item24_*`–`item29_*` figures | 35 s |

`<runid>` is a name the parameters produce — the particle and cell counts, the timestep, the final
time and the wave-number range — such as `BoT_np10000_nh16_dt0p1_T25p0_k0p1-0p5`.

### `scripts/bump_on_tail_training.jl`

Marches the ten training samples full order and writes their `Snapshots` — positions, velocities,
accelerations, potential coefficients and the three diagnostics of every saved time — to the runs
directory as `<runid>.h5`. The paper's command:

```bash
julia --project=scripts scripts/bump_on_tail_training.jl --np 10000 --nh 16 --order 3 \
    --dt 0.1 --final-time 25.0
```

The default is that command, at about 13 s and a 600 MB file.

### `scripts/bump_on_tail_projections.jl`

Builds the two reduced bases of the training set: the cotangent-lift basis of the positions and
velocities at the particle tolerance, and the POD basis of the accelerations at the field
tolerance. The paper's command, whose tolerances are those of the reference's
`VP-ROM-Projections.ipynb` cells 14 and 15:

```bash
julia --project=scripts scripts/bump_on_tail_projections.jl --particle-tol 1e-9 --field-tol 1e-5
```

At the paper's parameters — `--np 10000 --nh 16 --order 3 --dt 0.1 --final-time 25.0` with those
two tolerances — it retains `k = 209` particle modes and `k = 750` field modes on the machine named
above. It writes the particle basis and the field basis into `<runid>_projections.h5`, the reduced
positions and accelerations of every saved state into `<runid>_vectorfield.h5`, each with the root
attributes `commit`, `command` and `seed`, and the particle and field spectra to
`item11_spectra.png`. The default is the command above, at about 33 s.

### `scripts/bump_on_tail_testing.jl`

Marches the ten test samples — ten wave numbers drawn uniform in `[0.1, 0.5]` after
`Random.seed!(1234)`, sorted, each read as the time scale `χ = κ/0.3` — full order and reduced, and
writes their diagnostics. The paper's command:

```bash
julia --project=scripts scripts/bump_on_tail_testing.jl --dt 0.1 --final-time 25.0
```

The default step is `0.05`: the growth-rate fit reads the first four maxima of the field energy
over a window of five saved states, and at `0.1` two of the ten drawn time scales resolve so
coarsely that their fits turn negative. The script reads the particle basis that
`bump_on_tail_projections.jl` wrote, so run that script first; it finds it as the runs directory's
single `*_projections.h5`. It writes `growth_rates.csv`, whose header is
`kappa,gamma_hfm,gamma_rom`, one row per test sample, and the figures `item24_*`,
`item25_*`, `item26_*`, `item27_*`, `item28_growth_rate.png` and
`item29_particles_*.png`. The defaults take about 35 s.

## The analytic quadratic model

The model is the closed-form motion of `n` particles in the harmonic potential `β²x²/2`, whose
reduced forms the paper integrates at `k = 5` in a POD basis with RK4 and at `k = 20` in a
cotangent-lift basis with Lobatto IIIA–IIIB(2).

### `scripts/analytic_quadratic.jl`

The paper's command, from `reference/Vlasov_Analytic.jl`, `Reduced_Vlasov_ODE.jl` and
`Reduced_Vlasov_PODE_CotangentLift.jl`:

```bash
julia --project=scripts scripts/analytic_quadratic.jl --final-time 1000 --dt 0.001 \
    --dt-reduced 0.0001 --nsave 10 --nsave-reduced 1000
```

That is `dt = 0.001` for the full model and `dt = 0.0001` for the reduced one, both to `T = 1000`,
with the reduced runs saved every 1000th step — the 10 000 saved states of the paper. The default
is smaller, `--final-time 20 --dt 0.001 --dt-reduced 0.001 --nsave 10 --nsave-reduced 10`, at about
40 s.

It writes `analytic_full.h5`, the closed-form states with a `t` dataset and the root attribute
`dt`, and `analytic_ode_k5.h5` and `analytic_pode_k20.h5`, the reduced states of the two forms with
a `t` dataset each, and the spectrum of the snapshots to `item11_analytic_spectrum.png`.

## What the files hold

The figures are named `item<NN>_<what>.png` by the acceptance item of the paper they draw: item 11
the particle and field spectra, item 24 the trajectory errors against time and against `κ` at `T`,
`0.75T` and `0.5T`, item 25 the absolute potential error against time, item 26 the field energy,
item 27 the energy and its relative drift, item 28 the growth rate against `κ`, and item 29 the
particles at `t = 0` and at `T`. Every series is drawn with `plot_series` — the particle and field
spectra included — every particle frame with `plot_particles`, and the spectrum of the analytic
snapshots with `ReducedBasisMethods.spectrumplot`.

The HDF5 files of the bump-on-tail run are the input of the neural-network part of this work: the
projections file holds the two bases and the vectorfield file the reduced positions `X` and
accelerations `A` of every training sample, each of the size `(k, nₛ, nparam)`.
