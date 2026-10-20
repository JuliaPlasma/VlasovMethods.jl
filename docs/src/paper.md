```@meta
CurrentModule = VlasovMethods
```

# The reproduction scripts

The two example sets of the paper have a script each, and every script writes its data to the
runs directory and its figures to the results directory. The four scripts are:

| script | example set | what it does |
|:--|:--|:--|
| `scripts/bump_on_tail_training.jl` | bump-on-tail | the full-order march of the training sweep, and the snapshots it saves |
| `scripts/bump_on_tail_projections.jl` | bump-on-tail | the particle basis and the field basis of the sweep, and the reduced data of the HNN |
| `scripts/bump_on_tail_testing.jl` | bump-on-tail | the full-order and the reduced march of the test samples, and the diagnostics of the two |
| `scripts/analytic_quadratic.jl` | analytic quadratic | the closed form of the model, the reduced basis of its training trajectories, and the reduced march at the test parameter |

The three bump-on-tail scripts share `scripts/bump_on_tail_common.jl`, which holds the parameters
of the model, the wave numbers of a sweep, the field of a sample, the particle cloud, the name of
a run and the stride of the saved times; the analytic script includes it for the commit of the
checkout and for the stride of its saved times. The projections script reads the snapshots of the
training run, and the testing script reads the two files of the projections run, so the three run
in that order.

!!! note "The pages here run no script"
    The runs take too long for a documentation build, so this page names the commands and the
    files and embeds no figure. The figures are written as PNG files into the results directory.

## The commands

Every script takes `--runs-dir` and `--results-dir`, which redirect the two directories, and the
flags of its own. The defaults of the three bump-on-tail scripts are the paper's parameters; the
defaults of the analytic script are smaller, so that a run stays inside a few minutes, and the
third column names the flags of the paper's parameters.

| script | the default run | the paper's parameters | run time |
|:--|:--|:--|:--|
| `bump_on_tail_training.jl` | `julia --project=scripts scripts/bump_on_tail_training.jl` | the defaults: `--npart 10000 --final-time 25.0` | 29 s |
| `bump_on_tail_projections.jl` | `julia --project=scripts scripts/bump_on_tail_projections.jl` | the defaults: `--runid BoT_Np10000_k_010_050_ns10_T25` | 32 s |
| `bump_on_tail_testing.jl` | `julia --project=scripts scripts/bump_on_tail_testing.jl` | the defaults: `--runid BoT_Np10000_k_010_050_ns10_T25` | 32 s |
| `analytic_quadratic.jl` | `julia --project=scripts scripts/analytic_quadratic.jl` | `--npart 1000 --final-time 1000 --dt 1e-3 --rdt 1e-4` | 31 s |

The times are those of a cold process on the machine of the author, with Julia 1.13, and they
count the loading of the package. The run of the analytic script at the paper's parameters is far
longer and needs much more memory: the training matrix holds `2000 × 600006` doubles, about
9.6 GB, and the reduced march stores the `10⁷` steps of the 40-dimensional state twice over,
about 12.8 GB. A run of `10⁶` reduced steps at `1000` particles takes 112 s.

The scripts of one run share its name, which holds the parameters it names: the particle count,
the wave-number range of the sweep, the number of samples and the final time. A run of other
parameters has another name, and the projections and the testing script take it as `--runid`.

## The files

The training script writes the snapshots of the sweep, the positions, the velocities, the
accelerations, the potentials and the three diagnostics of every saved state, and the parameters
of the integrator:

```
runs/BoT_Np10000_k_010_050_ns10_T25.h5
```

The projections script writes the two files of the reduced model: the bases with the integrator
parameters of the training run, and the reduced data of the HNN, the positions and the
accelerations of every saved state projected onto the particle basis. Both carry the commit of
the checkout, the command and the seed of the run:

```
runs/BoT_Np10000_k_010_050_ns10_T25_projections.h5
runs/BoT_Np10000_k_010_050_ns10_T25_vectorfield.h5
```

The testing script writes the snapshots of the two marches, the full-order one and the reduced
one, and the analytic script writes the closed form and the reduced trajectory of the test
parameter:

```
runs/BoT_Np10000_k_010_050_ns10_T25_test.h5
runs/BoT_Np10000_k_010_050_ns10_T25_reduced.h5
runs/analytic_beta600_n200_T10p0.h5
```

## The figures

Every figure is drawn with [`plot_series`](@ref) or [`plot_particles`](@ref), and every script
writes the figures of the quantities of its example set. The names carry the run and the script,
so that the figures of the training run and of the test run do not overwrite each other.

| file | quantity | item |
|:--|:--|:--|
| `BoT_..._training_field_energy.png` | the field energy of every sample of the training sweep | 26 |
| `BoT_..._training_energy_drift.png` | the relative drift of the total energy of the training sweep | 27 |
| `BoT_..._training_particles_t0.png`, `BoT_..._training_particles_tT.png` | the phase space of the training sweep at the first and the last saved time | 29 |
| `BoT_..._particle_spectrum.png`, `BoT_..._field_spectrum.png` | the spectra of the particle basis and of the field basis | 11 |
| `BoT_..._testing_errors.png` | the errors of the reduced trajectory against the full-order one, against the time | 24 |
| `BoT_..._testing_error_x2.png`, `BoT_..._testing_error_xmax.png`, `BoT_..._testing_error_v2.png`, `BoT_..._testing_error_vmax.png` | the errors of the reduced trajectory against the full-order one, against the wave number at `T`, `0.75T` and `0.5T` | 24 |
| `BoT_..._testing_potential_error.png` | the error of the reduced potential against the full-order one | 25 |
| `BoT_..._testing_field_energy.png` | the field energy of the two marches | 26 |
| `BoT_..._testing_total_energy.png` | the total energy of the two marches | 27 |
| `BoT_..._testing_energy_drift.png` | the relative drift of the total energy of the two marches | 27 |
| `BoT_..._testing_growth_rate.png` | the growth rate of the two marches against the wave number | 28 |
| `BoT_..._testing_field_energy_best.png`, `BoT_..._testing_field_energy_worst.png` | the field energy of the best and of the worst test sample | 28 |
| `BoT_..._testing_particles_t0.png`, `BoT_..._testing_particles_tT.png` | the phase space of the reduced march at the first and the last saved time | 29 |
| `analytic_..._errors.png` | the errors of the reduced trajectory against the closed form | 24 |
| `analytic_..._energy_drift.png` | the relative drift of the total energy of the closed form and of the reduced march | 27 |

A figure reproduces the quantity of the paper's figure at the paper's parameters. It does not
reproduce its numbers: the paper's values are not pinned by a test, and the random draws of the
scripts differ from those of the notebooks the paper was written from.
