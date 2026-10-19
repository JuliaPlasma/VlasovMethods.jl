```@meta
CurrentModule = VlasovMethods
```

# The paper's two examples

The two example sets of the reduced-basis paper are scripts in `scripts/`, on the functions of this
package: the **self-consistent bump-on-tail** set, whose trajectories come from the full-order and
the reduced Vlasov–Poisson march, and the **analytic quadratic** set, whose trajectories are the
closed form of ``\ddot q = -\beta^2 q``. Each set writes the quantities of the paper's figures as
PNGs, and the quantities of the next section as HDF5 files.

The scripts run against the environment of `scripts/`, not against the package's test environment:

```bash
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts scripts/bump_on_tail_training.jl
```

`scripts/Project.toml` points at this package by its path, and takes `ReducedBasisMethods` from
the registry — a version that carries `spectrumplot`, which the spectra need. The first call of the
first run precompiles the environment; the runs themselves take minutes, and the measured times
are below.

Every script takes `--runs-dir` and `--results-dir`, which default to `runs/` and `results/` at the
repository root; both are ignored by git. A script writes nothing outside them. A script's header
gives its flags, and every default is the paper's own parameter, except where a note says so.

## The bump-on-tail set

The three scripts run in order. `bump_on_tail_training.jl` writes `<runid>.h5` into the runs
directory and the projections script reads it from there; the testing script reads the two basis
files. One `--runs-dir` serves all three.

### `bump_on_tail_training.jl`

The ten training χ of the paper's sweep, each a full-order march of one seeded initial condition,
saved into one `Snapshots`.

```bash
julia --project=scripts scripts/bump_on_tail_training.jl
```

Measured cold with the environment precompiled: **26.6 s**. Writes

- `<runid>.h5` — the `Snapshots` of the sweep, the `IntegratorParameters`, and the attributes
  `commit`, `command` and `seed`;
- `bump_on_tail_training_energy.png` — the full-order field energy `W(t)` of every sample, the
  training side of item 26.

### `bump_on_tail_projections.jl`

The reduced bases of the training sweep — a cotangent-lift EVD of the stacked positions and
velocities at `tol = 1e-4`, a POD of the accelerations at `tol = 1e-2` — their spectra, and the
input files of the learned reduced field.

```bash
julia --project=scripts scripts/bump_on_tail_projections.jl
```

Measured cold with the environment precompiled: **28.0 s**. Writes

- `<runid>_projections.h5` — the particle basis and, under `field`, the field basis;
- `<runid>_vectorfield.h5` — the reduced positions `X` and accelerations `A`, each
  `(kₚ, nₛ, nparam)`, with the attributes `commit`, `command` and `seed`;
- `bump_on_tail_spectrum_particles.png` and `bump_on_tail_spectrum_field.png` — the spectra of the
  two bases, item 11.

### `bump_on_tail_testing.jl`

The ten test κ of the paper's seeded draw, each marched in full and in reduced order — the reduced
field hyper-reduced by DEIM — and the paper's diagnostics of the two.

```bash
julia --project=scripts scripts/bump_on_tail_testing.jl
```

Measured cold with the environment precompiled: **29.9 s**. Writes

- `bump_on_tail_errors_time.png` — the RMS and `L∞` errors of positions and velocities per saved
  time, item 24;
- `bump_on_tail_errors_parameter.png` — the RMS position error against χ at `T`, `0.75T` and
  `0.5T`, item 24;
- `bump_on_tail_potential_error.png` — `‖Φ_ROM − Φ_HFM‖`, item 25;
- `bump_on_tail_energy.png` — the field energy of the full and the reduced march at the best and
  the worst sample of the growth rate, item 26;
- `bump_on_tail_energy_total.png` and `bump_on_tail_energy_drift.png` — the total energy `H = K + W`
  and its relative drift, item 27;
- `bump_on_tail_growth_rate.png` — the fitted `γ` of the full and the reduced sweep against χ,
  item 28;
- `bump_on_tail_phase_t0.png` and `bump_on_tail_phase_T.png` — the phase-space frame of the best
  sample at the first and the last saved state, item 29.

## The analytic quadratic set

### `analytic_quadratic.jl`

The closed-form trajectories of the six training β, the cotangent-lift and POD bases of them, and
the reduced ODE form at `k = 5` and PODE form at `k = 20` at the off-sample `β = 6.00`. The closed
form is evaluated at the saved times only, so no trajectory of `2n × 10⁶` states is stored.

The default stops at `T = 100` of the paper's `T = 1000`, at the same `dt = 1e-4`, so that a run is
minutes rather than hours. The paper's run is

```bash
julia --project=scripts scripts/analytic_quadratic.jl --final-time 1000.0
```

Measured cold at the default, with the environment precompiled: **121.6 s**. Writes

- `analytic_quadratic_spectrum.png` — the spectrum of the training basis, item 11;
- `analytic_quadratic_errors.png` — the RMS and `L∞` errors of the lifted PODE trajectory against
  the closed form, item 24;
- `analytic_quadratic_energy_drift.png` — the relative drift of the Hamiltonian of the ODE and
  PODE forms, the frozen-`ϕ` side of item 27.

## What this page does not run

The docs build runs no script and embeds no figure: a default run of the bump-on-tail set writes
hundreds of megabytes of HDF5 and takes minutes, and the paper's analytic run takes the paper's
`10⁷` steps. The figures here are therefore named, not drawn, and a number that a script prints is
a run its own header describes.
