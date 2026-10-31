# Release Notes

All notable changes to VlasovMethods.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and holds the releases from `0.3.0` on. 5 versions were
released before it, the most recent `v0.2.1` (2023-12-04), and none of them are written up here:
the record of that history is `git log` and the tags. It is named as a gap rather than
reconstructed, because a changelog assembled after the fact loses exactly the reasoning that
makes it worth keeping. The 69 commits from `v0.2.1` to this file's first commit (`f26be8b`,
"Add shared git hooks and seed CHANGELOG", 2026-08-31) ship in `0.3.0` without a write-up of
their own, except for what they removed and narrowed, which is listed in full below against
`v0.2.1`: four exports leave `names(VlasovMethods)` and eight `[compat]` entries lost versions.

## [Unreleased] — targeting 0.4.0

### New Features

- **`ScaledField`.** `ScaledField(field, χ)` wraps any `ElectricField` and divides the field of
  `efield!` and the value of `energy` by `χ²`, while `update!` and `coefficients` forward to
  `field` unchanged, because the deposit and the solve carry no amplitude. It is an
  `ElectricField`, so a reduced march wraps it as `ReducedElectricField(ScaledField(f, χ), Ψ, X₀)`
  and marches in the scaled time. `integrate_vp` and `reduced_integrate_vp` take the scale from
  `params.χ`, the one source of `χ`, and refuse a `ScaledField` built at another `χ` (next entry).
  The energy divides by `χ²` and not by `χ⁴`, because `wᵢ·E(xᵢ)/χ²` is the force of the energy
  `W/χ²`. A `χ` that is not finite and positive throws an
  `ArgumentError`, and each division converts `χ` to the element type of the value it divides, so
  a `Float32` field with a `Float64` `χ` stays `Float32`.

- **A march refuses a `χ` that it would discard.** `integrate_vp` and `reduced_integrate_vp` throw
  an `ArgumentError` that names both values when the field they are given is a `ScaledField`, or a
  `ReducedElectricField` over one, built at another `χ` than `params.χ`. A march takes its one `χ`
  from `params.χ`, for the timestep and the field together, so a stored `χ` that differs would be
  discarded in silence. A `ScaledField` built at `params.χ`, where an integer `χ` equals its
  float, and a plain field march as before. The message prints each value with `repr`, so a
  `Float32` `0.1f0` against a `Float64` `0.1` reads as two values. A `ScaledField` over a
  `ScaledField` is refused at any `χ`, the march's own included: it divides the field by both
  `χ²`, and the comparison of the outer `χ` alone would let it through.

- **`integrate_vp`.** The full-order leapfrog march, `reduced_integrate_vp` with no reduced basis
  and no cache: the same `Δt = IP.dt·χ`, the same save slots, and the same `X`, `V`, `Φ`, `W`, `K`
  and `M` written into the `Snapshots`. No function in `src/` filled a `Snapshots` from a
  full-order march before, so the full-order side of every ROM-against-HFM diagnostic had no
  producer. The march copies `P₀` and holds its own buffers, so it writes none of the caller's
  state, and on an identity basis it agrees with `reduced_integrate_vp` bit for bit.

- **`trajectory_errors`.** `trajectory_errors(Xh, Vh, Xr, Vr)` is the error of a full-order state
  against a reduced one, as the `NamedTuple` `(x2, xmax, v2, vmax)`: `x2` and `v2` the root mean
  square over the particle axis, the root sum of squares divided by `√nₚ`, so they do not change
  when the particles are duplicated, `xmax` and `vmax` the largest
  absolute difference, the `L∞` error. Each argument takes one saved state per `(ts, p)` as a
  `(nₚ, nₛ, nparam)` array or a `(nₚ, nₛ)` matrix, and each result has the trailing shape of its
  argument.

- **`potential_error`.** The error of a reduced potential against a full-order one: the root sum of
  squares of `Φh - Φr` over the coefficient axis, per saved time and parameter. It is an absolute
  error, with no divisor.

- **`field_energy_error`.** `abs.(Wh .- Wr) ./ Wh`, the relative error of a reduced field energy
  per saved time and parameter. A field energy is `½ϕᵀKϕ ≥ 0` by construction, so no guard is
  needed: a `Wh` of `0` gives `Inf`, or `NaN` where `Wr` is also `0`, and both propagate.

- **`energy_drift`.** `(H .- H[1]) ./ H[1]`, the energy at each saved time against the initial one.
  The energies of a conservative model are constant, so the drift is its discretisation error. An
  empty `H`, and an `H[1]` of `0`, throw an `ArgumentError`, because a relative drift needs a
  non-zero reference.

- **`saved_indices`.** The three saved states that a plot of one trend uses: the final one, the one
  nearest three quarters of the way and the one nearest halfway. `saved_indices(10001)` is
  `(10001, 7501, 5001)`. A `nₛ` below `2` throws an `ArgumentError`.

- **`growth_rate`.** `growth_rate(t, W; window = 2, maxima = 2:5)` is the field's growth rate `γ`,
  **half** the least-squares slope of `log W` over the selected local maxima, because
  `W ∝ |E|² ∝ e^{2γt}`. A maximum is a point no lower than the `window` points on each side whose
  first differing sample on its right is lower, so a flat maximum counts once, at its first index,
  and a rising shoulder is no maximum. The notebook this package's reduced-order work follows
  plotted the slope itself under the name `γ`, so the rates here are half the figures it drew.
  The fit runs in
  `eltype(W)` in closed form, so the package takes no `Statistics` dependency. A `t` and a `W` of
  different lengths throw a `DimensionMismatch`. A `window` below `1`, a `maxima` that repeats an
  index, selects fewer than two maxima, reaches past the maxima found or reaches an index below
  `1`, and a `W` that is not positive at a selected maximum, where `log` is not defined, each throw
  an `ArgumentError`. This closes K63 and K64 in `KNOWN_ISSUES.md`, which recorded a `BoundsError`
  at these inputs.

- **`best_and_worst`.** The indices of the smallest and of the largest of `abs.(γr .- γh)`: the
  parameters of a sweep whose growth rates agree and disagree most. A tie takes its first index.

- **`plot_series`.** `plot_series(t, series; labels = nothing, xlabel = "t", ylabel = "")` draws
  one line per element of `series` against `t`, so one function plots the error norms of
  `trajectory_errors` and `potential_error`, the field energies and the energy drift of a march
  against the saved times, and the growth rates of a sweep against the equation parameter.
  `labels` names the lines in a legend, one label per series; a count that differs throws a
  `DimensionMismatch`. The axis is linear, so set a log scale on the returned `Figure` to read a
  decaying error. It is defined in the package extension, so `CairoMakie` must be loaded.

- The `plot_particles` method on a `Snapshots` draws one saved state,
  `SS.X[1, :, ts, p]` against `SS.V[1, :, ts, p]`, so one frame of a march is one call.

- **`train_hnn`, the HNN of the reduced field.** `train_hnn(model, ps₀, data, target; η, epochs,
  batchsize, rng, batches = nothing)` trains a three-layer Hamiltonian network on the reduced
  positions and accelerations that `bump_on_tail_projections.jl` writes, with the step from
  `GeometricOptimizers` (`GradientMethod()`, `linesearch = η`), and returns the parameters and the
  loss of every epoch. The parameter gradient of the loss is written out by hand, as reverse mode
  through the network's input gradient, for `Chain`s of `Dense` layers of `tanh` or `identity`; an
  epoch allocates nothing, and the first call compiles in about a second at any width. It is a
  package extension on `GeometricOptimizers` and `AbstractNeuralNetworks`, two weak dependencies
  bounded `0.9.1` and `0.9.3`: loading this package does not load `GeometricOptimizers`, and
  `train_hnn` has no method until both are loaded. `AbstractNeuralNetworks` is loaded with
  `ReducedComplexityModeling`, which depends on it. The curve is checked against a reference
  trainer in `test/helpers/hnn_reference.jl`, which computes the same loss with nested forward-mode
  derivatives from ForwardDiff: from the same weights and batches at `m = 10`, the trainer matches
  its first epoch to `1e-12` and every one of its 1000 epochs to `1e-8`, and leaves it when the step
  size changes by a hundredth of itself. The driver `scripts/bump_on_tail_train_hnn.jl` reads the
  two files of the projections script and draws the curve at its default `m = 10`, and `hnn_data`
  in `scripts/bump_on_tail_common.jl` forms its training data.

- **The paper's two examples as scripts.** `scripts/bump_on_tail_training.jl`,
  `bump_on_tail_projections.jl` and `bump_on_tail_testing.jl` run the bump-on-tail example in three
  steps: the full-order training sweep of ten wave numbers, the cotangent-lift basis of the
  particles and the POD basis of the field, and the marches of ten test wave numbers by the
  full-order model and by two reduced models, with the DEIM field and with the full field of the
  lifted particles. `scripts/analytic_quadratic.jl` runs the analytic quadratic model, its full-order
  closed form and its ODE and PODE reduced forms. They are ported from the `bump_on_tail_*.jl`
  scripts and the `reference/` files of ReducedBasisMethods onto the marches, the diagnostics and
  the plots of this package. Every script writes into `--runs-dir` and `--results-dir`, and names
  each figure `item<NN>_<what>.png` after the item of the paper it draws, on the axes of the
  reference notebooks. The projections script also writes the particle basis and the reduced data
  `X` and `A` that the HNN reads. Every HDF5 file records the commit and the command that wrote
  it as attributes, and a bump-on-tail file also the seed. `scripts/bump_on_tail_common.jl` holds
  what the four bump-on-tail scripts share, and `test/integration/bump_on_tail_scripts.jl` tests
  it; `test/integration/bump_on_tail_runs.jl`, in the group `slow`, runs the three bump-on-tail
  scripts and the analytic model at small options. `scripts/Project.toml` gains
  `ReducedComplexityModeling`, at the package's bound `0.5`. The scripts call `spectrumplot` of
  `ReducedBasisMethods`, which every release in the package's bound `0.4` has.

- **The docs page "The paper's examples"** says how to run the five scripts, and gives the paper
  command, the measured run time and the files of each.

- **Check 3 runs to `T = 100` in the suite and to `T = 1000` from `scripts/analytic_quadratic.jl`.**
  The testset `"Reduced Equations at T = 100"`, in the group `slow`, marches the PODE form, the
  `k = 20` cotangent basis by `LobattoIIIAIIIB(2)` at `dt = 1e-4`, over `10^3` chunks of `1000`
  steps and asserts `maximum(drift) ≤ (β·dt)² = 3.6e-7`. It marched the ODE form too and ran to
  `T = 1000`, which took 68 of the 112 minutes of the job `Julia min - windows-latest` on PR #96
  (run 41179693589), whose limit is 120. The script now writes the root attributes `bound`,
  `max_drift` and `first_exceedance` on `analytic_ode_k5.h5` and `analytic_pode_k20.h5`, and prints
  them. Its paper command, `--final-time 1000 --saves 10000`, run on commit `be39325` on 2026-10-31
  with Julia 1.13.3 on Darwin arm64 (Apple M4 Max), gives `bound = 3.6e-7` for both forms, the ODE
  form (POD `k = 5`, `RK4`) `max_drift = 2.7063313296820274e-11` and `first_exceedance = 0`, and
  the PODE form `max_drift = 8.978506025917879e-8` and `first_exceedance = 0`. The script trains on
  the closed form at its own 10 001 saved times, spacing `0.1`, and the test on 1001 times, spacing
  `1`, so these are the paper command's numbers and not the test's to the bit.

- **`Picard_iterate_over_particles!`, the metriplectic Picard step in place.**
  `Picard_iterate_over_particles!(v, vn, vn_minus_one, dv_history, ti, t, Δt, abstol, reltol, mlb;
  maxiters = 1000)` writes the step into `v` and takes no `dv`, `m` or `β`. The model keeps the
  `SimpleSolvers` solver and its state for the element type of its particles, built with the model
  at the particle count, and the model's cache of that element type keeps the work arrays of the
  vector field. A call with other `abstol`, `reltol` or `maxiters` than those the solver was built
  for builds it again, a step whose `vn` has another element type than the particles builds a new
  solver on each call, and a step with `ti ≥ 4` allocates nothing when `vn` has the element type of
  the model's particles, where each step built a new solver and its unused `N×N` Jacobian cache
  before (K40, removed). The first three steps still allocate
  for their midpoint guess. `Picard_iterate_over_particles` keeps its signature and calls the
  in-place form, so its results are the same to the last bit:
  `test/data/lenard_bernstein_metriplectic_trajectory.txt` holds ten steps of a seeded state from
  commit `472dd42`, and both forms reproduce them with `==` on macOS aarch64, where the file was
  written, and to `rtol = 64·eps(Float64)` elsewhere. The in-place form throws a `DimensionMismatch`
  for a `v`, a `vn` or a `vn_minus_one` of another length than the model's particle count, and for
  a `dv_history` that is not `(N, 2)`, with `N` the particle count. It throws an `ArgumentError` for
  a `v` that shares memory with `vn`, `vn_minus_one` or `dv_history`, and for a `vn` or a
  `vn_minus_one` that shares memory with `dv_history`. The wrapper therefore refuses a `vn` or a
  `vn_minus_one` of another length than the model's particle count too. K42 now names only the
  wrapper, and K47 counts thirteen lines in four files, because the six commented-out lines of the
  old step body are gone.

### Bug Fixes

- **`initialize!(dist, BumpOnTail())` draws.** Both `draw!` methods of the bump-on-tail example
  bound `ParticleDistribution{1, 1}`, whose first parameter is the element type, so no 1D1V
  distribution matched them and the call raised a `MethodError`. They now take a 1D1V
  `ParticleDistribution` of any element type. The importance-sampling method also read a field
  `a` that `BumpOnTail` does not have; it reads the tail share `α`.

- **The grid right-hand side allocates nothing (K53).** `_grid_1d1v_rhs!`, the vector field of the
  grid 1D1V Vlasov–Poisson model, runs at every stage and every Newton iteration of a step, and it
  allocated the density, `ϕ`, `h` and the bracket's product on each call. The `params` that
  `GeometricIntegrator` builds hold buffers for them, in the element type of the distribution, and
  a state of that type writes into them through GeometricBrackets' `poisson_apply!`: a call
  allocates 0 bytes, where the allocating form allocated 8 800 bytes (Float32, 16 × 64 grid) up to
  264 096 bytes (Float64, 64 × 256). Per call it is as fast as or faster than the allocating form
  at every measured size, and the results do not change. A state of another element type, such as
  the `ForwardDiff.Dual` numbers of the implicit midpoint's Jacobian, keeps the allocating method.
  Two concurrent calls on one `params` are unsafe, because they share the buffers. The bound of
  `GeometricBrackets` rises from `0.2` to `0.3`, the first release with `poisson_apply!`.

- **`Picard_iterate_over_particles` checks its arguments before it writes.** Before this release, a
  history with one row too few let the solve run, wrote column 2 of the history, and only then
  threw a bare broadcast `DimensionMismatch`. A `vn` that was a view of the history, such as
  `view(dv_history, :, 1)`, was overwritten by the history update after the solve, so the caller's
  `vn` changed by up to 6.63 in a 64-particle probe. In the Hermite guess of `ti ≥ 4`, a one-entry
  `vn_minus_one` was broadcast over the particles, the solve wrote `v`, and in a 64-particle probe
  the call then threw an unrelated `DomainError`; lengths 63 and 65 threw a bare broadcast
  `DimensionMismatch`. A `vn_minus_one` that was a view of the history, such as
  `view(dv_history, :, 2)`, was overwritten by the history update at `ti = 1`, so the caller's
  `vn_minus_one` changed by 2.0 in a 64-particle probe. The wrapper now throws a `DimensionMismatch`
  for a `dv_history` that is not `(N, 2)`, and for a `vn_minus_one` of another length than the
  model's particle count, with `N` the model's particle count, and an `ArgumentError` for a `vn` or
  a `vn_minus_one` that shares memory with `dv_history`, all before it writes to any argument.

### Breaking Changes

- **`[compat] julia` rises from 1.11 to 1.12.** Pkg no longer resolves the package on Julia 1.11.
  Nothing in `src/` changes.

- **`[compat] ReducedComplexityModeling` is `"0.5"`, from `"0.1.1"`, and `ReducedBasisMethods` is
  `"0.4"`, from `"0.2.1"`.** ReducedBasisMethods 0.4 is the first release that binds
  ReducedComplexityModeling 0.5; 0.2 and 0.3 bind 0.1, so the two bounds move together, and
  `scripts/Project.toml` takes ReducedComplexityModeling `0.5` with them. `ParameterSpace`,
  `h5save` and `h5load` are ReducedComplexityModeling's, so a caller moves across its breaking
  releases 0.2 to 0.5. An HDF5 round trip of a `ReducedBasis` or of `Snapshots` keeps the written
  order of its parameters now, where ReducedComplexityModeling 0.1 returned them in alphabetical
  order.

- **`[compat] SimpleSolvers` rises from `0.14.1` to `0.14.2`.** SimpleSolvers 0.14.1 bounds
  `NeuralNetworkParameters` to 0.4, and ReducedComplexityModeling 0.5 requires 0.5, so the old
  floors could not be resolved together. GeometricOptimizers 0.9, the HNN extension's, requires
  SimpleSolvers 0.14.2 as well.

### Internal

- **`KNOWN_ISSUES.md` K44 leaves the file.** The reference generators stay in `test/helpers/`.
  Three of them write `test/data/` — `generate_collision_reference.jl`,
  `generate_metriplectic_reference.jl` and `generate_metriplectic_trajectory.jl`, each through
  a path relative to its own file, `joinpath(@__DIR__, "..", "data", …)`. A generator moved to
  `scripts/` would resolve that to `<root>/data/<name>`, where no `data/` directory exists, so
  `open` would throw a `SystemError`: it creates no intermediate directory.
  `vlasov_ode_reference.jl` writes nothing and stays because
  `test/particles/reduced_equations.jl` includes it.

- **`KNOWN_ISSUES.md` states which number a new entry takes.** It is the number after the highest
  `K<n>` that the file, `CHANGELOG.md` or their history has ever named, so the ID of an entry that
  has left the file, or that lived only on a merged branch, is not given again.
  "The next `K<n>`" read as the next number after the highest one left in the file, and that
  reuses a deleted ID.

## [0.3.0] — 2026-10-16

### Bug Fixes

- **The `Vlasov–Poisson` model is constructible and its splitting method matches it.** The first
  parameter of `DistributionFunction` is the element type, but seven declarations bound their
  `XD` to it, so `VlasovPoisson{Float64, 1, …}` carried a `Float64` dimension and
  `SplittingMethod(::VlasovPoisson)` matched no model. The declarations now bind
  `DistributionFunction{<:Any, XD, VD}`, so the two named parameters are the dimensions.

- **`run!(::SplittingMethod, h5file)` works with GeometricIntegrators 0.18.** One solution step
  is advanced step by step through the time span, with each step written to HDF5 as computed
  and the final state copied into the model's particles. The output is bit-identical to
  `GeometricIntegrators.integrate` on the same problem for Float64 and Float32.

- **`run!(::GeometricIntegrator, h5file)` works with GeometricIntegrators 0.18.** It read the
  problem's `tspan` field, which GeometricEquations no longer has, so every Lenard–Bernstein run
  threw `FieldError` before its first step. It now advances one solution step as
  `run!(::SplittingMethod, h5file)` does, and its output is bit-identical to
  `GeometricIntegrators.integrate` on the same problem.

- **The grid collision operator works.** `CollisionTensor` is indexed by three
  `CartesianIndex`es or three `Int`s, returning `ct.f(I, J, K)`; `ReducedCollisionTensor`
  projects it. `QuadraticCollisions(nx, nv, hx, hv, v)` has `v` as a type parameter
  `VT <: AbstractVector{DT}`, throws `DimensionMismatch` if `length(v) ≠ nv`, and
  `ArgumentError` if `nv < 3` (at `nv = 2` the neighbours coincide). Both
  `_get_MC̃_cubic` and `_get_MC̃_quadratic` take the velocity grid as their fifth argument,
  throw `DimensionMismatch` when `length(v) ≠ size(ci, 2)`, and return `eltype(V)` tensors.
  The `v`-stencils wrap the bounded grid periodically; conservation holds only while `f`
  vanishes at the ends, recorded in `KNOWN_ISSUES.md` (K27). K2 is closed.

- **`d(x, v)` with scalar arguments threw `MethodError`.** The `Vararg` call operator of every
  `DistributionFunction` used `view` on a `Tuple`, which does not support it. It now builds the two
  `SVector`s directly.

- **Twelve dependencies the package never used are gone.** `AdaptiveRejectionSampling`,
  `LaTeXStrings`, `LinearSolve`, `NaNMath`, `OffsetArrays`, `Plots`, `QuadratureRules`,
  `SciMLBase`, `SimpleSolvers`, `StatsBase`, `StatsPlots` and `Trapz` all left `[deps]`, together
  with their `[compat]` entries and their `using` lines in `src/VlasovMethods.jl`. Installing
  VlasovMethods no longer pulls in the `Plots` or `LinearSolve` stacks. `SimpleSolvers`, one of
  the twelve, returns below for the metriplectic solve; the other eleven stay gone. Nothing about
  the API changes: `using X` inside a module re-exports nothing, and all 68 exported names still
  resolve.

  Each removal was checked by matching the package's exported names against the *parsed syntax
  tree* of the 37 files `src/VlasovMethods.jl` includes, so that a name inside a comment, a string
  or a docstring does not count as a use. Four of the twelve are invisible to a plain `grep`:
  every name `LinearSolve` exports is also exported by `NonlinearSolve` (which this release also
  removes); `QuadratureRules` matched only a local variable named `weights`; `SimpleSolvers` and
  `NaNMath` were reached only from commented-out code, and `SimpleSolvers` returns below.
  `LaTeXStrings`, `OffsetArrays` and `Plots` were reached only from `src/electric_field.jl` and
  `src/visualisation.jl`, which no `include` names, and both were deleted; the
  `src/electric_field.jl` the reduced field brings back uses none of them. `SciMLBase` was
  reached only from `scripts/`, and the `scripts/` rewrite recorded under *Removals* leaves no
  caller of it anywhere.

- **`AdaptiveRejectionSampling` no longer constrains a runtime install.** It was used in one line
  of `test/projections/projections.jl` and nowhere in `src/`, yet it sat in `[deps]` rather than
  `[extras]` — so a test-only package capped `ForwardDiff` at `0.10` for everyone installing
  VlasovMethods. It was the only dependency doing so. The density the test samples is exactly the
  normal with mean `1/2` and standard deviation `1/(2π)`, so the test now draws from that
  directly, redrawing the ~0.17% of points that fall outside the domain.

- **`[compat]` bounds two dependencies that had none.** `SimpleSolvers` is bounded
  (CompatHelper #32, restored for the metriplectic solve; the bound is in the `GeometricBase` 0.15
  entry below), and `Parameters` gains `0.13`. Those are
  all that is left of the eight open CompatHelper requests: the other six still name a dependency
  this release deletes — `NonlinearSolve = "4"` is now among them, since `NonlinearSolve` itself
  leaves `[deps]`, while `SimpleSolvers` leaves that group because it returns.

- **The package is on the `GeometricBase` 0.15 line, and only that line (K35).** `[compat]` reads
  `GeometricBrackets = "0.2"`, `GeometricEquations = "0.21.5"`, `GeometricIntegrators = "0.18.6"`,
  `SimpleSolvers = "0.14.1"`, `SimpleSplines = "0.3.1"` and `PoissonSolvers = "0.6"`. Each floor is
  the first release on that line, so no floor admits a version that cannot resolve.
  `SimpleSolvers` 0.14.1 requires `GeometricBase` 0.15, and `GeometricBrackets` 0.2 requires
  `SimpleSolvers` 0.14.1 and `SimpleSplines` 0.3, which no `PoissonSolvers` below 0.6 admits.
  `PoissonSolvers` 0.6 and `SimpleSplines` 0.3 had to widen together for the same reason.
  The 0.14 line, the `SimpleSolvers` 0.13 bound and the `PoissonSolvers` 0.4 and 0.5 releases are
  no longer admitted, so an environment that pins one of them does not resolve. `src/` is
  unchanged. `GeometricBrackets` 0.2 rejects `PoissonTensor(Float32, nx, nv, ::Arakawa{Float64})`
  on construction with an `ArgumentError`, where 0.1 built it and `ReducedTensor` then threw a
  `MethodError`.

- **`ReducedBasisMethods` joins `[deps]`, and `MultiIndexArrays` rises to `0.1.2`.** `0.1.2` is
  the first `MultiIndexArrays` that exports `GridTensor`, and the floor of `ReducedBasisMethods`
  0.2, whose `ReducedTensor` the three reduced grid tensors are built on. An environment that
  pins `MultiIndexArrays` 0.1.1 does not resolve.

- **Three more `[compat]` floors name a version that resolves.** `HDF5` `0.16, 0.17` becomes
  `0.17`: no `HDF5` 0.16 release resolves with `GeometricIntegrators` 0.18.6, which requires
  `GeometricSolutions` 0.6.6, and that release admits only `HDF5` 0.17. `ParticleMethods` `0.1`
  becomes `0.1.2`, the first release that admits `HDF5` 0.17. `StaticArrays` `1` becomes
  `1.9.13`, the floor that `SimpleSolvers` 0.14.1 sets. The old floors could not be selected, so
  the advisory `Downgrade` job stopped at `forcedeps`. No environment that resolved before stops
  resolving.

- **`PoissonSolvers` and `StaticArrays` are no longer listed in `[extras]`.** Both are genuine
  `src/` dependencies and are already in `[deps]`, where the test environment picks them up.
  Listing them in both places was redundant. `OffsetArrays` was listed there too; the one file
  that used it, `test/electric_field.jl`, was in the `broken` group of `test/runtests.jl`
  (issue #56), and the `src/electric_field.jl` it covered was commented out. The
  `test/electric_field.jl` and `test/particles/electric_field.jl` the reduced field brings back
  run in `core` and use no `OffsetArrays`, which stays a test-only dependency of
  `test/Project.toml`.

- **`scripts/` has its own environment.** Until now the scripts ran against the package's own
  `Project.toml`, which forced every package they use — `GLMakie` among them — to be a dependency
  of `VlasovMethods` itself, or to go undeclared. `GLMakie` was in fact undeclared and worked only
  by accident of the shared environment.

- `src/models/vlasov_poisson.jl` is now Unicode NFC-normalised. It stored `ż` as a base letter
  plus a combining mark, inherited from macOS rather than chosen. Nothing about the compiled code
  changes — Julia's parser normalises identifiers to NFC — but a `grep` pattern or an editor
  search typed in NFC now matches, where before it silently matched nothing. The file is
  byte-equal to the NFC normalisation of its predecessor, and no string literal was affected.
  `scripts/charged_particles.jl`, normalised in the same commit, is deleted by the `scripts/`
  rewrite recorded under *Removals*.

- **`run!(::SplittingMethod, h5file)` steps over the solution's own time grid.** It sized its
  dataset by `ntime(equation)`, which rounds the step count up, and then indexed the range
  `t₀:Δt:t₁`, which is one point shorter whenever the span is not an exact multiple of the step: a
  run raised a `BoundsError` before its last step. The dataset now holds one slice per time point
  of the solution, the loop takes the steps `GeometricIntegrators.integrate` takes, and the HDF5
  chunk of a run with no particles no longer has a zero extent.

- **Charge deposition works, and its test is enabled.** `projection!` onto a
  `PoissonSolvers.Potential` could never run: it was written against `BSplineKit`, dispatching on
  `Potential{<:PeriodicBSplineBasis}` where that name resolves to the `SimpleSplines` type, so no
  method matched and `update_potential!` raised a `MethodError`. Its body then called
  `Splines.PeriodicVector`, which this package imports nowhere. It is rewritten on
  `evaluate_all!` and `basis_index`, and `test/projections/projections.jl` — commented out in
  `runtests.jl` — now runs, reconstructing a sampled density from a million particles to
  `2.5e-2`.

  It deposits onto a **`DirichletBasisSpline` as well as a periodic one**, which are the two
  bases a `Potential` can hold. A recombined basis has no single block width, so `evaluate_all!`
  pads the tail of its buffer with zeros whose indices run past `nbasis`; depositing those
  entries raised a `BoundsError` on every Dirichlet potential. They are skipped, and a nonzero
  value at an out-of-range index still raises rather than being folded onto another basis
  function.

  The test now also checks the deposited load vector **against a direct sum over every basis
  function**, on both bases. The density reconstruction alone does not: it still passes when the
  periodic contributions are dropped instead of wrapped, which is the property the rewrite
  turns on. The sampler is seeded, so the tolerance is asserted against the same draw each run.

- **`scripts/vlasov_poisson.jl` runs again.** It called `PeriodicBasisBSplineKit`, which
  `PoissonSolvers` 0.4 removed, so the only driver of the deposition path above was an
  `UndefVarError`. It calls `PeriodicBasisSpline` instead. The argument is unchanged in meaning:
  the old `nknots` built `nknots + 1` breakpoints, so it already counted cells.

- **`scripts/lenard_bernstein.jl` runs again.** Its animation loop called `Derivative(1) * f`, a
  name `BSplineKit` used to supply; that package is gone from `[deps]`, so the first frame raised
  an `UndefVarError`. `projection` returns a `SimpleSplines.Spline`, so the line becomes
  `derivative(f)` — what the `BSplineKit` migration table under *Breaking Changes* prescribes. A
  `SplineDerivative` is callable, so the `df.(xgrid)` plot below it is unchanged, and the
  rewritten driver draws the same curve.
  Six further `Derivative(1) * f` lines stood in `scripts/`, all commented out; they go with the
  scripts the rewrite recorded under *Removals* deletes or merges.

- **The package loads again.** `using VlasovMethods` failed outright. The `[compat]` bound
  `GeometricIntegrators = "0.16"` held `RungeKutta` at `0.5`, which still depends on
  `GenericLinearAlgebra`; that package guards a `LinearAlgebra.eigencopy_oftype` definition with
  `VERSION < v"1.14"`, and Julia 1.13's `LinearAlgebra` now defines the same method itself, so
  precompilation aborted with `Method overwriting is not permitted`. Raising the bound to
  `"0.18"` drops `GenericLinearAlgebra` from the manifest entirely. `QuadratureRules` had to move
  from `"0.1.6"` to `"0.2"` in the same step, because `GeometricIntegrators 0.18` requires it;
  the accessors this package uses (`nnodes`, and the `nodes`/`weights` fields of a
  `QuadratureRule`) are unchanged across that bump.

- **`GeometricIntegrators` was never imported under its own name.** `import A.B` binds only `B`,
  so the four call sites qualifying `GeometricIntegrators.integrate!`,
  `GeometricIntegrators.GeometricIntegrator`, `RK438` and `Strang` threw `UndefVarError` at run
  time — and had done so under `0.16` as well. `src/methods/splitting.jl`,
  `src/models/lenard_bernstein.jl` and `src/models/vlasov_poisson.jl` were affected.

- **The `Extrapolators` submodule no longer exists.** It was flattened in
  `GeometricIntegrators 0.18`; `extrapolate!`, `HermiteExtrapolation` and `MidpointExtrapolation`
  now live in `GeometricIntegratorsBase` and are re-exported at the top level. `import
  GeometricIntegrators.Extrapolators` only *warns* rather than failing, leaving the binding
  undefined, so the eight call sites in `src/methods/Landau_solver.jl` and
  `src/models/lenard_bernstein_metriplectic.jl` threw when reached. The three names are now
  imported directly and the call sites de-qualified, rather than a `const Extrapolators = …`
  shim standing in for a module that is gone.

- **The rescaled conservative operator ran the unrescaled one.** Its `GeometricIntegrator`
  wired `CLB_rhs_GI!` instead of `RCLB_rhs_GI!`, so every run built as a
  `RescaledConservativeLenardBernstein` integrated the conservative operator and
  `RCLB_rhs_GI!` was dead code. Separately, `compute_coefficients_rclb` — which solves the
  correct system for the rescaled parametrisation — was **never called**: both consumers called
  the conservative `compute_coefficients`, feeding `A₁ = -u/σ²` in as the multiplier of
  `f_s'/f_s` and `A₂ = 1/σ²` as the constant drift. Fixing the wiring is what first exposed the
  second bug.

- **The two-dimensional deposition allocated once per particle.** The `ntuple(D) do k … end`
  closure in `_deposit!` assigned to two captured variables, which boxes them, and `_deposit!`
  is the innermost loop of every step. Measured on a cubic 11-knot tensor basis with 2000
  particles: **4 611 136 B → 3 136 B** per `projection`, i.e. 2305.6 B per particle down to a
  fixed cost, with `Core.Box` gone from `code_typed`. The one-dimensional method never had it.
  At five Picard iterations per step this was ~23 MB of garbage per time step, in the operator
  whose headline is the `O(M²n⁴n_q⁴) → O(Q²)` speedup.

- **`compute_L!` evaluated the collision kernel five times per node pair.** `A^{cd}` swept the
  `Q²` pairs once, and the chunked second term then swept them again for each of the four
  `(c,d)` components, discarding three of the four returned components every time. The kernel is
  symmetric in `(c,d)`, so one sweep per block now fills all three independent components and
  the four contractions are read off them. The chained triple product also allocated a fresh
  `M×Q` and `M×M` per pass and now goes through scratch. Measured at `M = 169`, `Q = 2500`:
  **0.89 s → 0.49 s** and **49.1 MB → 36.4 MB** per call, agreeing with the previous result to
  `7.8e-16` relative.

- **`scripts/landau_profile.jl` still did not run.** It bound `pdist` and then read `dist` at
  five sites and `sdist2` at one, neither of which exists; the previous round fixed the same
  class of defect two lines away and stopped short. It also passed `rhs_full[:, :, 2]` — a copy
  — as an output argument, so the write went to a temporary. `scripts/landau_new.jl` had the
  same copy-as-output bug; `landau_newer.jl` already used a `view`. All three are among the
  scripts the rewrite recorded under *Removals* deletes or merges, because no Landau run gets
  past the positivity guard of K8; the profile they carried is
  `scripts/profile.jl --model landau`.

- **Three docstrings documented nothing.** A comment between a docstring and its definition
  detaches it silently: `(ent::CollisionEntropy)()` produced *no* documentation at all, and the
  file-header blocks in `projections/density.jl` and `projections/distribution.jl` were inert.
  The comment now precedes the docstring, and the two header blocks — which are file prose, not
  API documentation — are `#` comments. The formatter, the linters, the load test and
  `Pkg.test()` all pass either way; only a docs build sees it.

- **The metriplectic bracket docstring sat on the wrong function.** It states `L_αβ` with the
  `ε_h - u_h²` denominator, which is what `rhs_downstairs_factor!` computes and is the live
  right-hand side; it was attached to `rhs!`, which computes `(ε-u²)·L·dS`. The two now carry
  their own. `compute_dS!`'s signature line was also missing the `pdist` argument it takes.

- **`evaluate` is exported again.** It was dropped from the export list while two
  `evaluate(::SplineDistribution, …)` methods remained and the `DistributionFunction` call
  operators dispatched through it, so `using VlasovMethods` could not reach it and the tests
  worked around it with `VM.evaluate`.

- **`similar(AT, ::SplineDistribution)` no longer collides with five dependencies.** `AT` is an
  element type, but the argument was untyped, so the method read as
  `similar(::Any, ::SplineDistribution)` and was ambiguous against the `similar` of
  `Polynomials` (three), `RecursiveArrayTools` and `GeometricEquations`. Annotating it `::Type`
  takes `detect_ambiguities` from **7 to 1**; the remaining pair is the pre-existing
  `DistributionFunction{T,XD,0}` / `{T,0,VD}` call operators.

- **`compute_K!` left stale entries in `K1` and `K2`.** They are cache arrays written only at
  the entries the current particle positions overlap, and were never cleared, so as particles
  moved between cells `K` accumulated nonzeros from earlier Picard iterations and earlier time
  steps — making `K⁺` the pseudo-inverse of something that was not `w_α ∇φ_k(v_α)`.

- **Out-of-support contributions were aliased rather than discarded.** `compute_K!` bounded only
  the flattened index: with `i = 0`, `j = 3`, `M = 10` the flat index `(j-1)M + i = 20` passes a
  `1 ≤ k ≤ M²` test and decodes to `(10, 2)`, an unrelated basis function. Every index is now
  bounded per component. Relatedly, the gradient evaluators applied no periodic wrap while the
  spline evaluator did, so `∇φ` and the `f_s` it differentiates disagreed on a periodic basis;
  both now go through `basis_index`.

- **The particle weights were missing from the conservative coefficient system** — from all five
  sums. Both the matrix and the right-hand side are linear in `w`, so uniform weights cancel and
  no published result changes; what the old code annihilated was `Σ_α v̇_α` rather than
  `Σ_α w_α v̇_α`, so with the non-uniform weights `examples/bumpontail.jl` produces, momentum and
  energy conservation were lost. The same `w[1]`-for-every-particle substitution was in both
  metriplectic right-hand sides and in its entropy derivative, where `1/length(dS)` stood in for
  `w_α`.

- **The collision frequency `ν` was ignored** by both `Landau` and
  `MetriplecticLenardBernstein`: declared, stored, and never read, so any value other than the
  default `1.0` was silently discarded.

- **The plain Lenard-Bernstein operator was missing its `1/f_s`.** It computed
  `-ν (f_s' + v f_s)`, the collisional *flux*, where the advection coefficient is
  `-ν (f_s'/f_s + v)`. The equilibrium condition is the same either way, so the fixed point was
  right, but particles in the tails were barely advected.

- **Several names used in live code were never imported**, so the lines threw when reached:
  `Spline` in `projections/density.jl` (the derivative branch), `Splines.PeriodicVector` in
  `projections/potential.jl`, and `Integrals` / `HCubatureJL` in `project_Maxwellian`. The first
  and third are gone with the rewrite.

- **`CollisionEntropy` had no working implementation.** It routed a one-dimensional spline
  through a hard-coded two-dimensional quadrature that built `SVector{2}` sample points, took
  `nquad` as a positional default rather than a keyword, and applied `log` unguarded — which is
  why the entropy computations of `scripts/` were commented out. They leave with the rewrite
  recorded under *Removals*: two in scripts it deletes or merges, the third in
  `lenard_bernstein_metriplectic.jl`, which survives with the comments dropped. It now integrates
  on the basis's own Gauß-Legendre grid in any number of velocity dimensions.

- **`Cache`/`CacheType` referred to fields that do not exist** on `Landau`,
  `ConservativeLenardBernstein`, `RescaledConservativeLenardBernstein` and
  `MetriplecticLenardBernstein` (`pdist`/`sdist`, or a stale `clb` binding). Unreachable, since
  `CacheDict`'s parent is a `*Cache` and the `*Cache` methods are the ones that fire, but wrong
  as written. `eltype` is now defined for `LandauCache`, which was making the seeded cache entry
  key on `Any` and never be retrieved.

- **`MLBCache`'s midpoint buffer was hard-coded `Float64`** (`zeros(N)` rather than
  `zeros(T, N)`), truncating the midpoint under any wider element type.

- **Two unused `where {DT}` type parameters** are removed from `GeometricIntegrator` methods,
  and `lenard_bernstein.jl` no longer qualifies names as `GeometricIntegrators.` where
  `Integrators.` is meant.

- **`test/distributions/particle_distribution.jl` was flaky at roughly one run in four**, and had
  never run at all: `runtests.jl` included only the two spline test files. It asserts
  `mean(v) ≈ centre atol = 3.5/sqrt(np)`, which for a uniform distribution of width `w` — whose
  standard error of the mean is `w/sqrt(12 np)` — is about 3.0 standard errors, so each of the
  120 such assertions fails about 0.27 % of the time. Adding the file to the suite also made it
  sensitive to whatever RNG state the preceding testset left behind. Now seeded, so the outcome
  is deterministic and a failure is reproducible. The 3σ tolerance is left as it was; widening it
  would change what the tests assert. Its `using Statistics` also had to be declared in
  `[targets] test`.

- **`Cthulhu` is removed from `[deps]`.** Dead — the only occurrences anywhere were two
  commented-out `using Cthulhu` lines in the two Landau drivers the `scripts/` rewrite recorded
  under *Removals* deletes — with no `[compat]` entry, and it does not precompile against Julia
  1.14-DEV (`UndefVarError: ConstPropResult not defined in Compiler`), which was the sole cause of
  the advisory `nightly` CI job failing. Anyone using it interactively should add it to their own
  environment rather than to the package's dependencies.

- **Three Landau driver scripts contained invalid Julia and had never been runnable.**
  `const landau_rhs!(v̇, v, params) = …` is a syntax error — `const` takes an assignment, not a
  function definition — in `scripts/landau_new.jl`, `scripts/landau_newer.jl` and
  `scripts/landau_profile.jl`. Found by JuliaFormatter, which could not reparse its own output.
  The last of the three also called `collisions_rhs!`, which does not exist. All three are
  among the scripts the rewrite recorded under *Removals* deletes or merges.

- **The singular-system guards missed the `±Inf` half of the failure.** Both
  `compute_coefficients` and `compute_coefficients_rclb` tested `isnan(A₁) || isnan(A₂)`, which a
  vanishing determinant produces only when the numerator vanishes with it. With a nonzero
  numerator the division gives `±Inf`, the guard passed, and non-finite coefficients reached the
  right-hand side. Both now test `isfinite`; the message is unchanged.

- **The Landau vector field did the `O(Q²)` kernel sum before checking that `f_s > 0`.**
  `compute_L!` and `compute_J!` are independent, but the positivity check on the quadrature grid
  lives in `compute_J!`, which ran second — so a non-positive projection threw only after the
  expensive work. They are simply swapped. No results change.

- **Four driver scripts still declared `using QuadGK`** after the dependency was dropped from
  `Project.toml`, so each failed at load in the package environment. Every `quadgk` call site in
  them is commented out, so the import goes rather than the dependency coming back. Three of the
  four — `lenard_bernstein_conservative.jl`, `lenard_bernstein_metriplectic_scaling.jl` and
  `rescaled_lenard_bernstein_conservative.jl` — are among the scripts the rewrite recorded under
  *Removals* deletes or merges; the fourth, `lenard_bernstein_metriplectic.jl`, survives, so
  nothing under `scripts/` names `QuadGK` at all.

### Breaking Changes

- **`[compat]` bounds narrowed since 0.2.1.** Eight entries lost versions between `v0.2.1` and
  this release, so a manifest that resolved against one of the dropped versions has to move with
  them:

  | entry | `v0.2.1` | `0.3.0` |
  |:--|:--|:--|
  | `julia` | `1.8` | `1.11` |
  | `GeometricEquations` | `0.14` | `0.21.5` |
  | `GeometricIntegrators` | `0.13` | `0.18.6` |
  | `HDF5` | `0.15, 0.16, 0.17` | `0.17` |
  | `ParticleMethods` | `0.1` | `0.1.2` |
  | `PoissonSolvers` | `0.3` | `0.6` |
  | `ProgressMeter` | `1` | `1.3.3` |
  | `StaticArrays` | `1` | `1.9.13` |

  The dependencies that left or joined `[deps]` in the same interval have no row. Every name
  `v0.2.1` exported was its own or `Base`'s, so a dependency that left takes no name from a
  caller, and one that joined had no `v0.2.1` range to narrow. The names this release
  re-exports from its dependencies, such as `Spline` and `derivative`, are new exports; the
  `SimpleSplines` bullet of this section maps the `BSplineKit` calls they replace.

- **The metriplectic Lenard–Bernstein solve uses `SimpleSolvers`.** The one live
  `NonlinearSolve.NLsolveJL` call, in `Picard_iterate_over_particles`, is replaced by an
  unaccelerated `SimpleSolvers` Picard solve, so `NonlinearSolve` leaves `[deps]` and
  `SimpleSolvers` joins it at `0.13`; the dead `NonlinearProblem` construction it had supplied a
  name for in `Picard_iterate_Landau_nls!` goes with it. The solve no longer prints `@time` or an
  NLsolve trace, and it throws when it does not converge or meets a non-finite residual, giving
  the residual norm and the iteration count. The convergence test is
  `((x_settled || f_settled) && ‖F‖₂ ≤ abstol + reltol·‖F₀‖) || ‖F‖₂ ≤ abstol`, where
  `NLsolve` used `maximum(abs, F)`; the drivers set `reltol = 1e-50`, so it reduces to
  `‖F‖₂ ≤ abstol`, which is stricter than the old test and can turn a former success into a throw.
  `Picard_iterate_over_particles` now returns the solved velocity vector instead of an NLsolve
  solution object, so a caller that read `.u` or `SciMLBase.successful_retcode` changes, and it
  gains a `maxiters` keyword (default `1000`). The parameters `dv`, `m`, and `β` stay in
  the signature and are ignored. A one-step pin test in
  `test/integration/metriplectic_solve.jl` matches the pre-change result to `3e-16·√N` in `Float64`
  and to `3e-16·√N · eps(Float32)/eps(Float64)` in `Float32`.

- **`[compat] julia` rises from 1.10 to 1.11.** `GeometricBrackets` requires 1.11.

- **`projection_energy` is now `projection_second_moment`.** It returns `Σ_α v_α² f_s(v_α)`, the
  particle-sampled second moment of the spline distribution, with no factor `½`, so it is not an
  energy. The name was not exported; a caller of `VlasovMethods.projection_energy` renames the
  call.

- **The spline machinery is now `SimpleSplines`, and `BSplineKit` is gone.** `src/splines/` —
  `SplineND`, `TwoDSpline`, `NDSpline`, and the `gauss_quad*` / `eval_bfd` /
  `evaluate_der_2d*` / `remap_unit_interval` / `unique_knots` helpers — is deleted, together
  with `src/test.jl`. What replaces it:

  | was | is |
  |:--|:--|
  | `SplineND`, `TwoDSpline` | `SimpleSplines.TensorProductBasis` + `Spline` |
  | one basis shared by every velocity axis | one basis **per axis**, with its own degree, mesh, domain and boundary condition |
  | `cholesky(kron(M, M))`, dense | `KroneckerMass` — `D` one-dimensional solves, never assembled |
  | `L2projection!`, `O(M² n^{2D} n_q^{2D})` | `l2_projection!`, `D` sparse contractions |
  | `BSplineKit.Spline`, `Derivative(1) * fs` | `Spline`, `derivative(fs)` |
  | `BSplineKit.galerkin_projection` + `ldiv!` | `l2_projection!` |
  | `BSplineKit.evaluate_all` + hand-written `mod1` | `evaluate_all!` + `basis_index` |

  `SimpleSplines` is bounded by `[compat] SimpleSplines = "0.3.1"` and resolves from the registry.
  `FastGaussQuadrature` and `QuadGK` are dropped, having no remaining caller.

- **`SplineDistribution`'s boundary-condition argument takes a type, and its default changed
  to `Free()`.** It was `:Dirichlet`. The old symbols `:Dirichlet`, `:Periodic` and `:nothing`
  still map to the spaces they used to select, so no script changes meaning, but an
  unrecognised symbol now throws instead of falling through an `else` branch to the
  unconstrained basis — which is how `:nothing` used to work, and how a typo used to work too.

- **The lumped-mass option is removed.** The trailing `compute_mass_galerkin::Bool` is gone.
  Passing `false` assembled the mass matrix with a trapezoidal rule, which for a B-spline basis
  samples at the knots where `φ_i(t_j) = δ_ij` and therefore produced a **diagonal** matrix —
  a mass-lumped discretisation, not a cheaper assembly of `M_ij = ∫ φ_i φ_j`, and not the one
  either manuscript describes. Three Landau drivers passed `false`. A call with the old flag is
  now a `MethodError` rather than being silently reinterpreted.

- **`VlasovPoisson`'s field follows the particles.** `v_acceleration!` and `s_acceleration!`
  deposited from `model.distribution`, which the integrator never writes, so every step evaluated
  the field of the initial state. They now deposit from the state `z` the integrator steps, so
  the field is that of the current positions. This changes the numerical result of any
  `VlasovPoisson` run.

- **The `Vlasov–Poisson` right-hand side allocates nothing, and adds a buffer type parameter.**
  `VlasovPoisson` now holds a `work` buffer of the width the deposit and the field evaluation
  need — a fifth type parameter `WT` — and the two fields fill it rather than allocating one per
  evaluation. Code that names `VlasovPoisson` with four type parameters stops matching. Because
  the constructor allocates that buffer, a **particle or spline** `VlasovPoisson` on a basis with no
  deposit-buffer method — an `FFTWBasis` or a `FiniteDifferenceBasis` — now throws at
  construction rather than at the first deposit (K19).

- **A particle outside the velocity domain is an error.** It used to be dropped silently in the
  one-dimensional projection and printed a per-particle warning in the two-dimensional one.
  Either way its weight left `f_s` while the particle kept moving, and `f_s` then evaluated to
  zero at its position, making the `f_s'/f_s` of every operator here a division by zero — which
  the Lenard-Bernstein manuscript names in its closing remarks. `projection` now throws, naming
  the count and the domain.

- **A non-positive `f_s` is an error.** An `L²` projection of a particle distribution
  undershoots into negative values where the sampling is thin. The metriplectic operator wrote
  `1 + 0.5*log(f_s^2)`, which is `log|f_s|` — finite, plausible-looking, and wrong: `S = ∫ f
  log f` needs `f > 0`, and where `f_s < 0` the H-theorem *reverses*. The construction turned a
  visible failure into an invisible one and removed the only diagnostic that would have caught
  it. Positivity is now checked wherever a logarithm or a `1/f_s` appears. This is the problem
  the Landau manuscript's commented-out "Positivity-preserving projections" appendix names and
  does not solve.

- **`MaxwellianDistribution` is normalised for the dimension it is called in**,
  `(2π)^{-d/2} exp(-|v|²/2)`. It was fixed at `1/(2π)`, correct only for `d = 2`; in one
  dimension a "normalised" Maxwellian integrated to `0.399`.

- **`test/spline_tests.jl` and `test/spline_basis_tests.jl` are deleted.** They tested only
  `SplineND`, `TwoDSpline` and the helper functions, all of which are gone. The dimension counts
  they encoded per boundary condition — `n+p` clamped, `n` periodic, `n+p-2` Dirichlet — are now
  asserted in `SimpleSplines`, and `test/distributions/spline_distribution.jl` covers the
  `SplineDistribution` level.

- **The guards throw typed exceptions rather than `ErrorException`.** A caller can now tell the
  three failures apart, and a test can assert which one it got. A non-positive `f_s`, and a
  particle outside the velocity domain, throw `DomainError`; a singular conservative-coefficient
  system throws `ArgumentError`. Nine sites across `landau.jl`, `collision_entropy.jl`,
  `projections/distribution.jl` and the three Lenard-Bernstein models. Code catching
  `ErrorException` around a projection or an entropy needs updating.

- **`Picard_iterate_Landau_nls!` no longer takes `n`.** The argument was threaded into
  `params.n` and read by nothing once the quadrature became the basis's own; the drivers passed
  it under the comment `n = 1 # number of quadrature nodes`, which named a quantity that no
  longer exists. The parameter is dropped rather than kept as an ignored placeholder.

- **`SplineDistribution`'s `domain` argument is typed `Union{Tuple, AbstractVector}`.** Left
  untyped it made every array-valued domain an ambiguous call against the low-level
  `SplineDistribution(xdim, vdim, basis, quadrature, coefficients)` constructor, so
  `SplineDistribution(1, 1, 41, 4, [-10.0, 10.0])` was a `MethodError` — while the docstring
  reads `domain` only through `first`/`last`, so a vector or a range is the natural thing to
  pass. All three forms now work and agree.

- **`:Natural` is not a boundary-condition symbol.** It belonged to the deleted `SplineND` API,
  never to `SplineDistribution`, and the migration mapped it to `Free()` — while the exported
  `Natural` type means something else, giving a different spline space (23 basis functions and
  polynomial reproduction 3, against 21 and 1). Passing it is now an `ArgumentError`. Nothing
  in `src/`, `test/` or `scripts/` used it.

- **`PoissonSolvers` 0.4 is now required**, and 0.3 no longer works. Its spline backend moved to
  `SimpleSplines` as well, so one spline library is loaded where two were. Three consequences
  here: the potential's derivative is `ϕ(x, 1)` rather than `ϕ(x, Derivative(1))`, the basis and
  coefficients are reached with `basis(p)` and `coefficients(p)` rather than `p.basis` and
  `p.coefficients`, and `Potential`'s first type parameter is the solution type rather than the
  basis type.

- **`SimpleSplines` is resolved from the registry.** Its 0.1.0 is registered, so the `[sources]`
  table that existed to support it is gone. This removes the `[sources]` table that RegistryCI
  rejects outright. This covers the `CHANGELOG.md` half of `KNOWN_ISSUES.md` K33: the sentence
  that told a reader the dependencies without a `[compat]` bound were recorded in a section this
  file does not have is gone.

- **`DiffEqIntegrator` and the `DifferentialEquations` paths it ran are gone.**
  `DifferentialEquations` was never imported — `src/VlasovMethods.jl` carried its `import` only as
  a comment — so every `DifferentialEquations.ODEProblem` and `DifferentialEquations.solve` call
  threw `UndefVarError` when reached. The exported type, its `run!` and `run` methods, the file
  `src/methods/diffeq_integrator.jl`, and the `DiffEqIntegrator` constructor methods for
  `LenardBernstein{1,1}` and `ConservativeLenardBernstein{1,1}` are deleted. The one script that
  called the deleted constructor, `scripts/lenard_bernstein.jl`, was rewritten against
  `GeometricIntegrator` in the same release; the rewrite is recorded under *Removals*.

- **`run` is no longer exported.** Its only method was in the deleted `DiffEqIntegrator`
  type; removing the export keeps `names(VlasovMethods)` from listing `Base.run`. The one live
  call, in `scripts/lenard_bernstein.jl`, is rewritten under *Removals* onto `run!`.

- **`LB_rhs` is no longer exported.** It is still defined, so a caller that reached it through
  `using VlasovMethods` writes `VlasovMethods.LB_rhs`. The collision operators are built with
  `LB_rhs_GI!`, which the `GeometricIntegrator` constructors name, so no run changes.

- **`CLB_rhs` is no longer exported.** The same, under the conservative model:
  `VlasovMethods.CLB_rhs`, and `CLB_rhs_GI!` is what the constructors name.

- **The `Vlasov–Poisson` acceleration pair carries its model instead of reading it from the
  parameters.** `v_acceleration!` and `s_acceleration!` are replaced by the callable structs
  `VAcceleration` and `SAcceleration`, built as `VAcceleration(model)` and called with the
  arguments the two functions took; neither is exported, and `SplittingMethod` and
  `sodeensemble` are the two places in `src/` that build them. The model, its potential and its
  work buffer leave the problem's parameters, which now hold the weights of the deposit alone;
  `SplittingMethod` passes `(w = model.distribution.particles.w,)`, the live view of the
  model's own weights, so its trajectories are unchanged, bit-identical in `Float64` and
  `Float32`. `update_potential!(model, z)` becomes
  `update_potential!(model, z, w)`, whose `w` is the first row of weights the deposit reads. A
  caller of the two old names, or of the two-argument `update_potential!`, builds
  `VAcceleration(model)` or `SAcceleration(model)` once and calls it, or passes the weights.

### Removals

- **Four `scripts/` drivers are deleted and four merge into a survivor.** The fifteen `.jl` files
  the rewrite triaged under `scripts/` become eight — seven drivers and one shared parser — beside
  `electric_field_timing.jl`, which the rewrite does not touch. `bump_on_tail.jl` goes:
  `using VlasovMethods.BumpOnTail` names a module that does not exist and six further names are
  absent, and nothing under `scripts/` reproduces the bump-on-tail set. `charged_particles.jl`
  goes: it never loads VlasovMethods and calls an API `GeometricIntegrators` 0.18 lacks.
  `landau_new.jl` and `landau_newer.jl` go: each meets the positivity guard of K8 at step 0, so
  no Landau time-stepping run gets past it at any parameters tried, and one solver step lives in
  `test/integration/landau_solver.jl`. `landau_profile.jl` merges into `profile.jl` as its
  `--model landau`, `lenard_bernstein_conservative.jl` and
  `rescaled_lenard_bernstein_conservative.jl` merge into `lenard_bernstein.jl` as two further
  `--model` values, and `lenard_bernstein_metriplectic_scaling.jl` merges into
  `lenard_bernstein_metriplectic.jl` as its `--npart` list. This covers `KNOWN_ISSUES.md` K20,
  K22 and K37, which leave the file.

- **Five dead `src/` files are deleted.** `vlasov_poisson.jl`, `electric_field.jl` and
  `visualisation.jl` were included by nothing and had only commented-out `include` and `export`
  lines; `methods/lbm_solver.jl` was included by nothing and duplicated the `IM_rule!` method
  of `src/methods/Landau_solver.jl` (same signature and body); `examples/twostream.jl` was
  included but held a newline and no code. The `include` and `export` lines that named them in
  `src/VlasovMethods.jl` go with them, as does the commented `# include("sampling.jl")` and
  `# export` block at the end of that file. The `@safetestset "Electric Fields"` line of
  `test/runtests.jl` goes too, since its `test/electric_field.jl` covered the deleted source. A
  sixth orphan the audit named, `src/hdf5.jl`, exists in no commit, so there is nothing to
  delete. This covers `KNOWN_ISSUES.md` K12, which leaves the file.

- **`lorentz_force!` is deleted.** It evaluated the field through `update_potential!`, which
  projects `model.distribution` and never the stepped `z`, so the self-consistent field stayed
  frozen at its initial value. No code called it.

- **`CLB_rhs!` is deleted.** Deleting the `DiffEqIntegrator` constructors left it with no caller,
  and its body duplicated `CLB_rhs_GI!`; its sign-convention docstring moves to `CLB_rhs_GI!`.

- **The dead `IM_rule!` and the dead untyped `f!` are deleted.** Both lived in
  `src/methods/Landau_solver.jl` with no caller: `IM_rule!` is reached only from commented-out
  code, and the untyped `f!` duplicates the typed `f!` of the metriplectic model, which every
  call site dispatches to. This covers `KNOWN_ISSUES.md` K25 and K39, which leave the file.

- **`Distances` leaves `[deps]`.** Its one use, `euclidean` in the two residual prints of
  `Picard_iterate_Landau_nls!`, becomes `norm`, which `LinearAlgebra` already provides; the
  printed residual can differ in its last digits, from rounding alone.

### New Features

- **Every `scripts/` driver that writes files reads its output directories from the command line.**
  A new `scripts/options.jl` is the shared argument parser, and `electric_field_timing.jl` does
  not include it: that driver writes no file and takes one positional argument, the name of the
  path (`reduced` or `deim`) its cold subprocess times. The parser accepts `--runs-dir` and
  `--results-dir`, defaulting to `runs/` and `results/` at the repository root, creates both, and
  rejects an unrecognised flag, or one given twice, with an `ArgumentError` rather than ignoring
  it. Each such driver writes its numbers under `--runs-dir` and a figure or a report under
  `--results-dir`, so a run at the defaults adds no untracked path and `git status --short` is
  unchanged. The first run creates both ignored directories; `/runs` and `/results` join
  `.gitignore`. The surviving drivers draw with CairoMakie, and `scripts/Project.toml` loses
  `Plots`, `LaTeXStrings`, `GLMakie`, `ParticleMethods`, `SciMLBase`, `TerminalLoggers` and
  `Logging` and gains `CairoMakie`. Their defaults are measured to run: `lenard_bernstein.jl`
  steps `dt = 1e-2` over `t ∈ [0, 1]` for all three of its models, where the `dt = 0.1` over
  `t ∈ [0, 10]` it carried fails for every one of them, and `lenard_bernstein_metriplectic.jl`
  starts from the cell-centred lattice on `[-2, 2]`, which keeps the projection positive where
  the shifted double Maxwellian on `[-10, 10]` did not. `lenard_bernstein.jl` builds all three
  models on a `Free()` basis, the one condition under which the spline moments equal the
  particles'; `profile.jl --model clb` uses `Free()` too, and its default is one step
  (`--final-time 1e-3`, equal to `--tstep`).

- **The grid 1D1V Vlasov–Poisson scheme.** `VlasovPoisson(dist::GridDistribution,
  potential::Potential)` is a new method of the exported constructor, and
  `VlasovMethods.GeometricIntegrator(model, tspan, tstep, method = Integrators.ImplicitMidpoint())` steps it.
  The state is `vec(dist.values)` and the right-hand side is GeometricBrackets' Arakawa bracket,
  `∂ₜf = −[f, v²/2 + ϕ]` with `ϕ` the potential's own solve, so the scheme holds the discrete
  Hamiltonian `hx·hv·Σ f v²/2 + (hx/2)·Σ ρᵢϕᵢ` to round-off, and the `L¹` norm `hx·hv·Σ|f|` to
  round-off plus twice the negative mass it accumulates, over a horizon long enough to distinguish
  that from a passing tolerance.

  - The field is solved once, at construction, as the `nx × nx` Green's matrix of the potential's
    backend: `ϕ = Gρ` for `ρᵢ = hv·Σⱼ f_ij`, which makes `∂H/∂f = hx·hv·h` exact. `G` is stored in
    `VlasovPoisson`'s `work` field, whose type parameter `WT` widens from `AbstractVector` to
    `AbstractArray`: `work` is the deposit buffer of a particle distribution and this matrix for a
    grid one. Code that requires `work` to be a vector stops matching.
  - The constructor throws `ArgumentError` unless the potential's solver is a `PoissonSolverFFT`
    or a `PoissonSolverMatrixFree` — on a spline basis the density at the nodes is not a Galerkin
    load vector — `DimensionMismatch` unless it has one degree of freedom per `x`-node, and
    `ArgumentError` unless its `x`-grid is the distribution's.
  - The scheme is **periodic in `v`** as well as in `x`: Arakawa's stencil conserves at a boundary
    only where `h` is constant along it, and `h = v²/2 + ϕ(x)` is not constant on `v = v_b`. A run
    is valid only while `f` vanishes to `eps(T)·maximum(f)` at both ends of the velocity range,
    which the docstrings state and the tests assert at every step.
  - `run!(::GeometricIntegrator, h5file)` works with the grid model. At the end it writes the
    final state back into the model's distribution: the velocities of a `ParticleDistribution`,
    the row the particle models step, with the positions left as they were, and
    `distribution.values`, reshaped, for a `GridDistribution`.
  - `GeometricIntegrator` is a `VlasovMethod` and no longer a `ParticleMethod`, because it steps
    the grid model too. Code that dispatches on `ParticleMethod` no longer reaches it.
  - `scripts/landau_damping_1d1v.jl` fits the linear Landau damping rate to the peaks of
    `|E₁(t)|` and reproduces the textbook `γ ≈ −0.1533` for `k = 0.5` and unit thermal velocity.
    Run with `julia --project=scripts scripts/landau_damping_1d1v.jl`, it measures
    `γ = −0.153534` at `nv = 1025`, `dt = 0.05` and `γ = −0.153474` at `nv = 2049`, `dt = 0.025`
    — an error of `2.3e-4` and a change of `6.0e-5` — and requires the first to be within
    `1.5e-3` of `−0.1533` and the two to agree to `1.5e-4`. That second requirement is the reason
    the 1% tolerance is defensible, and `nv = 513` fails it.

- **`Snapshots`, the snapshots of a parameter sweep, and the two reduced bases built from them:
  `particle_basis` and `field_basis`.** `Snapshots{T}` holds the positions `X`, the velocities
  `V`, the accelerations `A`, the potential `Φ` and the diagnostics `W`, `K` and `M` of a sweep on
  a `ParameterSpace`, with `X`, `V` and `A` of size `(nd, nₚ, nₛ, nparam)` and `Φ` of size
  `(nh, nₛ, nparam)`. `Snapshots(particles, poisson, nₛ, paramspace)` reads `nd` and `nₚ` from the
  particle list and `nh` from its Poisson solver, so a 2D particle list gives `size(X, 1) == 2`,
  and throws an `ArgumentError` when `v` holds a different component count than `x`, which the
  cotangent lift cannot pair. `Snapshots(T, nd, nₚ, nh, nₛ, paramspace)` is the zero array form, for
  a sweep of a model with no potential; with `nh = 0` its empty `Φ` round-trips through HDF5.
  The inner constructor takes plain `Array`s only: a view or an `OffsetArray` throws a `MethodError`.
  `==` and `isequal` compare all eight fields, `paramspace` included, and `hash` agrees with them.
  `h5save` and `h5load` write and read the eight fields under the keys `X`, `V`, `A`, `Φ`, `W`,
  `K` and `M` and the `parameterspace` group, through the generics of `ReducedComplexityModeling`,
  which take an open file or a file path. `particle_basis(algorithm, s; parameters, …)` reduces the
  stacked `[Q; P]` of `X` and `V` with any of the five `ReducedBasisMethods` tags, and
  `field_basis(algorithm, s; …)` reduces `A` with one of the two POD tags. The tags are
  `ReducedBasisMethods`' and VlasovMethods does not export them. Both pass the `ParameterSpace` and
  every keyword other than `parameters`, such as `k` or `tolerance`, to `ReducedBasis` unchanged.
  A snapshot matrix is the `(nd·nₚ) × (nₛ·nparam)` reshape of the four-dimensional array, so its column
  `(p−1)·nₛ + t` is `vec(X[:, :, t, p])`. The reduced-basis driver this replaces reshaped to
  `(nₚ, nₛ·nparam)` — the layout of `nd == 1` alone — and fed the acceleration to the field basis
  through that same particle-shaped reshape, so a sweep of a particle list with two position
  components threw a `DimensionMismatch` instead of reducing. `Φ` loses the component axis of that
  layout, which held the filled first row and the zeros of `nd − 1` more.

- **`VlasovWithExternalFields`, one species in an external electromagnetic field, and the
  `ElectromagneticFields` extension that supplies such a field.** `VlasovWithExternalFields`,
  `UniformField`, `electric_field` and `magnetic_field` are exported.

  - The model holds the particles of a `ParticleDistribution{T,3,3}` and the field they move in,
    with `charge` and `mass` keywords stored as `T`; the two-argument constructor throws an
    `ArgumentError` unless the converted mass is positive. `SplittingMethod(model, tspan, tstep)`
    integrates it: a Strang splitting of the drift `ẋ = v` and the Lorentz force
    `v̇ = (q/m)(E + v×B)`, whose rotational step is a Boris push, so a constant magnetic field
    moves a particle along the analytic gyro-orbit to second order and holds `½|v|²` to round-off.
  - An empty particle list is a **vacuum run**, and the constructor does not warn: the model, its
    `SplittingMethod` and `run!` work over the empty state. This needs `GeometricBase` 0.15.1,
    whose `State` keeps the vector field of an empty state variable. `GeometricBase` is a direct
    dependency, and `[compat]` bounds it at `0.15.1`, so a `GeometricBase` 0.15.0 does not resolve
    with this version.
  - `electric_field(field, t, x)` and `magnetic_field(field, t, x)` are the field interface, each
    returning an `SVector{3,T}` in Cartesian components, with `T` the element type of `x`. There
    is no abstract field type and no fallback method: a field defines one method of each.
  - `UniformField(E, B)` is the field type this package ships, constant in `x` and in `t`.
  - `ElectromagneticFields` is a **weak** dependency. The extension
    `VlasovMethodsElectromagneticFieldsExt` adds the two methods for a `FieldFunctions` of a
    Cartesian chart, so `FieldFunctions(ABCEquilibrium(a, b, c))` goes to the constructor as it
    is. It is confined to Cartesian charts by dispatch — `B♯` and `E♯` are the Cartesian
    components only where the chart is the identity map — so a curvilinear field raises its
    `MethodError` at the first kick. `[compat]` is `ElectromagneticFields = "0.9.1"`, the first
    release on `GeometricBase` 0.15.

- **`GridDistribution`, the distribution function on a 1D1V phase-space grid.** It is a third
  `DistributionFunction{DT, 1, 1}`, beside `ParticleDistribution` and `SplineDistribution`, and
  both `GridDistribution` and `velocity_moments` are exported.

  - It stores node values as an `nx × nv` matrix, in the column-major order `_apply_∫dv!` uses.
    The `x`-grid is uniform and periodic, and its right end is not stored. The `v`-grid is
    uniform and bounded, and both of its ends are nodes.
  - `d(x, v)` interpolates bilinearly. It is periodic in `x` and zero outside the `v`-range.
  - `velocity_moments(dist)` returns `(density, momentum, energy)` at each `x`-node, by the
    rectangle rule of `_apply_∫dv!`, `_apply_∫vdv!` and `_apply_∫v²dv!`. `energy` is the
    kinetic energy density `½ ∫ v² f dv`. The rule is exact up to
    `hv/2` times the two end values, so it is first order where `f` does not vanish there.
  - `xdim` and `vdim`, unexported, are defined once on `DistributionFunction`, so the spline
    distribution answers them too. It had neither.

  The tests show that the grid moments of a Maxwellian that vanishes at the ends of the range
  agree with the spline's L²-projected moments to `1e-12`. Bilinear evaluation is second order.
  On a truncated range the moments are first order, and second order once the end term is
  subtracted.

- **The reduced phase-space tensors and velocity moments, imported from ReducedBasisMethods.**
  Every function and struct body is byte-identical to its source, except
  `src/gridbased/collisions.jl`, whose working form is described under *Bug Fixes* above, and
  the three tensors of `src/gridbased/reduced_tensors.jl`, described below.
  Three files arrive under `src/gridbased/`:

  - `reduced_tensors.jl` — `PotentialReducedTensor`, `VelocityReducedMatrix` and
    `FullyReducedTensor`, described in the entry for the three reduced grid tensors below.
  - `moments.jl` — `_apply_∫dv!` and its transpose, the weighted pair `_apply_∫dv_μ!` /
    `_apply_∫dvᵀ_μ!`, and the first and second moments `_apply_∫vdv!` and `_apply_∫v²dv!`.
    Only `_apply_∫dv!` is exported, as before. The Laplace and nullspace stencils that shared
    the source file went to `PoissonSolvers` instead.
  - `collisions.jl` — `CollisionTensor`, `QuadraticCollisions`, `ReducedCollisionTensor` and
    the two `_get_MC̃_*` assemblers. ReducedBasisMethods never included this file, so none of
    it was reachable there. It **is** included here, and these five names are defined in the
    module; the run-time faults they arrived with are repaired under *Bug Fixes* above.

  New dependencies: `GeometricBrackets`, for `Arakawa`, which `bracket_tensor` takes, and for
  `poisson_apply`; `MultiIndexArrays`, for `multiindex`, `linearindex`, `GridTensor`,
  `MultiIndexAxis` and `MultiIndexLazyArray`; and `ReducedBasisMethods`, for the `ReducedTensor`
  the three tensors are built on.

- **The three reduced grid tensors, on `ReducedBasisMethods`' `ReducedTensor`.** Each of
  `PotentialReducedTensor`, `VelocityReducedMatrix` and `FullyReducedTensor` takes a bracket's
  rank-3 `MultiIndexArrays.GridTensor`, which `bracket_tensor` builds from an `Arakawa` bracket on
  that bracket's own grid, and the projections of its first two indices.
  `PotentialReducedTensor` is the only struct, a lazy `AbstractArray{DT, 3}`;
  `VelocityReducedMatrix` and `FullyReducedTensor` are functions over it.
  `VelocityReducedMatrix` returns `view(prt, :, :, 1)` of the `PotentialReducedTensor` whose
  one-column projection carries `v²[n₂] / 2` for each node, and `FullyReducedTensor` returns the
  `PotentialReducedTensor` whose projection has the row `Px[n₁, :]` for each node. Each builds
  `ReducedTensor(tensor, Pi, Pj)` from them first, so that its checks on the rows, on the element
  types and on the stencil width run once and live in one place. With `R` that reduced tensor and
  `(n₁, n₂)` the grid indices of the node `K`, `PotentialReducedTensor[i, j, k] = ∑_K R[i, j, K]
  · Pk[K, k]` projects the third index onto a basis of `size(Pk, 2)` modes over the `N` nodes,
  `FullyReducedTensor[i, j, k] = ∑_K R[i, j, K] · Px[n₁, k]` projects it onto the potential's
  modes at the node's `x` index, and `VelocityReducedMatrix[i, j] = ∑_K R[i, j, K] · v²[n₂] / 2`
  contracts it with the squared velocity nodes. `Px` has one row per `x` node and not one per
  phase-space node, and `v²` has one entry per velocity node, so it holds the `nv` squared
  velocity nodes `dist.v .^ 2` and not one entry per grid node; the `½` is written `/ 2`, so a
  `Float32` sum stays in `Float32`. `FullyReducedTensor` is constructed from the `x`-node
  projection `Px` (K3). `getindex` and materialisation into a preallocated output (`copyto!`)
  allocate 0 bytes for all three. `PotentialReducedTensor` defines no `size` or `axes` past its
  own rank, so `size(rt, d)` and `axes(rt, d)` answer as `AbstractArray` does, and the view of
  `VelocityReducedMatrix` answers through `SubArray` (K15). No `Base.materialize` method exists:
  `collect` and `Array` materialise them (K16). All three are exported; `ReducedTensor` is not,
  and stays `ReducedBasisMethods`' type. A grid mismatch cannot arise through `bracket_tensor`,
  which reads the grid from its bracket (K17). A tensor written by hand is assumed to be on its
  bracket's grid, as the docstrings state.

- **`scripts/verify_conservation.jl`** measures the two conservation claims of the manuscripts
  and separates them, because they are not the same claim and the obvious reading is wrong.
  The particle sums `Σ w v̇` and `Σ w v v̇` vanish at round-off on **every** boundary condition,
  since `A₁` and `A₂` are solved from exactly those constraints. What the span requirement of
  `main.tex:226` controls is the *other* claim, `∫ vᵏ f_s dv = Σ w v^k`: measured on a cubic
  basis over `[-8,8]`, `Free()` reproduces `1`, `v`, `v²` to `3e-14`, `Periodic()` gets `v`
  wrong by `1.2` and `v²` by `0.12`, and `Dirichlet()` gets even the constant wrong by `4e-2`.
  So the boundary condition decides whether `f_s` has the moments the particles do, and hence
  whether the entropy and the H-theorem refer to the distribution being evolved.

- **`polynomial_reproduction`** is re-exported from `SimpleSplines`, and
  `check_conservation_basis` turns it into an assertion. `ConservativeLenardBernstein` warns
  once when built on a basis that cannot reproduce `v` and `v²`.

- **`project_function`** projects any function onto a distribution's basis, and
  **`project_Maxwellian`** now works in any number of velocity dimensions. It previously existed
  only for `VD == 2` and called `Integrals.solve` with `HCubatureJL`, neither of which was
  imported, so it raised `UndefVarError` whenever it was reached.

- **`L` is now checked against the double sum it rearranges.** `compute_L!` expands the two
  gradient differences of `eq:discrete-landau-matrix` into four terms and folds them in pairs
  using the symmetry of `U`. The suite previously asserted only that `L` was finite and
  symmetric, both of which survive a wrong folding; it now also evaluates the defining
  `O(Q²M²)` sum directly on a deliberately tiny basis (`M = 16`, `Q = 36`) and compares. This
  is what makes the identity claim in the docstring a verified one.

- **A deposition allocation test.** `test/distributions/spline_distribution.jl` asserts that the
  cost of `projection` is independent of the particle count, so a boxed closure in `_deposit!`
  cannot come back unnoticed — nothing else in the suite would see it, since the results stay
  correct and only the run time changes. Guarded on `--check-bounds=auto`, because `Pkg.test()`'s
  default `=yes` inflates allocation counts and would make the ceiling meaningless.

- **`sodeensemble`, the training data of a `Vlasov–Poisson` model as an `EnsembleProblem`.** It
  is exported. `sodeensemble(model, members, tspan, tstep)` takes one `ParticleDistribution` per
  member and returns a `GeometricEquations.SODEEnsemble`, which
  `GeometricIntegrators.integrate(ens, Strang())` solves into one
  `GeometricSolutions.EnsembleSolution`.

  - Each member keeps its own initial condition `(q = copy(m.particles.z),)` and its own weights
    `(w = copy(m.particles.w),)`, which are what the deposit reads: an importance-sampled draw
    gives each member different ones, so weights held in the shared model would be wrong for
    every member but one. They are copied, so a later change to a member's distribution does not
    reach the problem. A member whose weights and state name different particle counts throws an
    `ArgumentError`, and so does an empty member list. A member of another element type than the
    model throws a `MethodError`.
  - The problem carries the Hamiltonian
    `h(t, q, w) = ½Σᵢ wᵢvᵢ² + ½∫(∂ₓϕ)²dx` as its invariant, so `invariants(ens)` is a `NamedTuple`
    and not `NullInvariants()`. `hashamiltonian` stays `false`, the equation being an `SODE`. The
    field term is the Gauss quadrature of `(∂ₓϕ)²` at the quadrature points of the potential's
    basis, tabulated once in the `VlasovHamiltonian` constructor, so a call allocates nothing. Every
    term of that sum is non-negative, so the sum keeps the relative accuracy of its terms.
  - The weights alone are what HDF5 can store, so the solution round-trips with
    `GeometricBase.h5save` and `h5load` against the same problem, per-member parameters
    included. Their `EnsembleSolution` methods come from the HDF5 extension of
    `GeometricSolutions` 0.6.6, so `GeometricSolutions` joins
    `test/Project.toml` as a test-only dependency. `GeometricBase` is already a dependency of the
    package; `test/Project.toml` lists it so that the test can load `h5save` and `h5load` from it.
    Nothing in `src/` calls `GeometricSolutions`.

- **`plot_distribution` and `plot_particles` are exported and drawn by a CairoMakie package
  extension.** Their methods live in `ext/VlasovMethodsMakieExt.jl`, which loads with `CairoMakie`;
  `src/plotting.jl` declares the two names and documents them. `plot_distribution(dist::GridDistribution)`
  draws `dist.values` as one heatmap on the grid's own nodes, untransposed.
  `plot_distribution(dist::ParticleDistribution{T,1,1}; xdomain, vdomain, nx = 100, nv = 100)` draws
  the weighted particle histogram: each cell holds its summed weight divided by `hx·hv`, `x` is
  wrapped into the periodic `xdomain`, and a particle whose `v` lies outside `[va, vb)`, or whose
  `x` is not finite, is dropped rather than clamped into an end cell.
  `plot_particles(dist::ParticleDistribution{T,1,1}; xdomain, vdomain)` draws the wrapped `x` against
  `v`, with the axis limits set to the two domains. Both particle methods throw `ArgumentError` on a
  distribution with no particles. `CairoMakie` is a weak dependency, so the package installs and
  loads with no Makie; without the extension loaded, each name has no methods, and calling it throws
  `MethodError`. `[compat] CairoMakie` is `"0.15"`, the only version the extension is tested
  against.

- **`ElectricField`, the electric field of a particle Vlasov–Poisson problem, and `PoissonField`,
  the field of a `PoissonSolvers.Potential`.** Both are exported. `ElectricField` is the interface
  a reduced field lifts into. `update!(f, x, w, t)` deposits the positions `x` with the weights
  `w` and re-solves the potential; `efield!(f, e, x)` writes ``E(x) = -\partial_x\phi(x)`` at each
  position into `e`. The call forms are `f(e, x, w, t)` (both steps), `f(e, x)` (`efield!` alone)
  and `f(x, w, t)` (both steps, into an `e` the call allocates). `update!` and `efield!` are not
  exported; reach them as `VlasovMethods.update!` and `VlasovMethods.efield!`. `PoissonField`
  implements the interface on a B-spline potential, and its `update!` and `efield!` are
  type-stable and allocation-free in `Float64` and `Float32`. `energy(f)`, the field energy
  ``\tfrac12\int(\partial_x\phi)^2\,\mathrm{d}x``, is allocation-free too.

- **`ReducedElectricField` and `DEIMElectricField`, the reduced electric field and its hyper-reduced
  counterpart.** Both are exported. `ReducedElectricField(field, Ψ, X)` lifts the reduced
  coordinates `X` by `Ψ`, deposits and solves, evaluates the field at every lifted position, and
  projects back by `Ψ'`. At full basis size, lifted back by `Ψ`, it equals the full-order field to
  `√eps(T)`, the only size at which the two must agree exactly. `DEIMElectricField(field, Ψₓ, Ψₑ,
  p, X)` is the DEIM path. `p` holds the DEIM row indices of the field basis `Ψₑ`, from
  `ReducedBasisMethods.get_DEIM_indices`, and the field is lifted by `Ψₓ[p, :]` and projected by
  `ΨₓᵀΨₑ Ψₑ[p, :]⁻¹`, so no dense `n × kₑ` matrix `Πₑ` is formed. With `Πₑ = I[:, p]`, the path is
  ``E = \Psi_x^T \Psi_e (\Pi_e^T \Psi_e)^{-1} E(\Pi_e^T \Psi_x X)``, with `X` the reduced
  coordinates: the field is evaluated at the `kₑ` interpolation points `p` rather than at every
  particle. It agrees with the reduced path to the DEIM error bound of
  `ReducedBasisMethods.get_DEIM_interpolation_matrix`, and it is measurably cheaper: `0.86 ms`
  against `2.12 ms` per field evaluation at 5·10⁴ particles and `k = 20` (measured by
  `scripts/electric_field_timing.jl`). Both are type-stable and allocation-free in `Float64` and
  `Float32`, and their buffers follow the element types of the bases and of the coordinates. The
  three-argument constructor accepts a `Matrix` basis. Interpolation points are chosen by the
  `max|r|` selection of `ReducedBasisMethods` 0.2, so `ΠₑᵀΨₑ` stays nonsingular when a column's
  largest entry is negative. `DEIMElectricField(field, Ψₓ, Ψₑ, Πₑ, X)` takes the matrix of
  `get_DEIM_interpolation_matrix` instead of `p`, reads the row of each unit column, and calls the
  index form; the two give bit-identical fields. The `[compat]` floor of `ReducedBasisMethods` is
  `0.2.1`, the first release that exports `get_DEIM_indices`, the function that produces `p`.

- **`reduced_integrate_vp`, the reduced particle march of the analytic quadratic model, with the
  χ-scaled timestep; `AnalyticHamiltonian`, its Hamiltonian in the forms of
  `ReducedBasisMethods.reduced_ode` and `reduced_pode`; and `lift_solution`, which lifts one of
  their solutions back to `(q, p)`.** All five of `reduced_integrate_vp`, `IntegratorParameters`,
  `ReducedIntegratorCache`, `AnalyticHamiltonian` and `lift_solution` are exported.

  `reduced_integrate_vp(P₀, Ψ, params, efield, SS, IP, IC, p)` advances the reduced state of a
  `ParticleList` — `Ψ'·vec(P₀.x)` and `Ψ'·vec(P₀.v)` — `IP.nₜ` steps of `IP.dt` by a leapfrog: half
  an advection, the reduced acceleration of a `ReducedElectricField` at the new position, and
  the second half advection. The step is `Δt = IP.dt·T(params.χ)`, and the times handed to the
  field are unscaled: the kick of step `it` is at its midpoint `(it - ½)·IP.dt`, and a save of
  step `it` is at `it·IP.dt`. `χ` scales the timestep and nothing else; it is converted to the
  march's element type before the product, so a `Float32` march with a `Float64` `χ` steps in
  `Float32`. For a time-independent field, `(dt, χ = 2)` equals `(2dt, χ = 1)` state for state.
  The initial state is saved at `zero(T)` and then every `nsave = nₜ/(nₛ - 1)`-th step, which is
  10001 states for `nₜ = 10^7` and `nₛ = 10001`. The constructor of `IntegratorParameters`
  refuses, with an `ArgumentError`, `nₜ < 1`; `nₛ < 2`; and an `nₛ - 1` that does not divide
  `nₜ`. The HDF5 constructor applies the same checks to the stored counts. The march refuses,
  before the loop, snapshots that do not hold `nₛ` saved times; a component count other than one,
  which the scalar momentum `M = dot(w, v)` needs; a basis with a zero column count or more columns
  than rows; and a `χ` that is not `> 0`, which refuses `NaN` too. `params` without a `χ`
  raises the `NamedTuple`'s own error: there is no default. The step and the save are internal,
  type-stable and allocation-free.

  `AnalyticHamiltonian(β)` is `½|p|² + ½β²|q|²` for `x = [q; p]`. Called as `∇H!(g, t, x, params)`
  it writes the gradient whose first block is `β²q` and second `p` — in the argument order
  `reduced_ode` and `reduced_pode` call, and under a name that is not RBM's `hamiltonian_gradient!`
  — and called as `h(t, x, params)` it returns the value. It accepts any element type, which the
  Lobatto IIIA–IIIB(2) run of the PODE form needs: Newton's `JacobianAutodiff` passes it
  `ForwardDiff.Dual` states. `lift_solution(b, sol)` returns the two blocks of
  `lift(b, [sol.q; sol.p])`, each of `N` rows and one column per stored time; for a cotangent lift
  those blocks are `Φ*sol.q` and `Φ*sol.p`, and a complex SVD mixes the two solution halves.

  The march and the reduced forms are tested at `T = 1000` in the new `slow` group: the leapfrog at
  `k = 20` on the analytic model against its closed form, and the PODE form — PSD cotangent, `k =
  20`, Lobatto IIIA–IIIB(2) — which holds the Hamiltonian to `3.6e-7` over the horizon. The ODE
  form — POD, `k = 5`, RK4 — runs beside it, and its separation from the bound is logged, not
  asserted: its drift is round-off, so whether and when it leaves the bound depends on the
  platform. `src/particles/time_marching.jl`, which had
  never loaded, now does: its dead `PoissonSolvers` import and its two unused `LinearAlgebra` names
  go, it is included after `particles/electric_field.jl`, and the march decides its saves on the
  step index, where the file it replaces tested a counter it advanced only on a save, which wrote
  state 0 and nothing after it for `nsave ≥ 2`. `KNOWN_ISSUES.md` K1 then covers
  `src/particles/poisson.jl` alone.

### Documentation

- **The manual gains a page on projecting between the three representations of the distribution
  function.** `docs/src/projections.md` is the map of the directions between
  `ParticleDistribution`, `SplineDistribution` and `GridDistribution`: every direction that a
  method implements runs as an executed `@example` block, so the page fails its own build when a
  call drifts from the API, and every direction that no method implements is named as absent
  rather than left out. Deposition is the worked example of `test/projections/projections.jl` —
  the same million particles, seed, sampler and tolerance — and it ends in a CairoMakie figure of
  the particle histogram, the projected density and the analytic density on one axis. The page
  names `spline → particles` as absent, and K54 records why. `docs/Project.toml` gains CairoMakie,
  PoissonSolvers, SimpleSplines and Random, and a `[sources]` entry so that the examples run
  against the checkout instead of the registered package. The build warns that the navbar carries
  no repository link, which K55 records.

### Internal

- **`KNOWN_ISSUES.md` K1 covers the two files of `src/particles/` that stay out of the module.**
  `src/particles/snapshots.jl` is included now, and its `ParameterSpace` binding comes from
  `ReducedComplexityModeling`, so the clause that named both leaves the entry; the comment above
  the new `include` in `src/VlasovMethods.jl` names what remains. `src/particles/electric_field.jl`
  and `src/electric_field.jl` are included now as well, so K1 covers `poisson.jl` and
  `time_marching.jl` alone.

- **The `ProgressMeter` compat floor rises to `1.3.3`.** The CairoMakie test dependency brings in
  TiffImages, whose 0.10 releases require `ProgressMeter` 1.3.3 or later, and the Downgrade job
  resolves every floor to its oldest admissible version: with `ProgressMeter = "1"` it reports
  `Unsatisfiable requirements detected for package TiffImages` and runs no test.

- **The `MultiIndexArrays` compat floor rises to `0.1.2`.** `ReducedBasisMethods` 0.2, which this
  release requires, depends on `MultiIndexArrays` in `0.1.2 - 0.1`, and the Downgrade job resolves
  every floor to its oldest admissible version: with `MultiIndexArrays = "0.1.1"` its forcedeps
  check reports `MultiIndexArrays resolved to 0.1.2 but lower bound is 0.1.1` and fails.

- **CI uploads coverage from the `Julia 1 - ubuntu-latest` job instead of `Julia min`, and a
  test job saves the Julia cache only when it succeeds.**

- **Test suite reorganized to mirror `src/` structure.** Test files grouped under
  `test/distributions/`, `test/projections/`, and `test/gridbased/`; dependencies moved from
  `Project.toml` to `test/Project.toml`, where `OffsetArrays` returns as a test-only dependency.
  Pass totals unchanged (257, 53, 875, 3, 12). New `test/quality/aqua.jl` runs `Aqua.test_all`
  with ambiguity checking marked broken (issue #57); its `stale_deps` check runs unmarked and
  passes, which settles K13. Unreachable `test/electric_field_tests.jl` is now
  `test/electric_field.jl` in the broken group, with a fixed seed; it loads, then throws
  `UndefVarError` at line 19 before its first `@test` (issue #56). `test/profile.jl` moved to
  `scripts/profile.jl`.

- **A runtime dispatch on the projection path now fails the suite.** New `test/quality/jet.jl`,
  in the `core` group after `quality/aqua.jl`, runs `JET.report_opt` on the vector and the
  matrix method of `projection` onto a `SplineDistribution` — the calls that
  `test/distributions/spline_distribution.jl` asserts with `@allocated` — at the `Float64`
  argument types of those tests. Neither reports, so both are
  `@test`. Where JET does not work on the running Julia, the file records one `@test_skip`. `JET`
  joins `test/Project.toml` with no `[compat]` bound, and the five `[compat]` entries there for
  dependencies of `Project.toml` are removed, since the package's own bounds apply in the test
  environment: `GeometricBrackets = "0.1.1"`, `LinearAlgebra = "1"`,
  `PoissonSolvers = "0.4, 0.5, 0.6"`, `Random = "1"` and `SimpleSplines = "0.1, 0.2, 0.3"`. JET
  loads Revise, which prints EMFILE errors into the test log; K18 records it.

- **The collision operators are deduplicated.** The four `Cache`/`CacheType` method pairs that
  re-type a collision cache for a new element type are generated by one `@collision_cache` macro
  in `src/cache.jl` instead of being written once per cache, and the dead `Cache`/`CacheType`
  methods on the models themselves go (unreachable — `CacheDict`'s parent is always the cache).
  About 570 lines of commented-out and dead code go with them. Nothing about the arithmetic
  changes: a one-step pin of each of the five collision operators on a `StableRNGs`-seeded input
  matches the pre-change result to `rtol = 64*eps(Float64)`, in
  `test/integration/collision_operators.jl`. The testset "Re-type a collision cache to Float32"
  in that file re-types the cache of each of the four models the macro covers to `Float32`, and
  asserts the new spline's element type and that the particle distribution is shared. This covers
  `KNOWN_ISSUES.md` K46, which leaves the file.
