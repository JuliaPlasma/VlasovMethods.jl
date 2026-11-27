# Known issues

What is known to be wrong in VlasovMethods and is not fixed. An entry leaves this file when its
fix merges, and the CHANGELOG entry of that fix names its ID. IDs are never reused: a new entry
takes the number after the highest `K<n>` that this file, `CHANGELOG.md` or their history has
ever named, which `git log -p -- KNOWN_ISSUES.md CHANGELOG.md` shows, so the ID of an entry that
has left the file, or that lived only on a merged branch, stays taken.

### K8 · Positivity of `f_s` is not preserved, and the conservative Lenard-Bernstein operators refuse `f_s ≤ 0`.

- **location:** —
- **evidence:** The metriplectic paths run where the L² projection of a particle distribution
  changes sign, as the manuscripts do. `CollisionEntropy`, `compute_J!` of `Landau` and of
  `MetriplecticLenardBernstein` and `compute_entropy` take `log|f_s|`, and throw a `DomainError`
  only where `f_s` is zero or not finite at a quadrature node. `compute_L!` of `Landau` takes the
  mobility `|f_s|` and checks no node: a zero contributes zero mobility, and a non-finite value
  enters the matrix. `compute_coefficients` of `ConservativeLenardBernstein` and
  `compute_coefficients_rclb` of `RescaledConservativeLenardBernstein` throw a `DomainError` where
  `f_s ≤ 0` at a particle, because both operators assume `f_s > 0`: the pointwise ratio
  `f_s'/f_s` is undefined at `f_s = 0`, and the PD-LB manuscript states the Lenard-Bernstein
  H-theorem for a strictly positive `f`. So a run of either operator whose projection undershoots
  stops. On `f_s < 0` the entropy `∫ f log|f|` is not convex, so both manuscripts note that the
  uniqueness of the discrete equilibrium is not guaranteed; the entropy derivative and the
  H-theorem do not need `f_s > 0`. The real fix — a positivity-preserving projection — is the
  manuscripts' own open problem.
- **kind:** defect
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K9 · The rank hypothesis behind `K K⁺ = I` is no longer checked anywhere.

- **location:** —
- **evidence:** The step to `eq:particle_ode_final` needs `K` to have full row rank. The old code
  detected the failure with two SVDs per vector-field evaluation, printed, and proceeded
  regardless; the check is gone from the hot path and has not been given a home in a diagnostic
  script.
- **kind:** missing test
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K10 · No test covers a structure-preservation claim for Landau over a run.

- **location:** `scripts/verify_conservation.jl`
- **evidence:** `git grep -in 'landau' scripts/verify_conservation.jl` returns nothing, so the
  script exercises the conservative Lenard-Bernstein operator only. The `Landau` testsets under
  `test/` assert assembly, allocations, refusals, type parameters, a reference value and the Picard
  step. Two of them test a structure property: "𝕃 is positive semi-definite on a sign-changing
  projection" of `test/distributions/spline_distribution.jl` checks that `𝕃` is symmetric,
  positive semi-definite and annihilates the constants, and "the step converges and conserves" of
  `test/integration/landau_solver.jl` checks the momentum and the energy of one step on a
  lattice. No test checks an invariant or the entropy over a run. This entry is about the appendix
  model `Landau`: for the main-text model `PairwiseLandau`, "S_h decreases over 20 steps" of
  `test/methods/discrete_gradient.jl` checks the entropy, the momentum and the energy over a run.
- **kind:** missing test
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K18 · Revise prints EMFILE errors in the test log

- **location:** `test/quality/jet.jl`
- **evidence:** JET 0.12 loads Revise, and its file watcher runs out of file handles.
  `grep -c 'UNHANDLED TASK ERROR.*EMFILE'` on a `run-tests.jl full` log of the suite with
  `test/quality/jet.jl` counts 5 blocks, each an
  `IOError: FolderMonitor: too many open files (EMFILE)` stack trace, on Julia 1.13.1 with
  JET 0.12.2. The same count on a `run-tests.jl full` log of the suite without
  `test/quality/jet.jl` gives 0. The test totals do not change.
- **kind:** upstream
- **found:** 2026-10-02

### K19 · A particle or spline `VlasovPoisson` cannot be built with a grid `Potential` basis.

- **location:** `src/models/vlasov_poisson.jl`
- **evidence:** `PoissonSolvers.FFTWBasis` and `FiniteDifferenceBasis` are not
  `SimpleSplines.AbstractBSplineBasis`, so neither `local_width` nor `_local_buffers` has a method
  for them, and the charge deposit needs one.
  `VlasovPoisson(ParticleDistribution(1, 1, 10), Potential(FFTWBasis((0.0, 1.0), 16)))` throws
  `MethodError: no method matching _local_buffers(::FFTWBasis{Float64, …}, ::Type{Float64})`.
  The constructor allocates the deposit buffer, so the error is raised at construction.
  `Potential(PeriodicBasisSpline(…))` and `Potential(DirichletBasisSpline(…))` are the supported
  bases.
- **kind:** defect
- **found:** #60

### K24 · The Lenard-Bernstein right-hand sides keep an uncalled `LB_rhs!`.

- **location:** `src/models/lenard_bernstein.jl:80`
- **evidence:** `LB_rhs!(v̇, v, params, t)` keeps the argument order of an ODE solver that the
  package does not use, and nothing in `src/` calls it: `LB_rhs_GI!` (`:93`) computes the same
  right-hand side on the buffers of the model's cache. Outside `src/`, the tests and the reference
  generator call `LB_rhs!` as an oracle for `LB_rhs_GI!`, so it stays until they have another.
- **kind:** dead code
- **found:** #61

### K27 · The collision stencils wrap the bounded `v`-grid periodically.

- **location:** `src/gridbased/collisions.jl:66`
- **evidence:** The four `v`-stencils — `QuadraticCollisions`' call operator (`:66-67`),
  `ReducedCollisionTensor`'s `getindex` (`:119-120`) and the two `_get_MC̃_*` assemblers
  (`:174-175`, `:224-225`) — wrap with `mod1`, but the `v`-grid is bounded. At the wrap the second
  difference of `v` does not vanish, so the two ends couple. Every conservation claim of the two
  docstrings is conditional on `f` vanishing at the ends of the `v`-grid to `eps(T)`. On a grid
  that does not vanish there, the moments drift: on `nx = 4` with
  `v = range(-3, 3; length = 5)`, for the shifted Maxwellian of the invariants testset in
  `test/gridbased/collisions.jl`, the end value is `0.21` of the maximum. Over ten RK4 steps the
  relative drift of the cubic form is `3.7e-3` in the momentum and `2.9e-3` in the energy, and
  that of the quadratic form is `1.1e-3` in the energy, against `eps(T)` on the
  `range(-10, 10; length = 25)` grid of that testset. No code in the file treats the boundary of
  the `v`-grid.
- **kind:** defect
- **found:** #62. Carried from ReducedBasisMethods, and recorded when the operator was repaired.

### K28 · `CollisionTensor` repeats `PoissonTensor` of GeometricBrackets and does not tie `DT` to its values.

- **location:** `src/gridbased/collisions.jl:3`
- **evidence:** The review of #62 found equal entries and equal timings for `CollisionTensor` and
  `GeometricBrackets.PoissonTensor{DT}`, so `ReducedCollisionTensor` could take the latter and the
  type could go. Independent of that: the constructor `CollisionTensor(DT, nx, nv, f)` (`:8-10`)
  takes `DT` as a free argument and checks nothing about the values `f` returns, so the element
  type of the array and the type of its entries can differ. `size(ct)` (`:13`) builds a vector
  with `ones(Int, 3)` and splats it into a tuple, so it infers `Tuple{Vararg{Int}}`. Whether to
  replace the type is open.
- **kind:** defect
- **found:** #62

### K29 · The two `_get_MC̃_*` assemblers carry arguments they do not use.

- **location:** `src/gridbased/collisions.jl:157`
- **evidence:** `h₁` is never read in either assembler (`:157-192`, `:209-238`). `ci` is read only
  for `size(ci)` (`:158`, `:210`) and for the length check of `v`. The quadratic assembler
  does not read `∫vdv` either. `grep -rn '_get_MC̃_' src test scripts` finds the calls in
  `src/gridbased/collisions.jl` and `test/gridbased/collisions.jl` only. The nine-argument
  signature follows the advisor decision recorded in #62. Whether to shorten it is open.
- **kind:** dead code
- **found:** #62

### K32 · The `getindex` of `ReducedCollisionTensor` rebuilds three indices in the innermost loop.

- **location:** `src/gridbased/collisions.jl:107`
- **evidence:** The loop at `:117-135` builds `M`, `N`, `O` and their linear indices for every
  `(m1, m2, o2, n2)`, although `M`, `O`, `m` and `o` do not depend on `n2`, and the `v`-stencil
  of `m2` is the one that `QuadraticCollisions` computes (`:66-67`). The loop can be shorter. No
  benchmark compares the loop with a shorter form.
- **kind:** not verified
- **found:** #62

### K34 · Two testsets check words of a docstring.

- **location:** `test/gridbased/collisions.jl:294`
- **evidence:** `the docstrings name the invariants` (`:294-301`) and
  `the docstrings state the v-end requirement` (`:303-308`) test `occursin` on the text of
  `@doc`, so a rewording that keeps the facts fails them and a docstring that names the words
  but is wrong passes them. The plan of the part names both checks. Whether to keep them is open.
- **kind:** missing test
- **found:** #62

### K36 · A metriplectic Picard step whose particles leave the velocity domain stops with the model's `DomainError`.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:488`
- **evidence:** With `ν = 1e6` at `ti = 4`, all 64 particles leave the `-2.0 .. 2.0` velocity
  support and `_projection!` throws `DomainError: … 64 of 64 particles left the velocity domain …`
  from inside the solve. The catch at `:488` rethrows it, so the caller sees that `DomainError`
  and not the residual-and-count `ErrorException` the author decided for a solve that does not
  converge or meets a `NaN`. The error is loud, so the solve never returns in silence.
- **kind:** found late
- **found:** 2026-10-04

### K41 · `Picard_iterate_over_particles!` calls two `SimpleSolvers` names that are not public.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:489`
- **evidence:** `SimpleSolvers.status` (`:489`) and `SimpleSolvers.isconverged` (`:492`) are neither
  exported nor declared `public`; `SimpleSolvers` 0.14.1 keeps them unexported at
  `src/SimpleSolvers.jl:138`. `ExplicitImports.check_all_qualified_accesses_are_public` flags
  exactly these two. The fix is a `public` declaration in `SimpleSolvers`.
- **kind:** upstream
- **found:** 2026-10-04

### K42 · `Picard_iterate_over_particles` takes `dv`, `m` and `β` and ignores them.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:509`
- **evidence:** The wrapper does not read `dv`, `m` or `β` after its signature (`:509-512`); it
  passes the other arguments to `Picard_iterate_over_particles!`, which takes none of the three.
  They stay so that the call sites do not change, and the CHANGELOG says so. Removing them is an
  API change.
- **kind:** dead code
- **found:** 2026-10-04

### K48 · The collision test and its reference generator define the same fixture twice.

- **location:** `test/integration/collision_operators.jl:24`
- **evidence:** `setup_1d1v` (`:24-33`), `lb_params` (`:35`) and `setup_1d2v` (`:37-48`) of the
  test repeat `setup_1d1v` (`:32-41`), `lb_params` (`:43`) and `setup_1d2v` (`:76-87`) of
  `test/helpers/generate_collision_reference.jl`, and both files define `SEED`, `N` and `N1`. A
  change to one copy that is not made in the other makes the test compare against references
  from a different input. A shared file is not a simple fix: the generator runs on the tree
  whose arithmetic is the reference, which can be a tree without the shared file.
- **kind:** missing test
- **found:** #69

### K49 · The `which(f!, …)` assertion of the metriplectic collision test cannot fail.

- **location:** `test/integration/collision_operators.jl:95`
- **evidence:** `f!` has one method, `src/models/lenard_bernstein_metriplectic.jl:394`, and its
  last argument is `mlb::MetriplecticLenardBernstein`. So `which(VlasovMethods.f!, …)` at
  `:95-97` either returns that method, and the `@test` at `:98` passes, or throws a
  `MethodError`. The value pin at `:91` already runs `f!` through the Picard step. Whether to
  keep the assertion is open.
- **kind:** missing test
- **found:** #69

### K50 · `VlasovPoisson` does not warn on an empty particle list.

- **location:** `src/models/vlasov_poisson.jl:11`
- **evidence:** An empty particle list is a vacuum run: the deposit gives `ρ ≡ 0` and the solve
  `ϕ ≡ 0`, and no error occurs. The author decided that `VlasovPoisson` warns in that case, and
  the inner constructor at `:11-15` has no `@warn`. Three questions are open: the text of the message;
  whether `sodeensemble` also warns for an empty member; and whether `VlasovWithExternalFields`
  warns too. It does not: `src/models/vlasov_with_external_fields.jl:67` and the "vacuum run"
  bullet of `CHANGELOG.md` say that its constructor does not warn.
- **kind:** defect
- **found:** #71

### K51 · The grid `VlasovPoisson` checks the grid size and the grid nodes with two guards and two exception types.

- **location:** `src/gridbased/vlasov_poisson_1d1v.jl:82`
- **evidence:** `_green_matrix` throws `DimensionMismatch` when `length(b) != nx` (`:82`), and
  `ArgumentError` when `b.xgrid[1:nx] != dist.x` (`:86`). The second comparison also fails on a
  wrong length, so one guard with one exception type would do. The `VlasovPoisson` docstring, the
  CHANGELOG and the guard test name `DimensionMismatch`, so the merge changes a documented
  exception. Whether to merge the two guards is open.
- **kind:** defect
- **found:** #72

### K52 · The grid `VlasovPoisson` reads two internal fields of PoissonSolvers.

- **location:** `src/gridbased/vlasov_poisson_1d1v.jl:75`
- **evidence:** `_green_matrix` reads `potential.solver` (`:75`) and `b.xgrid` (`:86`). Neither
  field has a public accessor in PoissonSolvers, so a change of the field names in PoissonSolvers
  breaks the constructor. The fix is a public accessor upstream.
- **kind:** upstream
- **found:** #72

### K54 · Projecting a spline distribution back onto particles is not implemented.

- **location:** `src/projections/distribution.jl:200`
- **evidence:** The method raises
  `error("projecting a spline distribution back onto particles is not implemented")`. Two meanings
  are available for the direction and neither is chosen: (a) draw a new particle sample from
  `f_s`, the reverse of the projection onto the spline, which is the direction a reduced model
  would lift a reduced state back onto particles with; (b) keep the particle positions and set the
  weights `w_α ∝ f_s(v_α)`, which is a quadrature of `f_s` at the positions that happen to be
  there rather than a sample of it. Each needs a rule for where `f_s < 0`, which K8 records as
  open. The docs page on the projections names the direction as absent and writes the two
  meanings out.
- **kind:** defect
- **found:** #73

### K56 · `Parameter` of `ReducedComplexityModeling` rejects distinct bounds for a one-sample range.

- **location:** `ReducedComplexityModeling` v0.5.0, `src/parameters/parameter.jl:55`
- **evidence:** `Parameter(name, minimum::DT, maximum::DT, n::Int)` builds
  `LinRange(minimum, maximum, n)`, which throws for `n == 1` when the endpoints differ:
  `julia --startup-file=no --project -e 'using ReducedComplexityModeling; Parameter(:β, 5.95, 6.05, 1)'`
  prints `ArgumentError: range(5.95, stop=6.05, length=1): endpoints differ` on the RCM 0.5.0 this
  package resolves, while `Parameter(:β, 5.95, 6.05, 2)` gives
  `Parameter{Float64, Vector{Float64}}(:β, 5.95, 6.05, [5.95, 6.05])`. A one-sample sweep is
  supported through `Parameter(:β, [6.0])` or equal bounds, `Parameter(:β, 6.0, 6.0, 1)`:
  each gives a `ParameterSpace` of length 1 and `size(Snapshots(Float64, 2, 100, 16, 11, ps).X)`
  is `(2, 100, 11, 1)`.
- **kind:** upstream
- **found:** #75

### K61 · A non-finite reduced coordinate or particle position raises `InexactError` instead of returning a value.

- **location:** `src/particles/electric_field.jl`
- **evidence:** A `NaN` or an `Inf` in the reduced coordinates `X` reaches `evaluate_all!`,
  where the cell index is not an integer. Calling a `ReducedElectricField` with
  `X[1] = T(NaN)` or `X[1] = T(Inf)` raises `InexactError`: the recorded output is
  `Float32 Inf X => InexactError: Int64(NaN32)`, `Float32 NaN X => InexactError: Int64(NaN32)`,
  `Float64 Inf X => InexactError: Int64(NaN)` and `Float64 NaN X => InexactError: Int64(NaN)`.
  A `NaN` or an `Inf` position reaches the same error through a `PoissonField`. Neither field
  guards non-finite input, so a march that produces it stops at the step that deposits it.
- **kind:** found late
- **found:** #78, which adds no guard and no test for non-finite coordinates.

### K65 · `best_and_worst` names one parameter as both the best and the worst when a reduced growth rate is `NaN`.

- **location:** `src/particles/diagnostics.jl:220`
- **evidence:** `best_and_worst([1.0, 2, 3], [1.0, NaN, 3.5])` returns `(2, 2)`. A `NaN` in
  `abs.(γr .- γh)` leaves `argmin` and `argmax` at the same index. `NaN` and `Inf` propagate by
  design, with no guard, so the pair of indices is the caller's to read.
- **kind:** found late
- **found:** 2026-10-17

### K66 · `saved_indices` returns a wrapped three-quarter index for a very large `nₛ`.

- **location:** `src/particles/diagnostics.jl:116`
- **evidence:** `saved_indices(4611686018427387903)` returns
  `(4611686018427387903, -1152921504606846976, 2305843009213693952)`: `3 * (nₛ - 1)` overflows
  `Int` above about `2^62`. The docstring states this closed form, so no wider arithmetic is added.
- **kind:** found late
- **found:** 2026-10-17

### K69 · The legend of `plot_series` overflows its axis for many labelled series.

- **location:** `ext/VlasovMethodsMakieExt.jl:160`, `plot_series`
- **evidence:** at the default figure size of `600 × 450`, the legend of 20 labelled series is
  481 px high in an axis of 375 px, on `05431f10`; 10 series give 251 px.
  `scripts/bump_on_tail_testing.jl` draws its 30-series figures of items 26 and 27 without
  labels and adds a legend of one entry per model. Its figures of items 24 and 25 label one series
  per sample, so at the default of 10 samples no figure of the scripts shows the overflow.
- **kind:** found late
- **found:** 2026-10-23

### K76 · A numeric flag of `scripts/options.jl` is taken as given, so a bad value fails inside the run.

- **location:** `scripts/options.jl:48`
- **evidence:** a numeric default is read with `parse(T, value)` and never range-checked, so
  `julia --project=scripts scripts/vlasov_poisson.jl --tstep 0` exits 1 on `ArgumentError: range
  step cannot be zero`, raised by `Base`'s `(::Colon)` from `src/methods/splitting.jl:75`, and the
  message names neither the flag nor the value. The drivers that parse a `--tstep` are
  `lenard_bernstein.jl`, `lenard_bernstein_metriplectic.jl`, `profile.jl` and `vlasov_poisson.jl`.
- **kind:** defect
- **found:** 2026-10-26

### K77 · The Downgrade job cannot add CairoMakie to the manifest that its floor resolve leaves.

- **location:** `.github/workflows/CI.yml`, the `Downgrade - ubuntu-latest` job
- **evidence:** the `forcedeps` step passes every direct floor and freezes a manifest in which
  `ReducedComplexityModeling` 0.5 brings `KernelAbstractions` 0.9.45, LLVM 10.4.2 and `CEnum`
  1.0.0. The test step then adds the test-only `CairoMakie` with `ALLOW_RERESOLVE: false` and stops
  on `Unsatisfiable requirements detected for package FreeType`: `CEnum` 1.0.0 admits only FreeType
  2.1, and `CairoMakie` 0.15 needs FreeType 3.0–4.1.1 (PR #96, run 41146492438, job
  124372699434). The resolve of the test environment with CairoMakie in it has a solution: every
  other CI job resolves and passes. The test-only dependencies are not in the floor resolve, so the
  job cannot see this before the test step.
- **kind:** upstream
- **found:** #96

### K79 · `GeometricIntegrator(model, tspan, tstep, method)` for a 1D1V grid model does not infer a concrete type.

- **location:** `GeometricEquations` v0.21.6, `src/problems/equation_problem.jl:82`
- **evidence:** the constructor builds its `ODEProblem` with GeometricEquations
  (`src/gridbased/vlasov_poisson_1d1v.jl:129-131` on `8937e91`). `EquationProblem` sets its first
  type parameter with `superType = eval(typeof(equ).name.name)`, which inference does not see
  through. `Base.return_types` of `GeometricIntegrator(model, tspan, 0.05, RK4())` on a 64 × 128
  linear Landau model gives a `Union` of two `UnionAll` types over
  `EquationProblem{_A, …} where _A<:GeometricEquation`. An `ODEProblem` with no parameters infers
  `EquationProblem{_A, …}` too. The `@NamedTuple{arakawa::DataType, …}` in the `ODE` type
  parameter is also the run-time type, so it is not the cause. The object built at run time has
  concrete field types, so a caller that steps it behind a function barrier runs at full speed. A
  step loop in the same function as the call dispatches at run time: on 64 × 128 one `RK4()` step
  takes 2.72–2.76e-4 s, and of the 10.7 kB that a step allocates, about 2 kB is boxing and
  dynamic dispatch. A function barrier around the step loop changes the time per step by a factor
  of 0.99–1.06, within the noise of the machine. Measured in the experiment
  VlasovReductionBenchmark, which records it as its K4.
- **kind:** upstream
- **found:** 2026-10-31

### K81 · The Lenard-Bernstein right-hand sides of a `Float32` state differ from those of a `Float64` state by 1.6e-5 to 3.6e-5 on 2000 particles, and no cause is located.

- **location:** `src/projections/distribution.jl:28`
- **evidence:** the deposit `_deposit!` adds into the load vector with a plain sum in the state's
  element type: `rhs[j] += w * buf[t]` (`:28`) in one velocity dimension and
  `rhs[idx...] += val` (`:55`) in several, once per particle. The particle sums of
  `compute_coefficients` and `compute_coefficients_rclb` use compensated summation, and the
  right-hand sides did not change with it. Probe on `32dbd62`, and again on `6dc3a9f` with the
  same figures, Julia 1.13, 17 cubic `Free`
  basis functions on `(-2, 2)`, velocities uniform on `(-2, 2)` from `StableRNG(20261102)`, equal
  weights `1/N`, a `Float64` model: `max|v̇32 − v̇64| / max|v̇64|` of the `Float32` state against
  the `Float64` state is 2.7e-6 (`LenardBernstein`), 3.2e-6 (`ConservativeLenardBernstein`) and
  3.1e-6 (`RescaledConservativeLenardBernstein`) at N = 256, and 1.6e-5, 2.3e-5 and 3.6e-5 at
  N = 2000. The `Float64` state against a `BigFloat` state is 2.2e-15 to 5.6e-14 at those sizes.
  `LenardBernstein` has no coefficient sums and shows the same gap. Replacing the plain sum at
  `:28` with Kahan summation changes the `Float32` gaps by less than 1e-9 and leaves their
  digits as above, so the deposit of the one-dimensional case is not the cause. The gap stays
  when the `Float64` evaluation takes the velocities rounded to `Float32` (1.3e-6 to 1.4e-6 at
  N = 256, 1.4e-5 to 3.2e-5 at N = 2000), so the rounding of the input is not the cause either.
  The multi-dimensional deposit at `:55` is not probed. Which `Float32` step sets the gap is
  not measured.
- **kind:** not verified
- **found:** #101

### K82 · The package module relies on 79 implicit imports, and imports names that their modules do not declare public.

- **location:** `src/VlasovMethods.jl`
- **evidence:** with ExplicitImports 1.15.0 and `GeometricMachineLearning` loaded,
  `check_no_implicit_imports(VlasovMethods)` raises an `ImplicitImportsException`, and
  `explicit_imports(VlasovMethods)` lists 79 names for the package module: those of the `using`
  lines of `HDF5`, `LinearAlgebra`, `Parameters`, `ParticleMethods`, `PoissonSolvers`,
  `ProgressMeter`, `Random`, `SimpleSplines`, `Sobol`, `SpecialFunctions` and `StaticArrays`.
  `check_all_qualified_accesses_are_public(VlasovMethods)` names `Base.Callable`
  (`src/examples/bump.jl`), `SimpleSplines.ContractWorkspace` (`src/models/landau.jl`) and
  `HDF5.H5DataStore`, and `check_all_explicit_imports_are_public(VlasovMethods)` names
  `TanhActivation` of `ext/VlasovMethodsHNNExt.jl`. `test/quality/explicit_imports.jl` checks
  the stale imports and the owners of the package, the implicit, stale, owner and public checks
  of `MultiModalAutoEncoderExt` and the implicit and stale checks of `VlasovMethodsMakieExt`, and
  leaves the implicit and public checks of the package
  out, with a comment that names this entry.
- **kind:** defect
- **found:** 2026-11-03

### K83 · The odd latent coordinates of the tied grid and low-rank frontends agree with the spline frontend at a lower order than the even ones.

- **location:** `ext/MultiModalAutoEncoderExt.jl`
- **evidence:** `BilinearStiefelLayer` and `LowRankStiefelLayer` share one velocity basis between
  the even and the odd coordinates, and the spline factor applies `R⁺` to the even half and `R⁻`
  to the odd half, with `R⁺ ≠ R⁻` near `v = 0`. `tie_continuum_frontends!` therefore ties the
  velocity bases through the even modes. Probe at `p = 3`, `V = 8`, with the modes of
  `tie_errors` in `test/rom/multimodal_vlasov_frontends.jl` and `f = (1 + 0.1 cos x) v M(v)`:
  `‖p_spline − p_grid‖` is 2.8e-2, 7.8e-3, 2.0e-3 and 5.1e-4 at `N = 17, 33, 65, 129`, an order
  of 2, where the even coordinates of `(1 + 0.1 cos x) M(v)` reach a fitted order of 4.3. A tie
  of both parities needs a velocity basis per parity in the layers.
- **kind:** found late
- **found:** 2026-11-03

### K87 · The entropy paths refuse `f_s = 0`, although `f log|f|` has a limit there.

- **location:** `src/entropies/collision_entropy.jl`
- **evidence:** `CollisionEntropy` and `compute_entropy` of `MetriplecticLenardBernstein`
  (`src/models/lenard_bernstein_metriplectic.jl`) call `_check_log_defined` at every quadrature
  node, which throws a `DomainError` where `f_s` is zero. The integrand `f log|f|` tends to `0` as
  `f → 0`, so both could take the value `0` at such a node by continuity and accept it.
  `compute_J!` of `Landau` and of `MetriplecticLenardBernstein` projects `1 + log|f_s|` itself,
  which has no finite value at `f_s = 0`, so its refusal stays. The testsets "a zero projection
  throws" of `test/distributions/spline_distribution.jl` and
  `test/models/lenard_bernstein_buffers.jl` assert the current behaviour of the entropy paths.
- **kind:** defect
- **found:** #111

### K88 · The right-hand side of `MetriplecticLenardBernstein` throws under the Dougherty bracket where the computed `ε − u²` is not positive, as for all particles at one velocity, where its exact value is a finite zero.

- **location:** `src/models/lenard_bernstein_metriplectic.jl`
- **evidence:** `metriplectic_rhs!` calls `_check_velocity_spread` before `bracket_rhs!` for
  either bracket, and the check throws a `DomainError` where the computed `ε − u²` is at most
  zero or the variance is at most `eps(T) ε`. Under `DoughertyBracket`, the terms of
  `bracket_rhs!` carry the factors `ε − u²`, `ε − u v` and `v − u`, which for particles at one
  velocity `v = u` are all zero in exact arithmetic, so the exact right-hand side is a finite
  zero and the throw is stricter than the form needs. Under `JeyakumarBracket` the factor
  `n/(ε − u²)` is undefined there, so the throw is needed. The testsets "a degenerate velocity
  set throws" and "a computed ε − u² ≤ 0 throws" of
  `test/models/lenard_bernstein_metriplectic.jl` assert the throw for both brackets.
- **kind:** defect
- **found:** #112

### K91 · The anisotropy rate of the 1600-particle bi-Maxwellian does not converge under mesh refinement, and the discrete gradient step does not converge at 12 × 12, 24 × 24 and 32 × 32 cells.

- **location:** `test/examples/reference_cases.jl`, testset "the rate at 16 × 16 cells, and the
  8 × 8 control"
- **evidence:** the rate of `fitted_rate` (30 steps, `abstol = 1e-7`, the default `max_iters`) on
  the sample `rate_sample(1600)` with `reference_model(v, w; ncell, L = 3.0)` at `Δt = 1e-5`,
  Julia 1.13.3:

  | cells | fitted rate |
  |--:|:--|
  | 8 × 8 | 10.246 |
  | 12 × 12 | the step throws "did not converge in 100 evaluations", last residual 9.6e-6 |
  | 16 × 16 | 8.1073 |
  | 20 × 20 | 7.4254 |
  | 24 × 24 | the step throws, last residual 3.3e-5 |
  | 32 × 32 | the step throws, last residual 2.0e-6 |

  The closed rate is `8`. The 16 × 16 mesh is within the testset's 5 % clamp, and 20 × 20 is not,
  so the agreement at 16 × 16 is not shown to be a converged one. The control of the testset is
  the 8 × 8 row, and it shows only that the rate depends on the mesh. On `[-4, 4]²` at 8 × 8
  cells the rate is 1.859 at `Δt = 1e-4` and 1.9574 at `Δt = 1e-5`.
- **kind:** defect
- **found:** #119

### K94 · The figure drivers of the Lenard–Bernstein manuscript stop where the projection is not positive at a particle: figures 1, 2, 7 and 8 at the manuscript's parameters, figures 3 to 6 at some seeds, and figure 11 with the seed 1234, where the right-hand side of `LenardBernstein` divides by `f_s` without checking its sign.

- **location:** `src/models/lenard_bernstein_conservative.jl:172`
- **evidence:** `compute_coefficients` of `ConservativeLenardBernstein` throws a `DomainError`
  where `f_s ≤ 0` at a particle (`:172`), and the manuscripts' parameters reach it:
  one run of `julia --project=scripts scripts/pd_lb_fig1_2.jl` at its default `--seed 1`, with
  Julia 1.13.4 on 2026-11-20, prints `nbasis =     25: the operator stops: the projected
  distribution is non-positive at v = 3.9517739700081314`, and the same for `npart =   3162` at
  `v = -3.280449770438173` and `npart =  17783` at `v = -3.815087083306259`. The sweeps leave
  these points out of figures 1 and 2, record `stopped_at = 25` and `stopped_at = 3162`, and
  exit 0. Which points stop depends on the sample: over the seeds 1 to 20, `M = 25` stops at
  every seed, and 1 to 8 of the 13 values of `N` stop.
  `julia --project=scripts scripts/pd_lb_fig3_6.jl --seed 1` stops in the second step of its run
  of 1000 particles (`t = 0.0016`), where the projection is `−2.66e-4` at `v = 4.93`. Over the
  seeds 1 to 20, the run of 1000 particles to `t = 1` stops at 7 (1, 4, 7, 8, 9, 12, 14), the run
  of 200 particles at 4 (4, 9, 11, 19); the driver's default seed 2 completes both.
  `julia --project=scripts scripts/pd_lb_fig7_8.jl` stops in the step 1768 to `t = 1.4144`, where
  the projection is `−1.58e-7` at `v = −6.42`; it writes `pd_lb_fig7_8.hdf5` with
  `stopped_at = 1.4144`, writes no figure and exits non-zero.
  `LB_rhs_GI!` evaluates `f_s'/f_s + v` at each particle (`src/models/lenard_bernstein.jl:100`)
  and does not refuse `f_s ≤ 0`. Where an iterate of the Newton solve of the implicit midpoint
  rule meets a small or negative `f_s` in the tail, the ratio sends a particle out of `[−10, 10]`,
  and the projection throws its `DomainError` naming the domain rather than the sign of `f_s`.
  `julia --project=scripts scripts/pd_lb_fig11.jl --seed 1234` stops at step 7 (`t = 0.035`) with
  `DomainError with -10.0 .. 10.0: 1 of 1000 particles left the velocity domain`, while the seeds
  1, 2 and 3 complete the 200 steps, and the same run with the package's own `RK438` integrator
  and the seed 1234 completes them.
- **kind:** defect
- **found:** 2026-11-07

### K95 · The Picard solve of `MetriplecticLenardBernstein` does not converge at the parameters of the three tests of the metriplectic Lenard–Bernstein manuscript, so five of the six runs of `scripts/lenard_bernstein_metriplectic.jl` at the default seed 1 stop before their final time.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:492`
- **evidence:** at the default `--seed 1`, with Julia 1.13.4 on 2026-11-20,
  `julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case <case>
  --bracket <bracket>` stops with `the Picard solve of the metriplectic Lenard–Bernstein step did
  not converge` (`:492`), writes its record with `stopped_at`, writes no figure and exits
  non-zero: the shifted normal with the Jeyakumar form at `t = 0.0016` of 1, after 4 iterations at
  the residual `1.6e-3`; the double Maxwellian at `t = 0.1248` and `t = 0.2968` of 5, after 4
  iterations at `4.2e-4` and `1.6e-4`; and the uniform case at `t = 0.0399` and `t = 0.8046` of 1,
  after the 1000 iterations at `1.07e-10` and after 5 at `2.4e-5`. The shifted normal with the
  Dougherty form completes. Which runs stop depends on the sample: over the seeds 1 to 20 the
  shifted normal stops at 11 seeds with the Dougherty form and at 13 with the Jeyakumar form,
  between `t = 0.0008` and `t = 0.41`. A smaller step does not remove them: with
  `--quick --case double-maxwellian --bracket dougherty --npart 200`, the run completes at
  `h = 8e-4` and stops at `t = 0.0036` with `--tstep 4e-4`, and with `--bracket jeyakumar` and 200
  particles it stops at `t = 0.016`. The damped fixed-point solve of the driver's rate check,
  which freezes the bracket's moments, stalls on 2 of its 1250 steps under the Dougherty form and
  on 4 under the Jeyakumar form at the shifted normal's parameters and seed 1, with the largest
  residuals `2.36e-3` and `1.45e-3`. `julia --startup-file=no --project=scripts
  scripts/lenard_bernstein_metriplectic_rate_stalls.jl` probes each stalled step, steps 94 and 98
  (Dougherty) and steps 174, 196, 412 and 1068 (Jeyakumar), and finds its largest residual at a
  particle in the lower tail, `v = −1.17` and `−1.30`, and `−1.23`, `−1.30`, `−1.64` and `−2.12`,
  against the mean `2.00` and standard deviations from `1.06` to `1.25` of the last iterate, with
  2 to 7 particles within `±0.5`. Where a particle in a sparse tail sits where `f_s` is small, the
  score of `log|f_s|` changes fast across it, and the fixed-point map of the implicit midpoint
  step is not a contraction.
- **kind:** defect
- **found:** 2026-11-07

### K96 · `collisional_vectorfield!` of `Landau` returns `v̇_c = ν (K_c \ LJ)`, a pivoted-QR solve for each velocity component `c` alone, with the opposite sign of the appendix's `v̇ = −K⁺LJ` for the stacked matrix `K`, and at `γ = −3` its runs drift the energy by `1.7e-8` to `2.1e-7`, where the runs of `PairwiseLandau` keep it below `2e-15`.

- **location:** `src/models/landau.jl:613-617`
- **evidence:** `eq:particle_ode_final` and `eq:particle_ode_midpoint` of the manuscript's
  appendix give `v̇ = −K⁺ L J`; the code solves `v̇[c, :] = K_c \ LJ` for each component `c` with
  `qr!(…, ColumnNorm())` and scales by `ν`, with no minus sign. On the weighted lattice of
  `N = 60²` particles on `[−3.4, 3.4]²` with `16 × 16` quadratic cells on `[−3.5, 3.5]²`
  (`M = 324`), at `γ = 0`, `julia --project=scripts scripts/landau_compare.jl` measures the
  instantaneous `(d gap/dt) / gap` of `gap = T₁₁ − T₂₂` at `+8.0055` for `Landau` and at `−8.0095`
  for `PairwiseLandau`, against the closed form `−8` of `scripts/reference_forms.jl`: under
  `Landau` the anisotropy grows at the rate at which it should decay. On the same particles, after
  one evaluation, `rank(K₁) = 302` and `rank(K₂) = 306` of `M = 324`, while the stacked `[K₁ K₂]`
  has the rank `323 = M − 1` that `dN ≥ M − 1` gives; each `K_c K_c⁺` projects onto the range of
  `K_c` alone, so `‖K₁ v̇₁ + K₂ v̇₂ − LJ‖ / ‖LJ‖ = 0.73` where `K K⁺ = I` of the appendix would make
  it zero. The pivoted-QR solve of the rank-deficient `K_c` is not `pinv(K_c) LJ`, as the comment
  at `:599-601` says it is: `‖K₁ \ LJ − pinv(K₁) LJ‖ / ‖K₁ \ LJ‖ = 0.588`, and `0.431` against
  `‖pinv(K₁) LJ‖`; for `K₂` the two are `3.03` and `0.889`. The sign is the code's either way: each
  `K_c \ LJ` lies nearer `+pinv(K_c) LJ` than `−pinv(K_c) LJ`. `scripts/landau_cases.jl` therefore
  refuses `--rate-check` with `--method spline`.

  The energy `Σ w |v|² / 2` drifts under `Landau` at `γ = −3`: the runs
  `julia --project=scripts scripts/landau_cases.jl --method spline --case <case>`, 1000 particles,
  100 steps of `1e-8` on 8 cells per axis, each solved to the absolute residual `1e-10`, end at
  the relative energy drift `2.1e-7` for `bimaxwellian`, `5.1e-8` for `double-maxwellian` and
  `1.7e-8` for `bkw`, where the same runs with `--method pairwise` end at `1.5e-15`, `1.1e-16`
  and `3.4e-16`.
- **kind:** defect
- **found:** 2026-11-07

### K100 · Inference does not give a concrete type for `GeometricIntegrator(model, tspan, tstep[, method])` of `LenardBernstein`, `ConservativeLenardBernstein` and `RescaledConservativeLenardBernstein`, because the `ODEProblem` constructor of GeometricEquations and the Newton solver of SimpleSolvers infer free type parameters.

- **location:** `src/models/lenard_bernstein_conservative.jl:257`
- **evidence:** For each of the three models, `Base.return_types` of the constructor is a
  `UnionAll` over `ET <: EquationProblem` and `IT <: GeometricIntegrator`, or a `Union`, with or
  without the method argument and on `main` as on the branch; the object that the call returns
  has a concrete type. The constructors (`src/models/lenard_bernstein.jl:115`,
  `src/models/lenard_bernstein_conservative.jl:251`,
  `src/models/rescaled_lenard_bernstein_conservative.jl:167`) only call `ODEProblem` and
  `Integrators.GeometricIntegrator`, and both calls infer free parameters on a toy problem, with
  GeometricEquations 0.21.6, GeometricIntegrators 0.18.7 and SimpleSolvers 0.14.3:
  ```julia
  using GeometricEquations: ODEProblem
  using GeometricIntegrators: ImplicitMidpoint, Integrators
  f!(v, t, q, p) = (v .= -q; v)
  build(q) = ODEProblem(f!, (0.0, 1.0), 0.1, q; parameters = (ν = 1.0,))
  equ = build(rand(5))
  pt = only(Base.return_types(build, (Vector{Float64},)))
  it = only(Base.return_types(Integrators.GeometricIntegrator, (typeof(equ), ImplicitMidpoint)))
  println(isconcretetype(typeof(equ)), " ", pt.var, " ", it.var)
  ```
  prints `true _A<:GeometricEquations.GeometricEquation ST<:(SimpleSolvers.NewtonSolver{…, JT,
  …, LiSeT, …} where {JT<:(SimpleSolvers.JacobianAutodiff{Float64, …,
  ForwardDiff.JacobianConfig{Nothing, Float64, _A, …}, …} where {_A, …}), LiSeT<:…})`: the
  problem's equation parameter is free with or without the keyword `parameters`, and the solver's
  is free through the chunk size of `ForwardDiff.JacobianConfig(nothing, y, x)`
  (`SimpleSolvers/src/base/jacobian.jl:153`), which depends on the length of `x` at run time.
- **kind:** upstream
- **found:** #129

### K101 · On Julia 1.14.0-alpha2, Pkg refuses the multi-range `[compat]` entries of `test/Project.toml` under `force_latest_compatible_version`, so the jobs "Julia pre" and "Julia nightly" fail on every Dependabot pull request.

- **location:** `test/Project.toml:52`, `JET = "0.10.1, 0.11, 0.12, 0.13"`, and `:57`,
  `StableRNGs = "0.1, 1"`
- **evidence:** julia-runtest sets `force_latest_compatible_version = true` on a Dependabot pull
  request, and on Julia 1.14.0-alpha2 `Pkg.test(; force_latest_compatible_version = true)` stops
  with `inconsistency between compat values and string representation` before any test runs (PR
  #131, run 42444532599, job 128620570142,
  https://github.com/JuliaPlasma/VlasovMethods.jl/actions/runs/42444532599/job/128620570142).
  Neither job is a required check; `main` and other pull requests pass both, because the option is
  off there. One `[compat]` entry with more than one range is enough. The CI steps of a Dependabot
  pull request, instantiate and then
  `Pkg.test(; force_latest_compatible_version = true, allow_reresolve = true)`, under Julia
  1.14.0-alpha2 fail with that error on `main` at `97b1eeb` and on the head `0e11494` of #131;
  on `main` they pass with both entries at one range and fail with only `StableRNGs = "0.1, 1"`.
  The cause is in Pkg: `apply_force_latest_compatible_version!` narrows the value of a compat entry
  in place and leaves its string (`Operations.jl:4910` of the Pkg in 1.14.0-alpha2), and when a
  later step writes the project, `destructure` (`project.jl:294-295`) refuses an entry whose value
  differs from `semver_spec` of its string. It is reported as JuliaLang/Pkg.jl#4902, which is
  open. One more condition of the trigger is not known: a minimal package with
  `StableRNGs = "0.1, 1"` in its test project passes. The two entries span several ranges on
  purpose, so narrowing them is not a fix. The fix is a Pkg in a Julia 1.14 release whose
  `apply_force_latest_compatible_version!` keeps the string of an entry consistent with its value;
  the entry leaves this file when the pre job passes on a Dependabot pull request.
- **kind:** upstream
- **found:** #131

### K102 · `Landau` checks the kernel exponent `γ > -6` on the `γ` it is given, but `kernel` evaluates `T(γ)` in the element type `T` of the state, so a `Float32` model accepts a `γ` that rounds to the refused `-6`.

- **location:** `src/models/landau.jl:141`
- **evidence:** with the `Float32` fixture of `test/models/landau.jl`,
  `Landau(pdist, CollisionEntropy(sdist); γ = -5.999999999)` builds a model with
  `landau.γ = -5.999999999`, and `Float32(-5.999999999) == -6.0f0`, the exponent for which the
  double integral of `𝕃` diverges in two dimensions; `γ = -6.0` throws the `ArgumentError`. The
  check also has no upper bound: a finite `γ` above `-6` passes for any size. No test and no caller
  uses a `γ` near `-6`. The fix checks `T(γ)` with the `T` of the particle state.
- **kind:** defect
- **found:** #114

### K103 · The field `caches` of `CacheDict` has the abstract value type `Cache`.

- **location:** `src/cache.jl:8`
- **evidence:** `caches::Dict{DataType, Cache}` holds the caches of one model, keyed by element
  type, on `1e60db3`. The constructor (`:10-14`) stores the model's own cache under its
  `eltype`, and the models build it as `CacheDict(LBCache(dist, ent.dist))`
  (`src/models/lenard_bernstein.jl:58`, and the same in the other collision models and
  `src/gridbased/vlasov_poisson_1d1v.jl:160`). `getindex` (`:19-25`) builds a missing cache with
  `Cache(T, parent(c))` and asserts the result `::CacheType(T, parent(c))`, so a read such as
  `params.model.cache[ST]` (`src/models/lenard_bernstein.jl:94`) infers the concrete cache type
  from that assertion and not from the field. Each read still takes a `Dict` lookup of an abstract
  value and a type check. `fieldtype(VlasovMethods.CacheDict, :caches)` returns
  `Dict{DataType, VlasovMethods.Cache}`, and `isabstracttype(VlasovMethods.Cache)` returns `true`.
  A check of concrete field types, such as `concrete-fields.jl`, flags such a cache by design: it
  reports `VlasovMethods.CacheDict.caches::Dict{DataType, VlasovMethods.Cache}`. The `CacheDict`
  of GeometricIntegratorsBase has the same pattern. The fix is a design decision: for
  example, a structure with one concrete field per element type that the model uses, in place of
  the `Dict`, decided together with GeometricIntegratorsBase.
- **kind:** defect
- **found:** 2026-11-09

### K106 · `scripts/gate.jl` requires the quick runs of the shifted normal to reach their final time, which they do at the default seeds but not at every seed, so a change that only changes the sample can fail the gate.

- **location:** `scripts/gate.jl:507`
- **evidence:** the gate runs `pd_lb_fig3_6.jl --quick` (`:507`) and
  `lenard_bernstein_metriplectic.jl --quick --case shifted-normal` with both brackets (`:521`), and
  checks that each exits 0. The quick runs take 200 particles to `t = 0.05`, and where the
  projection is not positive at a particle (K94) or a Picard solve does not converge (K95), a run
  stops and exits non-zero. Over the seeds 1 to 20, a run of 200 particles of `pd_lb_fig3_6.jl`
  stops before `t = 0.05` at 4 seeds (4, 9, 11, 19), not at its default seed 2, and the quick
  metriplectic runs reach `t = 0.05` with both brackets at 9 seeds (1, 4, 6, 7, 12, 13, 14, 17,
  19), the default seed 1 among them. A gate that reads a stop with a consistent record as a
  pass, as it reads the sweeps of `pd_lb_fig1_2.jl`, would not depend on the sample.
- **kind:** defect
- **found:** 2026-11-20

### K107 · The quick rate check of the metriplectic shifted normal meets its tolerance of 5 % at the default seed but not at every seed.

- **location:** `scripts/gate.jl:317`
- **evidence:** the gate checks that the fitted rate of each quick shifted normal run, 200
  particles to `t = 0.05`, lies within `MPLB_RATE_TOLERANCE = 0.05` of its closed form (`:870`),
  and the driver exits non-zero at its default `--tolerance 0.05`. At the default seed 1 the fits
  differ by `2.5e-2` (Dougherty) and `3.6e-2` (Jeyakumar). Over the seeds 1 to 20 they differ by
  `4.7e-3` to `1.1e-1` and `3.0e-3` to `1.5e-1`, and 6 seeds miss 5 % with one bracket or both
  (2, 6, 9, 14, 16, 19); the controls stay `0.78` and `0.49` away or more. At the manuscript's
  1000 particles to `t = 1`, every seed from 1 to 20 meets 5 %, with differences below `1e-2`.
- **kind:** defect
- **found:** 2026-11-20

### K108 · `GridDistribution(values, xdomain, vdomain)` accepts a reversed x or v domain, whose negative step gives negative moment weights and makes `LowRankDistribution(dist, r)` throw a `DomainError`.

- **location:** `src/distributions/grid_distribution.jl:54`
- **evidence:** the constructor builds `range(first(domain), last(domain))` without a check of the
  order, so `step(dist.v) < 0` for `vdomain = (10, -10)` and `step(dist.x) < 0` for
  `xdomain = (4π, 0)`. On the `24 × 40` grid of the smooth fixture of
  `test/lowrank/distribution.jl`, `velocity_moments(dist).density[1]` is `-1.69` with
  `vdomain = (10.0, -10.0)`, against `1.69` with `(-10.0, 10.0)`. `LowRankDistribution(dist, 4)`
  then throws a `DomainError` (`sqrt was called with a negative real argument`), from
  `sqrt(h)` in `src/lowrank/grid_space.jl:107` for the reversed v domain, and from
  `sqrt(hx * hv)` in `src/lowrank/distribution.jl:114` for the reversed x domain. Whether a
  reversed domain is refused or reordered is not decided; the check belongs in the constructor of
  `GridDistribution`, so that every consumer of its grid gets it.
- **kind:** defect
- **found:** 2026-11-21

### K110 · The edge, rank-guard and run-driver testsets of the BUG tests run on 16, 24 or 32 nodes of `v ∈ [−10, 10]`, where f at `±10` does not stay below `1e-14` of its largest value.

- **location:** `test/lowrank/bug.jl:242`
- **evidence:** the fixtures of "a rank outside the bounds, $T" (`:168`), "the edges of the step,
  $T" (`:242`, `:311`, `:333`, `:362`, `:371`, `:384`) and "the run driver, $T" (`:442`). After
  5 steps of `τ = 0.05` from the rank-3 factorisation `LowRankDistribution(dist, 3)` of
  `(1 + 0.3 cos(x/2)) M(v)` on `16 × 24` nodes of `[0, 4π) × [−10, 10]` with the model's Green's
  matrix (the fixture of `:242` and `:442`), f at `±10` is `8.1e-14` (`BUG()`), `7.9e-14`
  (`AugmentedBUG()`) and `4.2e-12` (`MidpointBUG()`) of its largest value in Float64, and
  `8.3e-9`, `4.4e-6` and `1.2e-6` in Float32; at `t = 0` it is `2.1e-22`. The plan's velocity-domain rule asks for round-off, and
  its spec names `1e-14`; the fixtures of the numerical bounds run on 64 velocity nodes for this
  reason. Whether the boundary values change an assertion of these testsets is not checked.
- **kind:** not verified
- **found:** 2026-11-23

### K112 · The element-type assertions on the factors of a `LowRankDistribution` after a step cannot fail, so a promotion inside `lowrank_step!` is not seen there.

- **location:** `test/lowrank/projector_splitting.jl:115`
- **evidence:** `@test eltype(Y.X) === eltype(Y.S) === eltype(Y.V) === T` at lines 115 and 161,
  and the same assertion on `Ya` against `Float32` at line 428. The exactness testset at line 226 asserts no element type, and an assertion on `Y` there
  could not fail either. `LowRankDistribution` is a mutable struct whose factor fields have the
  concrete types of its parameters, so the assignment `Y.X, Y.S, Y.V, Y.nfixed = …` in the step
  converts a `Float64` matrix to the field's type. A `Float32` state fails on a promoted `τ` only
  because the fields of `LowRankFields1D1V` take `Matrix{T}` alone and raise a `MethodError`;
  with a generic field, such as the test's `PathField`, a promotion inside the step would pass
  these assertions. An assertion on a value the step computes, such as the `K` that a logging
  field receives, would see it.
- **kind:** missing test
- **found:** 2026-11-23

### K113 · Seven testsets of the bounded velocity closure assert no element type of their results.

- **location:** `test/gridbased/vlasov_poisson_1d1v.jl:570`
- **evidence:** "the end columns are not read" (`:570`), "energy, reflection and mass flux"
  (`:590`), "the edges of the bounded right-hand side" (`:615`) and "the initial state of a bounded
  run" (`:675`) of `test/gridbased/vlasov_poisson_1d1v.jl`, "the bounded DLR field: full rank and
  energy" (`:307`) and "the bounded workspace" (`:568`) of `test/lowrank/fields.jl`, and "run A
  over 10⁴ steps" (`:148`) of `test/verification/bounded_closure.jl` hold no `eltype` assertion.
  Each function they call has its element type asserted in both precisions in another testset of
  the same file: `eltype(ḟ) === T` for the bounded right-hand side (`:545`, `:564`), the fields
  and `_potential!` in "the bounded fields equal the dense bounded field", and
  `eltype(Fs[end]) === T` of run A (`test/verification/bounded_closure.jl:121`). The fix adds
  `@test eltype(...) === T` to each of the seven.
- **kind:** missing test
- **found:** 2026-11-23

### K114 · The aliasing error of the bounded grid right-hand side names the arguments of GeometricBrackets' stencil, not those of the right-hand side.

- **location:** `src/gridbased/vlasov_poisson_1d1v.jl:251`
- **evidence:** `_grid_1d1v_bounded_buffered!` throws
  `ArgumentError("Pc shares memory with û or c, which the stencil reads while it writes Pc")` when
  `ḟ` shares memory with `f`, the message that `poisson_apply!` gives the periodic right-hand side.
  `_grid_1d1v_bounded_rhs!` has the arguments `ḟ` and `f`, and no `Pc`, `û` or `c`. The bounded
  path reads its copy `g` of the unknowns and not `f`, so an aliased `ḟ` would not corrupt the
  stencil; the check exists so that both right-hand sides throw one error with one message. A
  message that names `ḟ` and `f` for both right-hand sides changes the periodic message too.
- **kind:** docs
- **found:** 2026-11-23

### K115 · The bounded low-rank fields add four range-restricted loops, `_cut_gram!` and `_cut_apply!` with and without a weight, that repeat `_shifted_gram!` and `_shifted_apply!`.

- **location:** `src/lowrank/fields.jl:151`
- **evidence:** `_cut_gram!(M, Y, s, js[, u])` (`:151`, `:164`) and
  `_cut_apply!(out, Y, s, js[, u], P)` (`:177`, `:189`) are each the loop of the periodic helper restricted to the rows `js` of
  `_unknowns`, about 50 lines, and each weighted method repeats its unweighted one with one more
  factor. In `_bounded_kfield!` (`:433`), two of the three `_shifted_apply!` calls of each offset
  add the zero matrix `Z`. Two simpler forms keep every clause: copies of `V` and `L` with zero
  end rows, `O(nv r)` bytes inside the asserted allocation bound, passed to the weighted
  `_shifted_gram!` with the weight `χ[j + dv]` and followed by zeroing rows `1` and `nv` of `L̇`;
  or one method per name with a weight argument, which removes about 24 lines at one multiply per
  term in the unweighted calls.
- **kind:** defect
- **found:** 2026-11-23

### K124 · The CI job "Julia 1 - ubuntu-latest - default - verification" takes about 70 minutes on a normal runner, so a runner about 1.6 times slower or more takes it past its timeout of 120 minutes.

- **location:** `test/runtests.jl:63`
- **evidence:** `split = true` gives the `verification` group a job of its own, but the group
  alone holds about 104 minutes of test time on two workers. In CI run 45340940249 (westus, a
  head without "Projector Splitting Rates") the job took 56 minutes: 14 minutes of setup and
  precompilation with an empty cache ("No cache found", the first run of the new job; a restored
  cache shortens it), then 42 minutes of tests, with "Low-Rank BUG Landau Rates" at 1725 s,
  "Low-Rank BUG Two-Stream Rates" at 1173 s, "Low-Rank BUG Two-Stream Rates, Float32" at 1179 s,
  "Bounded Velocity Closure Runs" at 318 s, "Reduced March at T = 1000" at 285 s and "Reduced
  Equations at T = 100" at 176 s. "Projector Splitting Rates" took 1366 s in CI run 45263801182
  (westus3), which puts the job at about 66 to 79 minutes. In CI run 45283090271 (eastus2) the
  test files ran about 2.4 times slower ("Projector Splitting Rates" 3466 s, "Landau Solver" 1259 s
  against 547 s); at that speed the tests of the group alone take about 124 minutes.
- **kind:** defect
- **found:** PR #147, from its CI

### K126 · Explicit RK4 on the bounded velocity closure is unstable at τ = 0.05 on the grids of the Landau runs, because the cut of h at the end nodes raises the spectral radius of the right-hand side by a factor that grows with nᵥ (5.2, 10.7 and 21.6 at nᵥ = 64, 128 and 256); the stable τ depends on nₓ, hᵥ and vmax.

- **location:** `src/gridbased/vlasov_poisson_1d1v.jl`
- **evidence:** linear Landau damping (k = 0.5, α = 0.01) with `closure = BoundedVelocityClosure()`
  and GeometricIntegrators' `RK4()` at τ = 0.05 on 64 × 128 nodes of [0, 4π) × [−6, 6]: E₁ is not
  finite from t = 0.55. The spectral radius ρ of the linearised right-hand side over all x-modes
  is 303.4 (periodic closure: 28.3) on 64 × 128 of [−6, 6], 1048.6 (periodic: 48.5) on 64 × 256
  of [−10, 10], and 44.0 (periodic: 8.45) on 16 × 64 of [−7.5, 7.5], so RK4 needs
  τ ≤ 2.83/ρ = 0.0093, 0.0027 and 0.064 there. The periodic ρ grows with vmax too: on 64 × 255 of
  [−12, 12], RK4 at τ = 0.05 gives a non-finite E₁ from t = 24.7 and runs at τ = 0.025. The tests
  of the bounded closure step with the implicit midpoint, which is unconditionally stable, except
  one that takes three RK4 steps of τ = 0.1 on 16 × 64 nodes of [−7.5, 7.5]
  (`test/gridbased/vlasov_poisson_1d1v.jl:682`). That τ is past the limit 0.064: a longer run
  from the test's own state (α = 0.5, `vlasov_poisson_1d1v.jl:101`) gives a non-finite state
  from t = 2.0, step 20, in Float64. The test does not see this: its assertions after the run
  (`:691-692`) check that the end columns are zero and that the interior has changed, not that
  the state is finite, and a non-finite interior also differs from the initial one. A low-rank
  step with RK4 substeps on bounded fields is expected to meet the same limit; no run measures
  it. `scripts/landau_turn_level.jl` prints ρ and the limit for each grid
  (group `linear`) and runs both unstable cases in Float64, the second from the test's state
  (group `closure`).
- **kind:** found late
- **found:** 2026-11-24

### K128 · The `f⁺` solve of `ParitySplitting` cannot meet its residual test on a state whose even part is zero or small against its odd part, so `run!` throws on the first step.

- **location:** `src/gridbased/parity.jl:224`
- **evidence:** the solve stops at `‖F‖₂ ≤ tol ‖F0 + f⁺‖₂`, a bound that scales with the even part
  of the state alone, while the residual `y - f⁺ - τ P(f⁻, h)` also holds the odd part `f⁻`.
  In Float32 on 16 × 32 nodes of [0, 4π) × [−7.5, 7.5], with `τ = 0.1` and the default `tol` and
  `max_iterations`: from `(1 + 0.2 cos(x/2)) e^{-v²/2} v`, whose even part is zero, so the bound is
  0, `run!` throws `MAXITERS, with ‖F‖ = 2.204519e-9 after 10 iterations`; from
  `(1 + 0.2 cos(x/2)) e^{-v²/2} (v + 1e-6)` it throws with `‖F‖ = 4.6566284e-10`. The bound is the
  one that the spec of the method decides; a bound with a term in `‖f⁻‖` or in `τ ‖P(f⁻, h)‖`
  would reach these states. No test runs an odd-dominated state.
- **kind:** found late
- **found:** 2026-11-25; the clause "The parity Strang step", whose step 2 decides the bound

### K129 · The `f⁺` solve of `ParitySplitting` stalls above its residual bound on a state scaled by 10² or more, and in Float32 on 64 × 128 nodes at τ = 0.2; the cause is not located.

- **location:** `src/gridbased/parity.jl:221`
- **evidence:** 40 steps of `τ = 0.05` from `s (1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π)` on 16 × 32 nodes of
  [0, 4π) × [−7.5, 7.5] with the default `tol` and `max_iterations` keep the energy to `3.4e-7`
  (Float32) and `1.0e-15` (Float64) at `s = 10`, and throw `MAXITERS` at `s = 10²`, with
  `‖F‖ = 1.49e-4` (Float32) and `5.82e-13` (Float64) after 10 iterations; at `s = 10⁴` in Float64,
  `max_iterations = 100` still ends at `‖F‖ = 6.07e-4`. From `(1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π)`
  on 64 × 128 nodes of the same domain, a Float32 run at `τ = 0.2` throws at step 10 with
  `‖F‖ = 1.73e-4`, and `max_iterations = 20` or `50` (`‖F‖ = 1.75e-4` and `1.53e-4`) or
  `tol = 16 eps(Float32)` (`1.72e-4`) do not change that. The tolerance, the forcing and the
  `GMRES` defaults are the ones that the spec of the method decides; whether the round-off of the
  residual or the convergence of `GMRES` sets this floor is not measured. The tests run 16 × 32
  and 16 × 33 at the scales 1 and 10⁻⁴.
- **kind:** not verified
- **found:** 2026-11-25; the clause "The parity Strang step", whose step 2 decides the solve

### K130 · The `GalerkinMidpoint` constructor accepts a time span that no time grid can hold, and `run!` then raises an `InexactError` from `TimeSeries` instead of an `ArgumentError` of a guard.

- **location:** `src/splinebased/vlasov_poisson_1d1v.jl:265`
- **evidence:** the constructor checks `tstep`, not `tspan`. With the Landau state on 16 × 100
  splines, `GalerkinMidpoint(model, (0.0, Inf), 0.05)`, `(0.0, NaN)` and `(NaN, 1.0)` construct,
  and `run!` raises `InexactError: Int64(Inf)` or `Int64(NaN)` from
  `GeometricIntegrators.TimeSeries`, in Float64 and Float32; a Float32 span `(0, 1e39)` gives
  `Int64(Inf)`, and `tstep = nextfloat(0)` gives `Int64(Inf)` in Float32 and `Int64(7.1e44)` in
  Float64. A reversed span raises the `ArgumentError` of `TimeSeries` at `run!`. No file is
  created in these cases. The spec of the method decides the checks of `tstep`, `tol` and
  `max_iterations`, and none of `tspan`.
- **kind:** found late
- **found:** 2026-11-26

### K131 · `GalerkinMidpoint` accepts `tol = Inf`, and `run!` then returns with the state unmoved and no error.

- **location:** `src/splinebased/vlasov_poisson_1d1v.jl:276`
- **evidence:** the guard refuses a negative and a `NaN` tolerance, as the spec of the method
  decides, so `tol = Inf` passes; the target `tol ‖A‖₂` is `Inf`, the start of every solve meets
  it, and one step of `τ = 0.05` from the Landau state returns `A` unchanged without an error. A
  Float32 `tol = 1e39` converts to `Inf` in the same way.
- **kind:** found late
- **found:** 2026-11-26

### K132 · `evaluate` on a Float32 `TensorSplineDistribution` with Float64 arguments returns a Float64.

- **location:** `src/distributions/tensor_spline_distribution.jl:92`
- **evidence:** on 16 × 100 Float32 splines with random coefficients,
  `typeof(evaluate(dist, 1.3, -7.1))` is `Float64`, and `Float32` for `1.3f0, -7.1f0`; a
  `BigFloat` coordinate, `evaluate(dist, big(1.0), 1)`, gives a `BigFloat` in both precisions.
  The spec of the distribution decides no result type for arguments of another precision, and no
  test passes arguments of a wider precision than the distribution's.
- **kind:** found late
- **found:** 2026-11-26

### K133 · The element-type assertions of the run, step and stopping testsets of the spline Galerkin model cannot fail, and the step testset asserts none.

- **location:** `test/splinebased/vlasov_poisson_1d1v.jl:438`
- **evidence:** `run!(::GalerkinMidpoint, h5file)` creates the datasets `z` and `t` with the
  element type `T` of the distribution (`src/splinebased/vlasov_poisson_1d1v.jl:361-362`), and
  writes the last state into the typed coefficient field of the distribution, so
  `eltype(z) === T`, `eltype(times) === T` and `eltype(dist.coefficients) === T` hold whatever
  type a step computes. A step that promoted would raise no error there: it would be rounded on
  the write. The step testset asserts no element type; the field, residual and Jacobian-vector
  testsets of the same file assert the element type of the computed matrices. The critic's probe
  of `_midpoint_step` on a Float32 model gives `Float32`, so no step promotes on the branch.
- **kind:** missing test
- **found:** 2026-11-26

### K134 · The macOS CI jobs of the test group `galerkin` take about 120 minutes, the timeout of the job, so the job "Julia 1 - macOS-latest - default - galerkin" is red on the weekly run, a tag and a run by hand.

- **location:** `test/runtests.jl:73`
- **evidence:** in the full-matrix CI run 45789228548 (`workflow_dispatch`), "Julia 1 -
  macOS-latest - default - galerkin" was cancelled at its timeout of 120 minutes, and "Julia min -
  macOS-latest - default - galerkin" passed in 119 min 51 s; the Windows jobs of the group took
  82–87 minutes and the Ubuntu jobs 60–68 minutes. The group holds the four files
  `test/verification/spline_galerkin_*.jl`: 10⁴ implicit midpoint steps in Float64 and in Float32,
  and the Landau and two-stream rate runs on 64 × 128 splines in both precisions. `galerkin` is in
  `once` of `test/ci-groups.toml`, so a pull request and a push run it on Ubuntu × Julia 1 only,
  and the red shows only where every group runs on every entry. The user accepts this red: the
  group is not split, the runs are not shortened, and the timeout does not change.
- **kind:** defect
- **found:** 2026-11-26, the full-matrix run of the spline Galerkin model's branch

### K135 · The default background of `ParitySplitting` overflows in its sum over x when the even part of finite initial values sums past `floatmax(T)` in a column, so the constructor throws an `ArgumentError` that calls the finite values not finite.

- **location:** `src/gridbased/parity.jl:161`
- **evidence:** `_parity_background(values, nothing)` halves each value before it adds the
  mirrored value, but adds the `nx` rows before it divides by `nx`. On 16 × 32 and 16 × 33 nodes
  of [0, 4π) × [−7.5, 7.5], from `(1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π)` with rows 1 and 2 of the
  columns 3 and `nv − 2` set to `floatmax(T)`, every value and the x-mean `sum(values ./ 16;
  dims = 1)` are finite, and `ParitySplitting(model, (0, 1), 0.1)` throws `ArgumentError: the
  initial values and the background of the parity splitting must be finite`, in Float32 and
  Float64. A division before the sum, `vec(sum(values ./ (2nx) .+ values[:, end:-1:1] ./ (2nx);
  dims = 1))`, gives a finite background, and the method constructs with it; on these values the
  Hamiltonian `h` and the bracket term `P(f, h)` that `run!` evaluates are not finite.
- **kind:** found late
- **found:** 2026-11-26

### K136 · No test of `ParitySplitting` passes a background with negative values, so a mutant of the evenness bound that drops the absolute value survives.

- **location:** `src/gridbased/parity.jl:147`
- **evidence:** the constructor tests the evenness of `f0` against `10 eps(T) maximum(abs, f0)`.
  The mutant `maximum(abs, f0)` → `maximum(f0)` passes all 260 tests of
  `test/gridbased/parity.jl`. On 16 × 32 nodes of [0, 4π) × [−7.5, 7.5], from
  `(1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π)`, the background `-e^{-v²/2}/√(2π)` constructs and is stored
  as given, in Float32 and Float64; the mutant's bound on it is `-2.9e-19` (Float32) and
  `-5.4e-28` (Float64), so the mutant refuses this even background. A test that constructs with
  `background = -M` catches the mutant.
- **kind:** missing test
- **found:** PR #160

### K137 · A step of `ProjectorSplitting` with a `LowRankFields1D1V{Float32}` field on a `Float64` state throws a `MethodError` from inside a substep, not an `ArgumentError` that names the element types.

- **location:** `src/lowrank/projector_splitting_midpoint.jl:349`
- **evidence:** a rank-2 `Float64` state on 8 × 16 nodes of `[0, 4π) × [−6, 6]` and the field
  `LowRankFields1D1V` of the same grid in `Float32`, one step of `τ = 0.05` on Julia 1.13.4.
  With `ImplicitMidpointSubsteps()`, `:lie` and `:strang` both throw
  `MethodError: no method matching _potential!(::Vector{Float32}, ::Matrix{Float64}, ::Matrix{Float64}, ::LowRankFields1D1V{Float32, Matrix{Float32}, PeriodicVelocityClosure})`.
  The RK4 substeps of the same `ProjectorSplitting(:lie)` and `ProjectorSplitting(:strang)`
  behave the same way, with
  `MethodError: no method matching _kfield!(::Matrix{Float64}, ::Matrix{Float64}, ::Matrix{Float64}, ::LowRankFields1D1V{Float32, Matrix{Float32}, PeriodicVelocityClosure})`
  from `src/lowrank/projector_splitting.jl:8`. In all four cases `X`, `S` and `V` stay equal to
  their values before the step, because each step writes the state only after its last substep.
  `_check_midpoint` checks the type of the field but not its element type against that of the
  state; a check `w isa LowRankFields1D1V{eltype(Y)}` there, and the same before the RK4
  substeps, gives an `ArgumentError` before any substep.
- **kind:** defect
- **found:** 2026-11-26

### K138 · `evaluate` on a `TensorSplineDistribution` raises an `InexactError` from SimpleSplines at a non-finite coordinate, and returns 0 at a velocity outside the domain of its velocity basis.

- **location:** `src/distributions/tensor_spline_distribution.jl:92`
- **evidence:** on 16 periodic cubic x-splines of [0, 4π) and 100 clamped cubic v-splines of
  [−25, 25], with the projection of `1 + v/3 + v²/5`, `evaluate(dist, NaN, 0)` and
  `evaluate(dist, Inf, 0)` throw `InexactError`, in Float64 and Float32, and
  `evaluate(dist, 1, 30)` returns `0.0` where the projected function is 191. The method passes the
  point to SimpleSplines' `evaluate` without a check, and the docstring says only that it returns
  `f(x, v)`; the spec of the distribution decides no behaviour outside the velocity domain and at
  a non-finite point, and no test reaches either case.
- **kind:** found late
- **found:** PR #161

### K139 · No test runs `GalerkinMidpoint` on a nonuniform x-mesh; the nonuniform testset checks the energy rate of the field by a central difference only.

- **location:** `test/splinebased/vlasov_poisson_1d1v.jl:255`
- **evidence:** the testset "nonuniform x-mesh" checks the potential, the energy and the rate
  `H(A + εȦ) − H(A − εȦ)` on 16 cells of [0, 4π) of widths between 0.52 and 1.48 times the mean,
  and calls no `run!`. A probe of 40 implicit midpoint steps of `τ = 0.05` from the Landau state
  on that mesh, with the default `tol`, keeps the energy to `6.5e-16` (Float64) and `1.4e-7`
  (Float32), so the step holds there; a step that breaks only on a nonuniform mesh passes the
  suite.
- **kind:** missing test
- **found:** PR #161

### K140 · The constructor of `LowRankParityDistribution` repeats the finiteness check, the evenness check and the symmetrisation of the background that the constructor of `ParitySplitting` makes.

- **location:** `src/lowrank/parity.jl:73`
- **evidence:** lines 73–79 of `src/lowrank/parity.jl` are lines 144–150 of
  `src/gridbased/parity.jl` with the two messages changed. One helper beside
  `_parity_background` in `src/gridbased/parity.jl` would hold both copies; the low-rank parity
  part does not change that file.
- **kind:** defect
- **found:** 2026-11-27; the files of the low-rank parity part, which leave
  `src/gridbased/parity.jl` unchanged

### K141 · The default background of `LowRankParityDistribution` overflows in its sum over x as that of `ParitySplitting` does, so the constructor throws an `ArgumentError` that calls finite values not finite.

- **location:** `src/gridbased/parity.jl:161`
- **evidence:** on 16 × 32 nodes of [0, 4π) × [−7.5, 7.5], from
  `(1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π) (1 + 0.3 v)` with rows 1 and 2 of the columns 3 and
  `nv − 2` set to `floatmax(T)`, every value is finite, and
  `LowRankParityDistribution(dist, (4, 3))` throws `ArgumentError: the values and the background
  of the low-rank parity distribution must be finite`, in Float32 and Float64. The constructor
  takes the background from `_parity_background(values, nothing)`, whose sum overflows.
- **kind:** found late
- **found:** 2026-11-27; K135, whose fix closes this one too

### K142 · A state of large amplitude makes the low-rank parity splitting unstable at τ = 0.05, and its cores turn non-finite within a few steps.

- **location:** `src/lowrank/parity.jl:422`
- **evidence:** the state `10⁴ (1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π) (1 + 0.3 v)` on 16 × 32 nodes
  of [0, 4π) × [−7.5, 7.5], at rank (6, 6), with `LowRankParitySplitting((6, 6))` and τ = 0.05:
  the status `all(isfinite, Y.plus.S) && all(isfinite, Y.minus.S)` is false after step 3 in
  Float32 and after step 5 in Float64. The dense transcription of the algorithm in
  `test/helpers/dense_parity.jl` grows the norm of `f⁻` to 1.6e12, 1.8e22, 9.0e54 and 1.8e206
  over the first four steps in Float64, and its SVD fails at step 5. The plasma frequency grows
  with the square root of the amplitude, to about 100 here, so `ω_p τ ≈ 5`; the step has no
  step-size control.
- **kind:** found late
- **found:** 2026-11-27; the clause "The decided edges of the step", whose non-finite core the
  status reports

### K143 · `_least_rank(σ, tol)` throws `ArgumentError: collection must be non-empty` from `first` for a tolerance that is `NaN` or negative, a message that names neither the tolerance nor the truncation.

- **location:** `src/lowrank/bug.jl:119`
- **evidence:** `VlasovMethods._least_rank([3.0, 1.0], NaN)` and
  `VlasovMethods._least_rank([3.0, 1.0], -1.0)` throw `ArgumentError: collection must be
  non-empty`, in Float32 and Float64: no rank `r` meets `norm(σ[r + 1:end]) ≤ tol`, so the
  generator is empty. The truncation of the low-rank parity splitting returns `nothing` before
  this call when the norm of its matrix is not finite, so a step does not reach it.
- **kind:** defect
- **found:** 2026-11-27; the fix of the overflow of the truncation's tolerance in the low-rank
  parity part

### K144 · A step of a rank-adaptive `LowRankParitySplitting` that ends with a false status can return a part at another rank than it had before the step.

- **location:** `src/lowrank/parity.jl:271`
- **evidence:** the state `(1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π) (1 + 0.3 v)` on 16 × 32 nodes of
  [0, 4π) × [−7.5, 7.5] at rank (6, 6), with `Y.minus.S[1, 1] = 1.34e154` (Float64) or `1.0f19`
  (Float32), and one step of τ = 0.05 with
  `LowRankParitySplitting((6, 6); tol = (1e-8, 1e-8), minrank = (2, 1))`: the status
  `all(isfinite, Y.plus.S) && all(isfinite, Y.minus.S)` is false, and the ranks go from (6, 6) to
  (3, 1) in Float64 and to (6, 2) in Float32. The first f⁻ substep truncates a finite core to a
  lower rank, and a later substep, whose core is not finite, writes its `NaN` core at the rank
  that its part has when that substep starts. The step docstring states this rule per substep;
  the clause on a non-finite state decides the ranks only for the default method, with a `NaN`
  or `Inf` entry.
- **kind:** found late
- **found:** 2026-11-27; the critic of the fix of the overflow of the truncation's tolerance in
  the low-rank parity part
