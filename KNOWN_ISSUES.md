# Known issues

What is known to be wrong in VlasovMethods and is not fixed. An entry leaves this file when its
fix merges, and the CHANGELOG entry of that fix names its ID. IDs are never reused: a new entry
takes the number after the highest `K<n>` that this file, `CHANGELOG.md` or their history has
ever named, which `git log -p -- KNOWN_ISSUES.md CHANGELOG.md` shows, so the ID of an entry that
has left the file, or that lived only on a merged branch, stays taken.

### K4 · The appendix Landau method takes the implicit midpoint rule, which is not a discrete gradient, so no discrete H-theorem is proven for it.

- **location:** `src/methods/Landau_solver.jl`
- **evidence:** The manuscript derives two methods. The main-text method is `PairwiseLandau`:
  the gradient form with the `G` operator of `eq:G_semi-disc` at the midpoint, whose structure
  *is* the momentum and energy conservation proof, and the Gonzalez discrete gradient of
  `discrete_gradient_step!`, which *is* the discrete H-theorem proof. The appendix method is
  `Landau`: the two-step `v̇ = K⁺LJ` with the implicit midpoint rule and `∇S` at the midpoint
  (`eq:particle_ode_midpoint`), the appendix's own choice, which is not a discrete gradient. The
  H-theorem of the time discretisation therefore holds for `PairwiseLandau` only. A gap between
  the appendix method and the main text's proof, not an error in either.
- **kind:** docs
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

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

- **location:** `src/models/lenard_bernstein.jl:60`
- **evidence:** `LB_rhs!(v̇, v, params, t)` keeps the argument order of an ODE solver that the
  package does not use, and nothing in `src/` calls it: `LB_rhs_GI!` (`:73`) computes the same
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

- **location:** `src/gridbased/vlasov_poisson_1d1v.jl:29`
- **evidence:** `_green_matrix` throws `DimensionMismatch` when `length(b) != nx` (`:29`), and
  `ArgumentError` when `b.xgrid[1:nx] != dist.x` (`:33`). The second comparison also fails on a
  wrong length, so one guard with one exception type would do. The `VlasovPoisson` docstring, the
  CHANGELOG and the guard test name `DimensionMismatch`, so the merge changes a documented
  exception. Whether to merge the two guards is open.
- **kind:** defect
- **found:** #72

### K52 · The grid `VlasovPoisson` reads two internal fields of PoissonSolvers.

- **location:** `src/gridbased/vlasov_poisson_1d1v.jl:22`
- **evidence:** `_green_matrix` reads `potential.solver` (`:22`) and `b.xgrid` (`:33`). Neither
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
  step cannot be zero`, raised by `Base`'s `(::Colon)` from `src/methods/splitting.jl:47`, and the
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

### K94 · The figure drivers of the Lenard–Bernstein manuscript stop where the projection is not positive at a particle: figures 1, 2, 7 and 8 at the manuscript's parameters, and figure 11 with the seed 1234, where the right-hand side of `LenardBernstein` divides by `f_s` without checking its sign.

- **location:** `src/models/lenard_bernstein_conservative.jl:144`
- **evidence:** `compute_coefficients` of `ConservativeLenardBernstein` throws a `DomainError`
  where `f_s ≤ 0` at a particle (`:144`), and the manuscripts' parameters reach it:
  `julia --project=scripts scripts/pd_lb_fig1_2.jl` prints `nbasis =     25: the operator stops:
  the projected distribution is non-positive at v = 4.055196391323677`, and the same for
  `npart =   5623` at `v = 3.799579553401782` and `npart =  31623` at `v = 4.432320203818066`; at
  `M = 25` the projected Maxwellian is `−2.4e-5` at two of the `10⁵` velocities. The sweeps leave
  these points out of figures 1 and 2, record `stopped_at = 25` and `stopped_at = 5623`, and
  exit 0. `julia --project=scripts scripts/pd_lb_fig7_8.jl` stops in the step 1768
  to `t = 1.4144`, where the projection is `−1.58e-7` at `v = −6.42`; it writes
  `pd_lb_fig7_8.hdf5` with `stopped_at = 1.4144`, writes no figure and exits non-zero.
  `LB_rhs_GI!` evaluates `f_s'/f_s + v` at each particle (`src/models/lenard_bernstein.jl:80`)
  and does not refuse `f_s ≤ 0`. Where an iterate of the Newton solve of the implicit midpoint
  rule meets a small or negative `f_s` in the tail, the ratio sends a particle out of `[−10, 10]`,
  and the projection throws its `DomainError` naming the domain rather than the sign of `f_s`.
  `julia --project=scripts scripts/pd_lb_fig11.jl --seed 1234` stops at step 7 (`t = 0.035`) with
  `DomainError with -10.0 .. 10.0: 1 of 1000 particles left the velocity domain`, while the seeds
  1, 2 and 3 complete the 200 steps, and the same run with the package's own `RK438` integrator
  and the seed 1234 completes them.
- **kind:** defect
- **found:** 2026-11-07

### K95 · The Picard solve of `MetriplecticLenardBernstein` does not converge at the parameters of the three tests of the metriplectic Lenard–Bernstein manuscript, so each of the six runs of `scripts/lenard_bernstein_metriplectic.jl` stops before its final time.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:492`
- **evidence:** `julia --project=scripts scripts/lenard_bernstein_metriplectic.jl --case <case>
  --bracket <bracket>` stops with `the Picard solve of the metriplectic Lenard–Bernstein step did
  not converge` (`:492`), writes its record with `stopped_at`, writes no figure and exits
  non-zero: the shifted normal at `t = 0.0048` (Dougherty) and `t = 0.2208` (Jeyakumar) of 1, the
  double Maxwellian at `t = 0.1248` and `t = 0.2968` of 5, and the uniform case at `t = 0.0399`
  and `t = 0.8046` of 1. Five of the six stop after 4 to 9 of the 1000 iterations, at residuals
  from `2.4e-5` to `4.2e-4`, and the uniform Dougherty run after the 1000 iterations at a residual
  of `1.07e-10`. A smaller step does not remove them: with
  `--quick --case double-maxwellian --bracket dougherty --npart 200`, the run completes at
  `h = 8e-4` and stops at `t = 0.0036` with `--tstep 4e-4`, and with `--bracket jeyakumar` and 200
  particles it stops at `t = 0.016`. The damped fixed-point solve of the driver's rate check,
  which freezes the bracket's moments, stalls on 1 of its 1250 steps under the Dougherty form and
  on 2 under the Jeyakumar form at the shifted normal's parameters, with the largest residuals
  `3.12e-3` and `2.07e-4`. `julia --startup-file=no --project=scripts
  scripts/lenard_bernstein_metriplectic_rate_stalls.jl` probes each stalled step, step 350
  (Dougherty) and steps 209 and 689 (Jeyakumar), and finds its largest residual at a particle in
  the lower tail, `v = −1.19`, `−1.56` and `−1.18` against the mean `1.98` and the standard
  deviation `1.19`, `1.07` and `1.19` of the last iterate, with 10, 2 and 10 particles within
  `±0.5`. Where a particle in a sparse tail sits where `f_s` is small, the score of `log|f_s|`
  changes fast across it, and the fixed-point map of the implicit midpoint step is not a
  contraction.
- **kind:** defect
- **found:** 2026-11-07

### K96 · `collisional_vectorfield!` of `Landau` returns `v̇_c = ν (K_c \ LJ)`, a pivoted-QR solve for each velocity component `c` alone, with the opposite sign of the appendix's `v̇ = −K⁺LJ` for the stacked matrix `K`, and at `γ = −3` its runs drift the energy by `1.7e-8` to `2.1e-7`, where the runs of `PairwiseLandau` keep it below `2e-15`.

- **location:** `src/models/landau.jl:589-593`
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
  at `:575-577` says it is: `‖K₁ \ LJ − pinv(K₁) LJ‖ / ‖K₁ \ LJ‖ = 0.588`, and `0.431` against
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

- **location:** `src/models/lenard_bernstein_conservative.jl:229`
- **evidence:** For each of the three models, `Base.return_types` of the constructor is a
  `UnionAll` over `ET <: EquationProblem` and `IT <: GeometricIntegrator`, or a `Union`, with or
  without the method argument and on `main` as on the branch; the object that the call returns
  has a concrete type. The constructors (`src/models/lenard_bernstein.jl:95`,
  `src/models/lenard_bernstein_conservative.jl:223`,
  `src/models/rescaled_lenard_bernstein_conservative.jl:144`) only call `ODEProblem` and
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

- **location:** `src/models/landau.jl:117`
- **evidence:** with the `Float32` fixture of `test/models/landau.jl`,
  `Landau(pdist, CollisionEntropy(sdist); γ = -5.999999999)` builds a model with
  `landau.γ = -5.999999999`, and `Float32(-5.999999999) == -6.0f0`, the exponent for which the
  double integral of `𝕃` diverges in two dimensions; `γ = -6.0` throws the `ArgumentError`. The
  check also has no upper bound: a finite `γ` above `-6` passes for any size. No test and no caller
  uses a `γ` near `-6`. The fix checks `T(γ)` with the `T` of the particle state.
- **kind:** defect
- **found:** #114
