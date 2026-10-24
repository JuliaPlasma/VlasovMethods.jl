# Known issues

What is known to be wrong in VlasovMethods and is not fixed. An entry leaves this file when its
fix merges, and the CHANGELOG entry of that fix names its ID. IDs are never reused; a new entry
takes the next `K<n>`.

### K1 · `src/particles/poisson.jl` is present but not included.

- **location:** `src/particles/poisson.jl`
- **evidence:** The file names an undefined binding, so including it breaks the load:
  `poisson.jl` wants `PoissonSolverPBSplines`, which neither the module nor its dependencies
  define.
- **kind:** dead code
- **found:** 2026-09-17. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of this is a regression; it is either a numerical-methods decision or work the
  migration deliberately did not take on.

### K4 · The implemented Landau scheme is not the one the main text derives.

- **location:** —
- **evidence:** The manuscript builds the gradient form with the `G` operator — whose structure
  *is* the momentum and energy conservation proof — and a Gonzalez discrete gradient, which *is*
  the discrete H-theorem proof. What runs is the **appendix** two-step `v̇ = K⁺LJ` with plain
  implicit midpoint and `∇S(midpoint)`, which is not a discrete gradient. Neither structural proof
  transfers to the code as written. `G` is never formed, and no code implements the gradient
  form. A gap between paper and code, not an error in either.
- **kind:** docs
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K5 · The Landau Picard solver does not iterate to convergence.

- **location:** `Landau_solver.jl`
- **evidence:** `Landau_solver.jl` runs exactly five iterations and prints the residual without
  testing it. `tol`, `ftol`, `β`, `m` and `chunksize` are accepted and unused. Every conservation
  property in the appendix is a property of the *exactly* solved implicit system, so momentum and
  energy drift at the size of that printed residual.
- **kind:** defect
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K6 · The accepted Landau solution is one Picard step behind its stored derivative.

- **location:** —
- **evidence:** The loop's last action recomputes `v̇` at the midpoint from the newest guess
  without updating the guess, so the stored state and derivative do not correspond and the next
  step's Hermite extrapolation is fed an inconsistent pair.
- **kind:** defect
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K7 · The Landau kernel's coincident-point value is a regularisation, not a limit.

- **location:** —
- **evidence:** `kernel` returns zero at `|u| = 0` so that a product quadrature sharing nodes does
  not produce `Inf`. The error does not vanish under refinement, and because both factors use the
  same Gauß–Legendre nodes it fires on every *diagonal cell pair* rather than on a set of measure
  zero. The kernel is integrable in two dimensions; what is needed is a singularity-aware rule or
  offset grids.
- **kind:** defect
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K8 · Positivity of `f_s` is detected, not solved.

- **location:** —
- **evidence:** Every `log f_s` and `1/f_s` now throws rather than continuing, which is strictly
  better than the old `0.5·log(f_s²)` returning `log|f_s|`. But a run whose projection undershoots
  now stops, and the real fix — a positivity-preserving projection — is the manuscripts' own open
  problem.
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

### K10 · No test covers any structure-preservation claim for Landau.

- **location:** `scripts/verify_conservation.jl`
- **evidence:** `git grep -in 'landau' scripts/verify_conservation.jl` returns nothing, so the
  script exercises the conservative Lenard-Bernstein operator only. The one `Landau` testset under
  `test/` is `test/distributions/spline_distribution.jl:277`, "Landau: K, J and L assemble": it
  asserts assembly, not a conservation invariant.
- **kind:** missing test
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K11 · The cumulant-scaling experiment is not implemented.

- **location:** —
- **evidence:** The appendix fixes `A₀ = 0`, `A₁ = 1` by hand; no code does. The name occurs once
  in the package, in the docstring of the rescaled conservative operator, where it says the
  operator is *not* that experiment; no computation of it exists anywhere.
- **kind:** missing test
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K14 · `fatou lint` reports `unused-binding` locals and no `unused-import`.

- **location:** `src/gridbased/moments.jl:18`
- **evidence:** `fatou lint --force-exclude --output concise .` reports three `unused-binding`
  locals — `src/gridbased/moments.jl:18`, `test/distributions/spline_distribution.jl:74` and
  `:76` — and no `unused-import`.
  `ExplicitImports` agrees with the imports the module file keeps.
- **kind:** dead code
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

### K23 · The one-argument `update_potential!(model)` has no caller in `src/`.

- **location:** `src/models/vlasov_poisson.jl:18`
- **evidence:** It deposits from `model.distribution`, which the integrator never writes. The
  splitting fields call the three-argument method at `:70` and `:92`, and the Hamiltonian at
  `:131`. `grep -rn 'update_potential!' src test` finds the one-argument call only at
  `test/integration/vlasov_poisson.jl:142`, where `deposit(model, x₀)` at `:147` builds the
  same reference.
- **kind:** dead code
- **found:** #61

### K24 · The Lenard-Bernstein right-hand sides keep an indirection and two plotting copies.

- **location:** `src/models/lenard_bernstein.jl:33`
- **evidence:** `LB_rhs!(v̇, v, params, t)` keeps the argument order of an ODE solver that the
  package does not use; its one caller is the wrapper `LB_rhs_GI!` at `:44-46`.
  `LB_rhs` (`:49`) and `CLB_rhs` (`src/models/lenard_bernstein_conservative.jl:184`) repeat the
  bodies of their `_GI!` functions, and neither has a caller: `grep -rn -E '\bC?LB_rhs\b' src
  test scripts` finds them at those two definitions alone, the `LB_rhs!` definition of
  `src/models/lenard_bernstein.jl:33` with its one call at `:45`, and the `LB_rhs!` call of the
  collision-operator test and of its reference generator.
- **kind:** dead code
- **found:** #61

### K26 · The `LB_rhs!` docstring and comment describe past code.

- **location:** `src/models/lenard_bernstein.jl:26`
- **evidence:** The admonition "The division by `f_s` was missing" (`:26-31`) says what "both
  right-hand sides computed", and only one right-hand side is in the file. The comment at
  `:34-35` explains a line that "threw". The facts belong in the present tense, and the history
  in `CHANGELOG.md`.
- **kind:** docs
- **found:** #61

### K27 · The collision stencils wrap the bounded `v`-grid periodically.

- **location:** `src/gridbased/collisions.jl:66`
- **evidence:** The four `v`-stencils — `QuadraticCollisions`' call operator (`:66-67`),
  `ReducedCollisionTensor`'s `getindex` (`:132-133`) and the two `_get_MC̃_*` assemblers
  (`:187-188`, `:237-238`) — wrap with `mod1`, but the `v`-grid is bounded. At the wrap the second
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

- **location:** `src/gridbased/collisions.jl:170`
- **evidence:** `h₁` is never read in either assembler (`:170-205`, `:222-251`). `ci` is read only
  for `size(ci)` (`:171`, `:223`) and for the length check of `v`. The quadratic assembler
  does not read `∫vdv` either. `grep -rn '_get_MC̃_' src test scripts` finds the calls in
  `src/gridbased/collisions.jl` and `test/gridbased/collisions.jl` only. The nine-argument
  signature follows the advisor decision recorded in #62. Whether to shorten it is open.
- **kind:** dead code
- **found:** #62

### K30 · `size(rt, i)` and `axes(rt, i)` throw for `i > 3` on `ReducedCollisionTensor`.

- **location:** `src/gridbased/collisions.jl:103`
- **evidence:** `size(rt, i) = size(rt)[i]` (`:103`) and `axes(rt, i) = Base.OneTo(size(rt, i))`
  (`:104`) index past a three-element tuple, so `size(rt, 4)` and `axes(rt, 4)` throw a
  `BoundsError`. The `AbstractArray` fallbacks give `1` and `Base.OneTo(1)`. The fix is to delete
  the two lines.
- **kind:** defect
- **found:** #62

### K31 · `_stencil_indices_v` has no caller.

- **location:** `src/gridbased/collisions.jl:109`
- **evidence:** `grep -rn '_stencil_indices_v' src test scripts docs` returns the definition
  (`:109-118`) only.
- **kind:** dead code
- **found:** #62

### K32 · The `getindex` of `ReducedCollisionTensor` rebuilds three indices in the innermost loop.

- **location:** `src/gridbased/collisions.jl:120`
- **evidence:** The loop at `:130-148` builds `M`, `N`, `O` and their linear indices for every
  `(m1, m2, o2, n2)`, although `M`, `O`, `m` and `o` do not depend on `n2`, and the `v`-stencil
  of `m2` is the one that `QuadraticCollisions` computes (`:66-67`). The loop can be shorter. No
  benchmark compares the loop with a shorter form.
- **kind:** not verified
- **found:** #62

### K33 · A comment names something a reader cannot find.

- **location:** `src/gridbased/collisions.jl:2`
- **evidence:** The comment at `:2` says "A N × N × N × tensor" and does not name `N`.
- **kind:** docs
- **found:** #62

### K34 · Two testsets check words of a docstring.

- **location:** `test/gridbased/collisions.jl:267`
- **evidence:** `the docstrings name the invariants` (`:267-274`) and
  `the docstrings state the v-end requirement` (`:276-281`) test `occursin` on the text of
  `@doc`, so a rewording that keeps the facts fails them and a docstring that names the words
  but is wrong passes them. The plan of the part names both checks. Whether to keep them is open.
- **kind:** missing test
- **found:** #62

### K36 · A metriplectic Picard step whose particles leave the velocity domain stops with the model's `DomainError`.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:422`
- **evidence:** With `ν = 1e6` at `ti = 4`, all 64 particles leave the `-2.0 .. 2.0` velocity
  support and `projection` throws `DomainError: … 64 of 64 particles left the velocity domain …`
  from inside the solve. The catch at `:422` rethrows it, so the caller sees that `DomainError`
  and not the residual-and-count `ErrorException` the author decided for a solve that does not
  converge or meets a `NaN`. The error is loud, so the solve never returns in silence.
- **kind:** found late
- **found:** 2026-10-04

### K41 · `Picard_iterate_over_particles!` calls two `SimpleSolvers` names that are not public.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:423`
- **evidence:** `SimpleSolvers.status` (`:423`) and `SimpleSolvers.isconverged` (`:426`) are neither
  exported nor declared `public`; `SimpleSolvers` 0.14.1 keeps them unexported at
  `src/SimpleSolvers.jl:138`. `ExplicitImports.check_all_qualified_accesses_are_public` flags
  exactly these two. The fix is a `public` declaration in `SimpleSolvers`.
- **kind:** upstream
- **found:** 2026-10-04

### K42 · `Picard_iterate_over_particles` takes `dv`, `m` and `β` and ignores them.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:443`
- **evidence:** The wrapper does not read `dv`, `m` or `β` after its signature (`:443-446`); it
  passes the other arguments to `Picard_iterate_over_particles!`, which takes none of the three.
  They stay so that the call sites do not change, and the CHANGELOG says so. Removing them is an
  API change.
- **kind:** dead code
- **found:** 2026-10-04

### K43 · `Picard_iterate_Landau_nls!` prints from library code.

- **location:** `src/methods/Landau_solver.jl:35`
- **evidence:** The calls at `:35`, `:47` and `:55` print the residual of each iteration and an
  empty line. They are the only residual report of the function (K5), so removing them removes
  the one output that shows the drift.
- **kind:** defect
- **found:** 2026-10-04

### K44 · The reference generator of the metriplectic test is in `test/helpers/`.

- **location:** `test/helpers/generate_metriplectic_reference.jl`
- **evidence:** The file is a script that a person runs to write the reference data, and not a
  helper that a test loads. A script of this kind belongs in `scripts/`.
  `test-layout.jl --check` passes with the file where it is, and
  `test/integration/metriplectic_solve.jl` names this path.
- **kind:** docs
- **found:** 2026-10-04

### K45 · Three collision caches have no `eltype`, so `CacheDict` builds their first cache twice.

- **location:** `src/cache.jl:12`
- **evidence:** `CacheDict(p)` stores its parent under `_cachehash(eltype(p))`. Only `LandauCache`
  defines `Base.eltype` (`src/models/landau.jl:53`); for `CLBCache`, `MLBCache` and `RCLBCache`,
  `which(eltype, Tuple{C}).module` is `Base`, whose fallback returns `Any`. The first
  `cache[Float64]` then misses and builds a new cache with `Cache(Float64, parent)`. The result is
  correct, but the gap decides which spline an operator overwrites. Landau's `cache[Float64]` is
  the parent, so it projects into `entropy.dist`. The three Lenard–Bernstein operators project
  into the new cache's own copy of the spline. An `eltype` for the three caches therefore changes
  what they write, not only how often a cache is built: for each, `model.cache[Float64] ===
  parent(model.cache)` is `false`, and for Landau it is `true`.
- **kind:** defect
- **found:** 2026-10-10

### K47 · Thirteen lines of commented-out code remain in four collision model files.

- **location:** `src/models/landau.jl:385`
- **evidence:** The lines are `src/models/landau.jl:385`; `src/models/lenard_bernstein.jl:57`
  and `:67`; `src/models/lenard_bernstein_conservative.jl:196`, `:206-208` and `:212`;
  `src/models/rescaled_lenard_bernstein_conservative.jl:125`, `:138-140` and `:144`. Each is an
  old assignment, call or `return` behind `#`, such as `# int = Integrators.Integrator(equ,
  Integrators.RK438())`. #69 deleted the other dead code of these files and left these lines,
  because they are outside its hunks. The fix is to delete the thirteen lines.
- **kind:** dead code
- **found:** #69

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

- **location:** `test/integration/collision_operators.jl:94`
- **evidence:** `f!` has one method, `src/models/lenard_bernstein_metriplectic.jl:339`, and its
  last argument is `mlb::MetriplecticLenardBernstein`. So `which(VlasovMethods.f!, …)` at
  `:94-96` either returns that method, and the `@test` at `:97` passes, or throws a
  `MethodError`. The value pin at `:90` already runs `f!` through the Picard step. Whether to
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

- **location:** `src/projections/distribution.jl:188`
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

### K55 · The manual's navbar carries no repository link.

- **location:** `docs/make.jl:9`
- **evidence:** Every build of the manual warns
  `┌ Warning: Unable to determine the repository root URL for the navbar link.`, and the same
  message names the two remedies: pass a `Remotes.Remote` object to `repo`, or set `repolink` on
  the `HTML` format. The `makedocs` call passes the string
  `"https://github.com/JuliaPlasma/VlasovMethods.jl/blob/{commit}{path}#{line}"`, which serves the
  source links but not the navbar.
- **kind:** docs
- **found:** #73

### K56 · `Parameter` of `ReducedComplexityModeling` rejects distinct bounds for a one-sample range.

- **location:** `ReducedComplexityModeling` v0.1.1, `src/parameters/parameter.jl:33`
- **evidence:** `Parameter(name, minimum::DT, maximum::DT, n::Int)` builds
  `LinRange(minimum, maximum, n)`, which throws for `n == 1` when the endpoints differ:
  `julia --startup-file=no --project -e 'using ReducedComplexityModeling; Parameter(:β, 5.95, 6.05, 1)'`
  prints `ArgumentError: range(5.95, stop=6.05, length=1): endpoints differ` on the RCM 0.1.1 this
  package resolves, while `Parameter(:β, 5.95, 6.05, 2)` gives
  `Parameter{Float64, Vector{Float64}}(:β, 5.95, 6.05, [5.95, 6.05])`. A one-sample sweep is
  supported through `Parameter(:β, [6.0])` or equal bounds, `Parameter(:β, 6.0, 6.0, 1)`:
  each gives a `ParameterSpace` of length 1 and `size(Snapshots(Float64, 2, 100, 16, 11, ps).X)`
  is `(2, 100, 11, 1)`.
- **kind:** upstream
- **found:** #75

### K57 · `plot_distribution` checks neither the orientation of the two domains nor the cell counts.

- **location:** `ext/VlasovMethodsMakieExt.jl:19-28`
- **evidence:** A reversed `xdomain` gives a histogram with a negative "density" and no error, and
  `nx = 0` fails with an indexing error; every other edge case of the two methods raises
  `ArgumentError`. Probe, CairoMakie 0.15.15 on Julia 1.13.1 (the review of #74, finding 3):
  `xdomain = (1.0, 0.0)` with one particle of weight 1 and `nx = nv = 4` returns a density of
  `-8.0`; `nx = 0` throws `BoundsError: attempt to access 0×100 Matrix{Float64} at index [0, 51]`.
  One guard covers both cases: `nx ≥ 1 && nv ≥ 1 && xa < xb && va < vb || throw(ArgumentError(…))`.
- **kind:** defect
- **found:** #74

### K58 · A particle whose `v` is `NaN` is dropped from the histogram without a word.

- **location:** `src/plotting.jl:34-36`
- **evidence:** A particle whose `v` is `NaN` is also dropped, because `va ≤ v` is false for it.
  Neither the docstring nor the CHANGELOG mentions this case; both name only `v` outside
  `[va, vb)` and a non-finite `x`. Probe (the review of #74, finding 6): particles `(0.5, NaN,
  w = 1)` and `(0.5, 0.0, w = 0.25)`, and `sum(F)·hx·hv == 0.25` passes.
- **kind:** defect
- **found:** #74

### K59 · The grid heatmap's docstring says the heatmap is drawn "on the grid's own nodes".

- **location:** `src/plotting.jl:19-20`
- **evidence:** Makie treats the nodes as cell centres, so the drawn `v` range extends `hv/2` past
  `vdomain` at each end and every `x` cell is centred on its node. Probe (the review of #74,
  finding 5): `GridDistribution(reshape(1:21, 7, 3), (0, 1), (-1, 1))` draws `v` cell edges
  `[-1.5, -0.5, 0.5, 1.5]` and `x` edges `[-0.0714…, …, 0.9286…]`. `plot.args[]` equals
  `(dist.x, dist.v)`, as the test asserts.
- **kind:** docs
- **found:** #74

### K60 · The CairoMakie extension relies on implicit imports.

- **location:** `ext/VlasovMethodsMakieExt.jl:3-4`
- **evidence:** `print_explicit_imports(VlasovMethodsMakieExt; report_non_public = true)` reports
  `relying on implicit imports for 8 names`: `using CairoMakie: CairoMakie, Axis, Figure,
  heatmap!, scatter!` and `using VlasovMethods: VlasovMethods, GridDistribution,
  ParticleDistribution`. GeometricBrackets' extension uses the same pattern.
- **kind:** defect
- **found:** #74

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

### K62 · A comment says that the clamped cubic basis reproduces 1, v and v² to 1e-14, and v² is measured at 2.8e-14.

- **location:** `src/models/lenard_bernstein_conservative.jl:39`
- **evidence:** `julia --startup-file=no --project=scripts scripts/verify_conservation.jl` measures
  the cubic basis on `(-8.0, 8.0)` with 31 knots and prints, under Claim 2,
  `Free()  3  6.661e-16  1.776e-15  2.842e-14` for `f = 1`, `f = v` and `f = v²`. The other
  figures of the comment agree with the same table: `Periodic()` gives `1.234e+00` for `v` and
  `1.164e-01` for `v²`, and `Dirichlet()` gives `4.066e-02` for `1`. The CHANGELOG states the
  bound as `3e-14`.
- **kind:** docs
- **found:** #80

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

### K68 · No 1D1V distribution matches the `draw!` of the `Bump` example.

- **location:** `src/examples/bump.jl:28`
- **evidence:** `initialize!(ParticleDistribution(1, 1, 100), VlasovMethods.Bump())` raises
  `MethodError: no method matching draw!(::ParticleDistribution{Float64, 1, 1, …}`, on
  `05431f10` and on this branch. The method binds `ParticleDistribution{1, 1}`, whose first
  parameter is the element type.
- **kind:** found late
- **found:** 2026-10-23

### K69 · The legend of `plot_series` overflows its axis for many labelled series.

- **location:** `ext/VlasovMethodsMakieExt.jl:128`, `plot_series`
- **evidence:** at the default figure size of `600 × 450`, the legend of 20 labelled series is
  481 px high in an axis of 375 px, on `05431f10`; 10 series give 251 px.
  `scripts/bump_on_tail_testing.jl` draws its 30-series figures of items 26 and 27 without
  labels and adds a legend of one entry per model. Its figures of items 24 and 25 label one series
  per sample, so at the default of 10 samples no figure of the scripts shows the overflow.
- **kind:** found late
- **found:** 2026-10-23

### K74 · The metriplectic Picard step overwrites a `vn` that is a view of `dv_history[:, 1]`.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:431`, `Picard_iterate_over_particles!`
- **evidence:** a `ti = 4` call of `Picard_iterate_over_particles` with
  `vn = view(dv_history, :, 1)`, 64 particles on `range(-1.99, 1.99)`, returns, and the history
  update after the solve changes the caller's `vn` by up to `6.628111373904642`, on `3175a68`
  and on this branch. The step refuses an output `v` that aliases `vn` or `dv_history`, and does
  not check `vn` against `dv_history`.
- **kind:** found late
- **found:** 2026-10-24
