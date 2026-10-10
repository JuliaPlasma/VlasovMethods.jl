# Known issues

What is known to be wrong in VlasovMethods and is not fixed. An entry leaves this file when its
fix merges, and the CHANGELOG entry of that fix names its ID. IDs are never reused; a new entry
takes the next `K<n>`.

### K1 · `src/particles/` is present but not included.

- **location:** `src/particles/`
- **evidence:** The four files imported from ReducedBasisMethods — `electric_field.jl`,
  `poisson.jl`, `snapshots.jl`, `time_marching.jl` — each name a binding that no longer exists, so
  including any one breaks the load: `poisson.jl` wants `PoissonSolverPBSplines`;
  `time_marching.jl` imports `PBSpline`, `stiffnessmatrix`, `eval_deriv_PBSBasis` and
  `rhs_particles_PBSBasis` from `PoissonSolvers` (neither 0.5 nor 0.6 has any);
  `electric_field.jl` wants `ElectricField`, which nothing under `src/` defines, and
  `snapshots.jl` wants `ParameterSpace`. They were moved unrepaired on
  purpose, so the relocation stays reviewable.
- **kind:** dead code
- **found:** 2026-09-17. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K3 · `FullyReducedTensor` cannot be constructed.

- **location:** —
- **evidence:** Its inner constructor asserts `size(Pk, 1) == size(tensor, 3)`, but the parameter
  is named `Pα`, so `Pk` is undefined; and its `getindex` reads `rt.projection_k[k, α]` with `α`
  unbound. Imported unrepaired from ReducedBasisMethods, where it had the same defects.
- **kind:** defect
- **found:** 2026-09-17. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

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

- **location:** `scripts/lenard_bernstein_conservative.jl`
- **evidence:** The appendix fixes `A₀ = 0`, `A₁ = 1` by hand; no code does. The cumulant
  computation is commented out in `scripts/lenard_bernstein_conservative.jl` and the experiment
  survives as a stale `run_name`.
- **kind:** missing test
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K14 · `fatou lint` reports `unused-binding` locals and no `unused-import`.

- **location:** `scripts/bump_on_tail.jl:47`
- **evidence:** `fatou lint --force-exclude --output concise .` reports seven `unused-binding`
  locals — `scripts/bump_on_tail.jl:47`, `scripts/lenard_bernstein_metriplectic_scaling.jl:38`
  and `:49`, `src/gridbased/moments.jl:18`, `src/gridbased/reduced_tensors.jl:210`,
  `test/distributions/spline_distribution.jl:74` and `:76` — and no `unused-import`.
  `ExplicitImports` agrees with the imports the module file keeps.
- **kind:** dead code
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K15 · `size(rt, i)` and `axes(rt, i)` throw for `i > ndims(rt)` on three reduced tensors.

- **location:** `src/gridbased/reduced_tensors.jl:82`
- **evidence:** `PotentialReducedTensor` (`:82-83`), `VelocityReducedMatrix` (`:137-138`) and
  `FullyReducedTensor` (`:188-189`) define `size(rt, i) = size(rt)[i]` and
  `axes(rt, i) = Base.OneTo(size(rt, i))`, which index past the tuple:
  `size(PotentialReducedTensor(t, Pi, Pj, Pk), 4)`, `axes(…, 4)`,
  `size(VelocityReducedMatrix(t, Pi, Pj, v²), 3)` and `axes(…, 3)` each throw a `BoundsError`.
  The `AbstractArray` fallbacks give `1` and `Base.OneTo(1)`, as `size(rt, 4)` and `axes(rt, 4)`
  do for `ReducedTensor`, which defines neither. The fix is to delete the six lines.
- **kind:** defect
- **found:** #51

### K16 · A commented-out `Base.materialize` for `ReducedTensor` and `PotentialReducedTensor` remains.

- **location:** `src/gridbased/reduced_tensors.jl:113-115`
- **evidence:** The three lines are comments, so no method exists; `collect(rt)` and `Array(rt)`
  already materialise any of these `AbstractArray`s. The fix is to delete the three lines.
- **kind:** dead code
- **found:** #51

### K17 · `ReducedTensor` accepts a tensor whose grid is not the grid of its `Arakawa`.

- **location:** `src/gridbased/reduced_tensors.jl:24`
- **evidence:** `PoissonTensor(DT, nx, nv, f)` in GeometricBrackets does not check that `f` is on
  the `nx × nv` grid, and the `ReducedTensor` constructor does not check it either. The 3 × 3
  stencil then loses coefficients or counts one twice, and no error occurs:
  `ReducedTensor(PoissonTensor(Float64, 2, 3, Arakawa(3, 3, 0.5, 0.5)), Pi, Pj)` deviates from the
  dense sum by up to 0.50, and `PoissonTensor(Float64, 5, 4, Arakawa(3, 3, …))` by up to 2.51. The
  fix is a grid check in the `PoissonTensor` constructor of GeometricBrackets.
- **kind:** upstream
- **found:** #53

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

### K19 · `VlasovPoisson` cannot be built with a grid `Potential` basis.

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

### K20 · `scripts/lenard_bernstein.jl` calls a constructor the package does not define.

- **location:** `scripts/lenard_bernstein.jl:30`
- **evidence:** Line `:30` calls `DiffEqIntegrator(model, tspan, tstep)`, a type the package does
  not define, so the script throws `UndefVarError` when it reaches it. Line `:37` then calls
  `VlasovMethods.run(integrator)`, which resolves to `Base.run`: the package defines no `run`
  method. The `scripts/` rewrite repairs or removes both calls.
- **kind:** defect
- **found:** 2026-10-03

### K22 · `scripts/bump_on_tail.jl` calls two names the package does not define.

- **location:** `scripts/bump_on_tail.jl:31`
- **evidence:** `:31` calls `VPIntegratorParameters(dt, nₜ, nₜ+1, nₕ, nₚ)` and `:58` calls
  `integrate_vp!(P, efield, params, IP, IC)`. Both were defined only in the root
  `src/vlasov_poisson.jl`, which no `include` reached, so neither name exists in the package.
  `grep -rn -E 'VPIntegratorParameters|integrate_vp!' src` returns nothing.
- **kind:** defect
- **found:** #61

### K23 · The one-argument `update_potential!(model)` has no caller in `src/`.

- **location:** `src/models/vlasov_poisson.jl:18`
- **evidence:** It deposits from `model.distribution`, which the integrator never writes. The
  splitting fields call the two-argument method at `:65` and `:83`.
  `grep -rn 'update_potential!' src test` finds the one-argument call only at
  `test/integration/vlasov_poisson.jl:122`, where `deposit(model, x₀)` at `:127` builds the
  same reference.
- **kind:** dead code
- **found:** #61

### K24 · The Lenard-Bernstein right-hand sides keep an indirection and two plotting copies.

- **location:** `src/models/lenard_bernstein.jl:33`
- **evidence:** `LB_rhs!(v̇, v, params, t)` keeps the argument order of an ODE solver that the
  package does not use; its one caller is the wrapper `LB_rhs_GI!` at `:44-46`.
  `LB_rhs` (`:49`) and `CLB_rhs` (`src/models/lenard_bernstein_conservative.jl:184`) repeat the
  bodies of their `_GI!` functions. `grep -rn -E '\bC?LB_rhs\b' src test scripts` finds
  `LB_rhs` only at `scripts/lenard_bernstein.jl:47` and `:57`, after that script fails at `:30`
  (K20), and `CLB_rhs` only in commented script lines.
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

### K33 · Two comments and the CHANGELOG name something a reader cannot find.

- **location:** `src/gridbased/collisions.jl:2`
- **evidence:** The comment at `:2` says "A N × N × N × tensor" and does not name `N`.
  `CHANGELOG.md:226` says the open items "are recorded under *Open Issues*", and the CHANGELOG
  has no such heading: they are in `KNOWN_ISSUES.md`.
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

- **location:** `src/models/lenard_bernstein_metriplectic.jl:300`
- **evidence:** With `ν = 1e6` at `ti = 4`, all 64 particles leave the `-2.0 .. 2.0` velocity
  support and `projection` throws `DomainError: … 64 of 64 particles left the velocity domain …`
  from inside the solve. The catch at `:300` rethrows it, so the caller sees that `DomainError`
  and not the residual-and-count `ErrorException` the author decided for a solve that does not
  converge or meets a `NaN`. The error is loud, so the solve never returns in silence.
- **kind:** found late
- **found:** 2026-10-04

### K37 · The two metriplectic scripts read a solve-result object the solve does not return.

- **location:** `scripts/lenard_bernstein_metriplectic.jl:81`
- **evidence:** `Picard_iterate_over_particles` returns the solved velocity `Vector`, but
  `scripts/lenard_bernstein_metriplectic.jl` still reads `sol_object.u` (`:81`, `:91`) and calls
  `SciMLBase.successful_retcode(sol_object)` (`:84`), and
  `scripts/lenard_bernstein_metriplectic_scaling.jl` reads `sol_object.u` (`:96`, `:107`); each
  throws once it reaches those lines. `scripts/lenard_bernstein_metriplectic.jl:70` also sets
  `abstol = 1e-15`, which the solve now compares with `‖F‖₂` and not `maximum(abs, F)`. That is
  below the round-off floor of `‖F‖₂`: one step of the uniform cloud of
  `test/integration/metriplectic_solve.jl` with `abstol = 1e-15` throws at `N = 64`, `256` and
  `1000` (`residual = 1.21e-15`, `1.26e-15`, `1.74e-15` after 49, 29 and 26 iterations), and
  returns with the scaling script's `3e-16·√N`. The `scripts/` rewrite (P46) repairs both.
- **kind:** defect
- **found:** 2026-10-04

### K40 · Each metriplectic Picard step builds a fresh `SimpleSolvers` solver and its unused Jacobian cache.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:293`
- **evidence:** `Picard_iterate_over_particles` constructs a fresh
  `SimpleSolvers.NonlinearSolver(Picard(), …)` on every call. The constructor allocates an `N×N`
  `solver.cache.j` and a `ForwardDiff.JacobianConfig` the unaccelerated Picard step never reads:
  measured `126288 B` at `N = 64` and `96629344 B` at `N = 2000` per step. The solver's type also
  leaves the chunk size open, so `solve_with_status!` is one dynamic dispatch per call. The return
  type is concrete and `@inferred` passes. Avoiding the allocation needs the solver to be reused
  across steps (an API change) or an upstream change to the Picard constructor's cache.
- **kind:** defect
- **found:** 2026-10-04

### K41 · `Picard_iterate_over_particles` calls two `SimpleSolvers` names that are not public.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:301`
- **evidence:** `SimpleSolvers.status` (`:301`) and `SimpleSolvers.isconverged` (`:304`) are neither
  exported nor declared `public`; `SimpleSolvers` 0.14.1 keeps them unexported at
  `src/SimpleSolvers.jl:138`. `ExplicitImports.check_all_qualified_accesses_are_public` flags
  exactly these two. The fix is a `public` declaration in `SimpleSolvers`.
- **kind:** upstream
- **found:** 2026-10-04

### K42 · `Picard_iterate_over_particles` takes `dv`, `m` and `β` and ignores them.

- **location:** `src/models/lenard_bernstein_metriplectic.jl:253`
- **evidence:** The function body does not read `dv`, `m` or `β` after the signature (`:253-256`).
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
