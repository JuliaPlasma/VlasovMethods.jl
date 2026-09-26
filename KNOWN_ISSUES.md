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
  `electric_field.jl` wants `ElectricField` from `src/electric_field.jl` — which this module still
  keeps commented out — and `snapshots.jl` wants `ParameterSpace`. They were moved unrepaired on
  purpose, so the relocation stays reviewable.
- **kind:** dead code
- **found:** 2026-09-17. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K2 · `src/gridbased/collisions.jl` compiles but cannot be used.

- **location:** `src/gridbased/collisions.jl`
- **evidence:** The file is included and the module precompiles, so `CollisionTensor` and the rest
  are defined — the faults are at run time, not load time.
  `Base.getindex(ct::CollisionTensor, i, j, k)` returns `ct[I, J, K, L]` with `L` never bound, and
  both `_get_MC̃_*` assemblers read a global `v` that no longer exists. `ReducedCollisionTensor`'s
  `getindex` evaluates `rt.tensor[M, N, O]` with three `CartesianIndex{2}` arguments;
  `Base.to_indices` flattens them to six `Int`s, so the indexing throws `BoundsError` and never
  reaches the four-`CartesianIndex` method of `CollisionTensor` — binding `L` alone does not repair
  it. The `QuadraticCollisions` inner constructor calls `new{DT}(nx, nv, hx, hv, v)` — five values
  for six fields — leaving `factor` uninitialised. ReducedBasisMethods never included the file, so
  none of this was reachable there and none of it is new.
- **kind:** defect
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
  transfers to the code as written. `G` is never formed; the gradient form survives only as a
  commented-out `Landau_rhs`. A gap between paper and code, not an error in either.
- **kind:** docs
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K5 · The Landau Picard solver does not iterate to convergence.

- **location:** `Landau_solver.jl`
- **evidence:** `Landau_solver.jl` runs exactly five iterations and prints the residual without
  testing it. `tol`, `ftol`, `β`, `m` and `chunksize` are accepted and unused, and `probN` is
  constructed and never solved — it is left in place, with the commented-out `NonlinearSolve`
  calls it belongs to, rather than deleted. Every conservation property in the appendix is a
  property of the *exactly* solved implicit system, so momentum and energy drift at the size of
  that printed residual.
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
- **evidence:** `scripts/verify_conservation.jl` covers the conservative Lenard-Bernstein operator
  only, and `test/electric_field_tests.jl` remains commented out of `runtests.jl`.
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

### K12 · Seven `src/` files are included by nothing:

- **location:** `src/`
- **evidence:** `electric_field.jl`, the root `vlasov_poisson.jl` (distinct from
  `models/vlasov_poisson.jl`), `visualisation.jl`, `methods/lbm_solver.jl`, `hdf5.jl`, and two
  whose `include` lines are commented out. `src/hdf5.jl` is untracked in git.
- **kind:** not verified
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K13 · Dependencies that are no longer used are still declared.

- **location:** `src/`
- **evidence:** `NaNMath` appears only in a commented import; `Plots`, `LaTeXStrings`,
  `StatsPlots`, `StatsBase`, `SciMLBase` and `AdaptiveRejectionSampling` have no occurrences by
  name in `src/`. Five non-stdlib dependencies still carry no `[compat]` entry — `LinearSolve`,
  `NaNMath`, `NonlinearSolve`, `SimpleSolvers` and `Trapz` — and General's AutoMerge blocks on
  every one of them, so this is what stands between the package and registration. Four driver
  scripts also `using` `GLMakie`, `Printf` and `Profile`, none of which are declared.
  `ExplicitImports.jl` *has* now been run and reports no stale or improper explicit imports, so
  what remains is `[deps]` hygiene rather than dead `import` lines — `Aqua.test_stale_deps` is the
  check that settles it.
- **kind:** not verified
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.

### K14 · `fatou lint` reports 10 warnings, of which one is deliberate and nine are a known false positive.

- **location:** `src/VlasovMethods.jl`
- **evidence:** The nine are `unused-import` on `src/VlasovMethods.jl`, where the rule does not
  follow `include` and so flags the module file's load-bearing imports; `ExplicitImports`
  contradicts all nine. The tenth is `probN` below. An earlier version of this changelog and of the
  pull request described `fatou lint` as clean, which was not reproducible.
- **kind:** upstream
- **found:** 2026-09-07. Carried over from the audit that accompanied the `SimpleSplines`
  migration. None of these are regressions; each is either a numerical-methods decision or work
  the migration deliberately did not take on.
