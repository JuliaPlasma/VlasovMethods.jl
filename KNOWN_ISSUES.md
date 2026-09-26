# Known issues

What is known to be wrong in VlasovMethods and is not fixed. An entry leaves this file when its
fix merges, and the CHANGELOG entry of that fix names its ID. IDs are never reused; a new entry
takes the next `K<n>`.

## K1

- **location:** `src/particles/electric_field.jl`, `src/particles/poisson.jl`,
  `src/particles/snapshots.jl`, `src/particles/time_marching.jl`
- **problem:** `src/particles/` is present but not included.
- **evidence:** The four files imported from ReducedBasisMethods — `electric_field.jl`,
  `poisson.jl`, `snapshots.jl`, `time_marching.jl` — each name a binding that no longer exists, so
  including any one breaks the load: `poisson.jl` wants `PoissonSolverPBSplines`;
  `time_marching.jl` imports `PBSpline`, `stiffnessmatrix`, `eval_deriv_PBSBasis` and
  `rhs_particles_PBSBasis` from `PoissonSolvers` (neither 0.5 nor 0.6 has any);
  `electric_field.jl` wants `ElectricField` from `src/electric_field.jl` — which this module still
  keeps commented out — and `snapshots.jl` wants `ParameterSpace`. They were moved unrepaired on
  purpose, so the relocation stays reviewable.
- **kind:** dead code
- **found:** 2026-09-17, in the audit that accompanied the `SimpleSplines` migration

## K2

- **location:** `src/gridbased/collisions.jl:33`, `src/gridbased/collisions.jl:46`,
  `src/gridbased/collisions.jl:139`
- **problem:** `src/gridbased/collisions.jl` compiles but cannot be used.
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
- **found:** 2026-09-17, in the audit that accompanied the `SimpleSplines` migration

## K3

- **location:** `src/gridbased/reduced_tensors.jl:122`, `src/gridbased/reduced_tensors.jl:154`
- **problem:** `FullyReducedTensor` cannot be constructed.
- **evidence:** Its inner constructor asserts `size(Pk, 1) == size(tensor, 3)`, but the parameter
  is named `Pα`, so `Pk` is undefined; and its `getindex` reads `rt.projection_k[k, α]` with `α`
  unbound. Imported unrepaired from ReducedBasisMethods, where it had the same defects.
- **kind:** defect
- **found:** 2026-09-17, in the audit that accompanied the `SimpleSplines` migration

## K4

- **location:** `src/models/landau.jl`
- **problem:** The implemented Landau scheme is not the one the main text derives.
- **evidence:** The manuscript builds the gradient form with the `G` operator — whose structure
  *is* the momentum and energy conservation proof — and a Gonzalez discrete gradient, which *is*
  the discrete H-theorem proof. What runs is the **appendix** two-step `v̇ = K⁺LJ` with plain
  implicit midpoint and `∇S(midpoint)`, which is not a discrete gradient. Neither structural proof
  transfers to the code as written. `G` is never formed; the gradient form survives only as a
  commented-out `Landau_rhs`. A gap between paper and code, not an error in either.
- **kind:** docs
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K5

- **location:** `src/methods/Landau_solver.jl:152`, `src/methods/Landau_solver.jl:124`
- **problem:** The Landau Picard solver does not iterate to convergence.
- **evidence:** `Landau_solver.jl` runs exactly five iterations and prints the residual without
  testing it. `tol`, `ftol`, `β`, `m` and `chunksize` are accepted and unused, and `probN` is
  constructed and never solved — it is left in place, with the commented-out `NonlinearSolve`
  calls it belongs to, rather than deleted. Every conservation property in the appendix is a
  property of the *exactly* solved implicit system, so momentum and energy drift at the size of
  that printed residual.
- **kind:** defect
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K6

- **location:** `src/methods/Landau_solver.jl`
- **problem:** The accepted Landau solution is one Picard step behind its stored derivative.
- **evidence:** The loop's last action recomputes `v̇` at the midpoint from the newest guess
  without updating the guess, so the stored state and derivative do not correspond and the next
  step's Hermite extrapolation is fed an inconsistent pair.
- **kind:** defect
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K7

- **location:** `src/models/landau.jl`
- **problem:** The Landau kernel's coincident-point value is a regularisation, not a limit.
- **evidence:** `kernel` returns zero at `|u| = 0` so that a product quadrature sharing nodes does
  not produce `Inf`. The error does not vanish under refinement, and because both factors use the
  same Gauß–Legendre nodes it fires on every *diagonal cell pair* rather than on a set of measure
  zero. The kernel is integrable in two dimensions; what is needed is a singularity-aware rule or
  offset grids.
- **kind:** defect
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K8

- **location:** `src/models/landau.jl`
- **problem:** Positivity of `f_s` is detected, not solved.
- **evidence:** Every `log f_s` and `1/f_s` now throws rather than continuing, which is strictly
  better than the old `0.5·log(f_s²)` returning `log|f_s|`. But a run whose projection undershoots
  now stops, and the real fix — a positivity-preserving projection — is the manuscripts' own open
  problem.
- **kind:** defect
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K9

- **location:** `src/models/landau.jl:546`
- **problem:** The rank hypothesis behind `K K⁺ = I` is no longer checked anywhere.
- **evidence:** The step to `eq:particle_ode_final` needs `K` to have full row rank. The old code
  detected the failure with two SVDs per vector-field evaluation, printed, and proceeded
  regardless; the check is gone from the hot path and has not been given a home in a diagnostic
  script.
- **kind:** missing test
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K10

- **location:** `test/runtests.jl:16`, `scripts/verify_conservation.jl`
- **problem:** No test covers any structure-preservation claim for Landau.
- **evidence:** `scripts/verify_conservation.jl` covers the conservative Lenard-Bernstein operator
  only, and `test/electric_field_tests.jl` remains commented out of `runtests.jl`.
- **kind:** missing test
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K11

- **location:** `scripts/lenard_bernstein_conservative.jl:9`
- **problem:** The cumulant-scaling experiment is not implemented.
- **evidence:** The appendix fixes `A₀ = 0`, `A₁ = 1` by hand; no code does. The cumulant
  computation is commented out in `scripts/lenard_bernstein_conservative.jl` and the experiment
  survives as a stale `run_name`.
- **kind:** missing test
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K12

- **location:** `src/VlasovMethods.jl:170`, `src/VlasovMethods.jl:175`, `src/VlasovMethods.jl:183`
- **problem:** Seven `src/` files are included by nothing.
- **evidence:** `electric_field.jl`, the root `vlasov_poisson.jl` (distinct from
  `models/vlasov_poisson.jl`), `visualisation.jl`, `methods/lbm_solver.jl`, `hdf5.jl`, and two
  whose `include` lines are commented out. `src/hdf5.jl` is untracked in git.
- **kind:** not verified
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K13

- **location:** `Project.toml`
- **problem:** Dependencies that are no longer used are still declared.
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
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration

## K14

- **location:** `src/VlasovMethods.jl`, `src/methods/Landau_solver.jl:124`
- **problem:** `fatou lint` reports 10 warnings, of which one is deliberate and nine are a known
  false positive.
- **evidence:** The nine are `unused-import` on `src/VlasovMethods.jl`, where the rule does not
  follow `include` and so flags the module file's load-bearing imports; `ExplicitImports`
  contradicts all nine. The tenth is `probN` below. An earlier version of this changelog and of the
  pull request described `fatou lint` as clean, which was not reproducible.
- **kind:** upstream
- **found:** 2026-09-07, in the audit that accompanied the `SimpleSplines` migration
