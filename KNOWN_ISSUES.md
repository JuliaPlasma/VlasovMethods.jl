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

- **location:** `src/models/lenard_bernstein.jl:37`
- **evidence:** `LB_rhs!(v̇, v, params, t)` keeps the argument order of an ODE solver that the
  package does not use; its one caller is the wrapper `LB_rhs_GI!` at `:48-50`.
  `LB_rhs` (`:53`) and `CLB_rhs` (`src/models/lenard_bernstein_conservative.jl:210`) repeat the
  bodies of their `_GI!` functions. `grep -rn -E '\bC?LB_rhs\b' src test scripts` finds
  `LB_rhs` only at `scripts/lenard_bernstein.jl:47` and `:57`, after that script fails at `:30`
  (K20), and `CLB_rhs` only in commented script lines.
- **kind:** dead code
- **found:** #61

### K25 · `IM_rule!` in `src/methods/Landau_solver.jl` has no caller.

- **location:** `src/methods/Landau_solver.jl:5`
- **evidence:** `grep -rn 'IM_rule!' src test scripts` finds the definition at `:5-12` and
  otherwise only comments: `Landau_solver.jl:1` and `:16`, and
  `src/models/lenard_bernstein_metriplectic.jl:284`, `:394` and `:445`.
- **kind:** dead code
- **found:** #61

### K26 · The `LB_rhs!` docstring and comment describe past code.

- **location:** `src/models/lenard_bernstein.jl:30`
- **evidence:** The admonition "The division by `f_s` was missing" (`:30-35`) says what "both
  right-hand sides computed", and only one right-hand side is in the file. The comment at
  `:38-39` explains a line that "threw". The facts belong in the present tense, and the history
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
