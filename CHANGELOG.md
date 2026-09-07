# Release Notes

All notable changes to VlasovMethods.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and deliberately holds no entries. 5 versions were
released before it, the most recent `v0.2.1`, and none of them are written up here: the
record of that history is `git log` and the tags. It is named as a gap rather than
reconstructed, because a changelog assembled after the fact loses exactly the reasoning that
makes it worth keeping. The `[Unreleased]` target below is provisional — confirm it when the
first entry is written.

## [Unreleased] — targeting 0.3.0

### Breaking Changes

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

  `SimpleSplines` is not yet registered in General, so it carries no `[compat]` bound and
  **this package cannot be registered until it is**; resolving needs a scratch environment with
  both developed. `FastGaussQuadrature` and `QuadGK` are dropped, having no remaining caller.

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
  asserted in `SimpleSplines`, and `test/spline_distribution_tests.jl` covers the
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

### New Features

- **`scripts/verify_conservation.jl`** measures the two conservation claims of the manuscripts
  and separates them, because they are not the same claim and the obvious reading is wrong.
  The particle sums `Σ w v̇` and `Σ w v v̇` vanish at round-off on **every** boundary condition,
  since `A₁` and `A₂` are solved from exactly those constraints. What the span requirement of
  `main.tex:226` controls is the *other* claim, `∫ vᵏ f_s dv = Σ w v^k`: measured on a cubic
  basis over `[-8,8]`, `Free()` reproduces `1`, `v`, `v²` to `1e-14`, `Periodic()` gets `v`
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

- **A deposition allocation test.** `test/spline_distribution_tests.jl` asserts that the cost of
  `projection` is independent of the particle count, so a boxed closure in `_deposit!` cannot
  come back unnoticed — nothing else in the suite would see it, since the results stay correct
  and only the run time changes. Guarded on `--check-bounds=auto`, because `Pkg.test()`'s
  default `=yes` inflates allocation counts and would make the ceiling meaningless.

### Bug Fixes

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
  same copy-as-output bug; `landau_newer.jl` already used a `view`.

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
  why every entropy computation in `scripts/` is commented out. It now integrates on the basis's
  own Gauß-Legendre grid in any number of velocity dimensions.

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

- **`test/particle_distribution_tests.jl` was flaky at roughly one run in four**, and had never
  run at all: `runtests.jl` included only the two spline test files. It asserts
  `mean(v) ≈ centre atol = 3.5/sqrt(np)`, which for a uniform distribution of width `w` — whose
  standard error of the mean is `w/sqrt(12 np)` — is about 3.0 standard errors, so each of the
  120 such assertions fails about 0.27 % of the time. Adding the file to the suite also made it
  sensitive to whatever RNG state the preceding testset left behind. Now seeded, so the outcome
  is deterministic and a failure is reproducible. The 3σ tolerance is left as it was; widening it
  would change what the tests assert. Its `using Statistics` also had to be declared in
  `[targets] test`.

- **`Cthulhu` is removed from `[deps]`.** Dead — the only occurrences anywhere are two
  commented-out `using Cthulhu` lines in scripts — with no `[compat]` entry, and it does not
  precompile against Julia 1.14-DEV (`UndefVarError: ConstPropResult not defined in Compiler`),
  which was the sole cause of the advisory `nightly` CI job failing. Anyone using it
  interactively should add it to their own environment rather than to the package's dependencies.

- **Three Landau driver scripts contained invalid Julia and had never been runnable.**
  `const landau_rhs!(v̇, v, params) = …` is a syntax error — `const` takes an assignment, not a
  function definition — in `scripts/landau_new.jl`, `scripts/landau_newer.jl` and
  `scripts/landau_profile.jl`. Found by JuliaFormatter, which could not reparse its own output.
  The last of the three also called `collisions_rhs!`, which does not exist; it is now
  `collisional_vectorfield!` like the other two.

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
  them is commented out, so the import goes rather than the dependency coming back.

## Open Issues

Carried over from the audit that accompanied the `SimpleSplines` migration. None of these are
regressions; each is either a numerical-methods decision or work the migration deliberately did
not take on.

- **The implemented Landau scheme is not the one the main text derives.** The manuscript builds
  the gradient form with the `G` operator — whose structure *is* the momentum and energy
  conservation proof — and a Gonzalez discrete gradient, which *is* the discrete H-theorem proof.
  What runs is the **appendix** two-step `v̇ = K⁺LJ` with plain implicit midpoint and
  `∇S(midpoint)`, which is not a discrete gradient. Neither structural proof transfers to the code
  as written. `G` is never formed; the gradient form survives only as a commented-out
  `Landau_rhs`. A gap between paper and code, not an error in either.

- **The Landau Picard solver does not iterate to convergence.** `Landau_solver.jl` runs exactly
  five iterations and prints the residual without testing it. `tol`, `ftol`, `β`, `m` and
  `chunksize` are accepted and unused, and `probN` is constructed and never solved — it is left
  in place, with the commented-out `NonlinearSolve` calls it belongs to, rather than deleted. Every
  conservation property in the appendix is a property of the *exactly* solved implicit system, so
  momentum and energy drift at the size of that printed residual.

- **The accepted Landau solution is one Picard step behind its stored derivative.** The loop's
  last action recomputes `v̇` at the midpoint from the newest guess without updating the guess, so
  the stored state and derivative do not correspond and the next step's Hermite extrapolation is
  fed an inconsistent pair.

- **The Landau kernel's coincident-point value is a regularisation, not a limit.** `kernel`
  returns zero at `|u| = 0` so that a product quadrature sharing nodes does not produce `Inf`. The
  error does not vanish under refinement, and because both factors use the same Gauß–Legendre
  nodes it fires on every *diagonal cell pair* rather than on a set of measure zero. The kernel is
  integrable in two dimensions; what is needed is a singularity-aware rule or offset grids.

- **Positivity of `f_s` is detected, not solved.** Every `log f_s` and `1/f_s` now throws rather
  than continuing, which is strictly better than the old `0.5·log(f_s²)` returning `log|f_s|`. But
  a run whose projection undershoots now stops, and the real fix — a positivity-preserving
  projection — is the manuscripts' own open problem.

- **The rank hypothesis behind `K K⁺ = I` is no longer checked anywhere.** The step to
  `eq:particle_ode_final` needs `K` to have full row rank. The old code detected the failure with
  two SVDs per vector-field evaluation, printed, and proceeded regardless; the check is gone from
  the hot path and has not been given a home in a diagnostic script.

- **No test covers any structure-preservation claim for Landau.** `scripts/verify_conservation.jl`
  covers the conservative Lenard-Bernstein operator only, and `test/projections_tests.jl` and
  `test/electric_field_tests.jl` remain commented out of `runtests.jl`.

- **The cumulant-scaling experiment is not implemented.** The appendix fixes `A₀ = 0`, `A₁ = 1` by
  hand; no code does. The cumulant computation is commented out in
  `scripts/lenard_bernstein_conservative.jl` and the experiment survives as a stale `run_name`.

- **`src/projections/potential.jl` is dead.** It dispatches on
  `PoissonSolvers.Potential{<:PeriodicBSplineBasis}`, and that name now resolves to the
  `SimpleSplines` type, so it can never match a `PoissonSolvers` basis. It was already unreachable
  before the migration, and its body calls a `Splines.PeriodicVector` that was never imported. Left
  in place rather than deleted.

- **Seven `src/` files are included by nothing:** `electric_field.jl`, the root `vlasov_poisson.jl`
  (distinct from `models/vlasov_poisson.jl`), `visualisation.jl`, `methods/lbm_solver.jl`,
  `hdf5.jl`, and two whose `include` lines are commented out. `src/hdf5.jl` is untracked in git.

- **Dependencies that are no longer used are still declared.** `NaNMath` appears only in a
  commented import; `Plots`, `LaTeXStrings`, `StatsPlots`, `StatsBase`, `SciMLBase` and
  `AdaptiveRejectionSampling` have no occurrences by name in `src/`. Nine non-stdlib
  dependencies still carry no `[compat]` entry. Four driver scripts also `using` `GLMakie`,
  `Printf` and `Profile`, none of which are declared. `ExplicitImports.jl` *has* now been run
  and reports no stale or improper explicit imports, so what remains is `[deps]` hygiene rather
  than dead `import` lines — `Aqua.test_stale_deps` is the check that settles it.

- **`fatou lint` reports 10 warnings, of which one is deliberate and nine are a known false
  positive.** The nine are `unused-import` on `src/VlasovMethods.jl`, where the rule does not
  follow `include` and so flags the module file's load-bearing imports; `ExplicitImports`
  contradicts all nine. The tenth is `probN` below. An earlier version of this changelog and of
  the pull request described `fatou lint` as clean, which was not reproducible.

- **This package cannot be registered until `SimpleSplines` is.** `SimpleSplines` is not in
  General, so it is resolved through a `[sources]` table — which RegistryCI rejects outright — and
  its only version, `1.0.0-DEV`, is not a parseable compat bound. The `[sources]` table also forces
  `julia = "1.11"` rather than the tree's LTS floor of 1.10. Delete the table and add
  `SimpleSplines = "<version>"` the moment it is registered.
