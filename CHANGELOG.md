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

- **Three Landau driver scripts contained invalid Julia and had never been runnable.**
  `const landau_rhs!(v̇, v, params) = …` is a syntax error — `const` takes an assignment, not a
  function definition — in `scripts/landau_new.jl`, `scripts/landau_newer.jl` and
  `scripts/landau_profile.jl`. Found by JuliaFormatter, which could not reparse its own output.
  The last of the three also called `collisions_rhs!`, which does not exist; it is now
  `collisional_vectorfield!` like the other two.

### Breaking Changes

## Open Issues
