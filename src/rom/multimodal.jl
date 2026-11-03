# The Vlasov frontends of the multi-modal symplectic autoencoder are a package extension. This
# file declares the functions and documents them; `ext/MultiModalAutoEncoderExt.jl` holds their
# methods and the types they return, so that loading the package does not load
# `GeometricMachineLearning`, and each function has no methods until the extension loads.

@doc raw"""
    build_spline_parity_factors(basis_x, basis_v, mx, mv) -> BilinearStiefelLayer

Return the `BilinearStiefelLayer` of `GeometricMachineLearning` that encodes the coefficients
of a tensor-product spline on `basis_x ⊗ basis_v` with ``m_x`` spatial and ``m_v`` velocity
modes, with `metric_factors = (R_x, R_v)`, ``M_x = R_x^\mathsf{T} R_x`` and
``M_v = R_v^\mathsf{T} R_v`` for the mass matrices of the two bases.

Defined in a package extension: load `GeometricMachineLearning` to make it available.

``R_x`` is the upper Cholesky factor of ``M_x``. ``R_v`` is a square-root factor and not
triangular: ``R_v = S^\mathsf{T} (R^+ \oplus R^-) S``, with ``S`` the orthogonal parity
split of the layer (index ``j`` with its mirror ``N_v + 1 - j``, scaled by ``1/\sqrt 2``) and
``R^\pm`` the upper Cholesky factors of the even and odd blocks of ``S M_v S^\mathsf{T}``. A
triangular factor of the whole ``M_v`` does not commute with the reflection ``J`` of the
velocity coefficients, so an even latent state would decode to a function that is not even;
``R_v`` commutes with ``J``, and the layer applies ``R^+`` to the even half and ``R^-`` to the
odd half.

With ``\Phi_x`` and ``\Phi_v`` the Stiefel bases of the layer, the spatial modes
``R_x^{-1} \Phi_x`` and the even and odd velocity modes
``R_v^{-1} S_\pm^\mathsf{T} \Phi_v`` are orthonormal in ``M_x`` and in ``M_v``: the latent
coordinates are ``L^2``-orthonormal projections of the spline. They are not canonical
coordinates of the Vlasov bracket: linearised at an equilibrium ``f_0``, the bracket of the even
and the odd coordinates is ``-D_x \otimes K_v``, not the identity. [`dual_parity_modes`](@ref)
gives odd modes for which it is the identity.

# Throws

- `ArgumentError` when the velocity domain is not symmetric, ``v_{\min} = -v_{\max}`` to
  ``10\,\mathrm{eps}(T)`` relative to ``v_{\max}``, or the velocity knot vector is not
  symmetric under ``v \mapsto -v``, and when `basis_v` has an odd number of functions. A
  clamped basis of degree ``p`` on ``N`` cells has ``N + p``, so ``N + p`` must be even.
"""
function build_spline_parity_factors end

@doc raw"""
    build_grid_parity_factors(grid_x, grid_v, mx, mv) -> BilinearStiefelLayer

Return the `BilinearStiefelLayer` of `GeometricMachineLearning` that encodes the values of a
phase-space function on the uniform grid `grid_x × grid_v` with ``m_x`` spatial and ``m_v``
velocity modes, with `metric_factors = (Diagonal(√w_x), Diagonal(√w_v))`.

Defined in a package extension: load `GeometricMachineLearning` to make it available.

Each node carries the width of one cell, ``w_x = \Delta x`` and ``w_v = \Delta v``, the weights
of the periodic trapezoidal rule and of the midpoint rule. The layer therefore applies
``\sqrt{\Delta x\, \Delta v}`` to every value, and its coordinates approximate the ``L^2``
inner products of the function with the modes.

# Throws

- `ArgumentError` when `grid_v` is not symmetric under ``v \mapsto -v`` to
  ``10\,\mathrm{eps}(T)`` relative to its largest node, and when it has an odd number of nodes.
"""
function build_grid_parity_factors end

@doc raw"""
    tie_continuum_frontends!((spline_layer, ps_spline), (grid_layer, ps_grid),
        (dlra_layer, ps_dlra), E_x, E_v) -> (ps_grid, ps_dlra)

Set the Stiefel bases of the grid and the low-rank frontends from those of the spline
frontend, so that a phase-space function given as a spline, as its grid values or as a
low-rank factorisation of them maps to the same latent state up to the quadrature error of
the grid.

Defined in a package extension: load `GeometricMachineLearning` to make it available.

`spline_layer` is a layer of [`build_spline_parity_factors`](@ref), `grid_layer` one of
[`build_grid_parity_factors`](@ref), and `dlra_layer` a `LowRankStiefelLayer` with the metric
factors of `grid_layer`; `ps_*` are their parameters, and the function overwrites `ps_grid`
and `ps_dlra` with the same bases. `E_x` and `E_v` are the nodal evaluation matrices,
`E_x[i, a]` the value of the spatial basis function `a` at grid node `i`, and `E_v` the same
for the velocity basis.

A spline mode with orthonormal coordinates ``a`` is the function with coefficients
``R^{-1} a``, and its grid coordinates are its nodal values times ``\sqrt{w}``. The grid bases
are the polar factors of these images, ``\Phi_x^{\mathrm{grid}} = \mathrm{polar}(W_x^{1/2}
E_x R_x^{-1} \Phi_x)`` and the even half of the same map in velocity; the polar factor keeps
them orthonormal and differs from the image by the quadrature error of the Gram matrix.

The layers share one velocity basis ``\Phi_v`` between the even and the odd coordinates.
The velocity bases are tied through the even modes: the odd spline modes differ from the odd
extension of the even ones near ``v = 0``, because ``R^+ \ne R^-``, so the odd coordinates of
the three frontends agree only to that difference, which falls more slowly under refinement
than the quadrature error of the even ones.

# Throws

- `ArgumentError` when `dlra_layer` does not carry the metric factors of `grid_layer`.
- `DimensionMismatch` when `E_x` or `E_v` does not map the spline basis to the grid.
"""
function tie_continuum_frontends! end

@doc raw"""
    dual_parity_modes(x_modes, even_v, odd_v, df0; nq = 24) -> Matrix

Return the coefficients of the odd modes ``h = e B^{-1}`` dual to the even product modes
``g_a = \varphi_i \psi^+_j`` for the Vlasov bracket linearised at the homogeneous
equilibrium ``f_0(v)``: with ``Q_a = \langle g_a, \delta f\rangle`` and
``P_b = \langle h_b, \delta f\rangle``, ``\{Q_a, P_b\} = \delta_{ab}`` and
``\{P_b, P_{b'}\} = 0``.

Defined in a package extension: load `GeometricMachineLearning` to make it available.

`x_modes` is a vector of `Spline`s on one periodic basis, `even_v` and `odd_v` vectors of
`Spline`s on one velocity basis, even and odd in ``v``, and `df0` the derivative
``f_0'(v)``. The odd product modes are ``e_c = \varphi_k \psi^-_l``, and
``B = -D_x \otimes K_v`` with ``(D_x)_{ik} = \int \varphi_i' \varphi_k \,\mathrm{d}x`` and
``(K_v)_{jl} = \int f_0'(v)\, \psi^+_j \psi^-_l \,\mathrm{d}v``. The product index runs over
``x`` fastest, the order of `vec(Q_red)` in `BilinearStiefelLayer`. Column ``b`` of the result is ``h_b``,
the coefficients on the tensor basis, again ``x`` fastest. ``K_v`` is assembled with `nq`
Gauß–Legendre points per cell, enough for a Maxwellian on cells of unit width to round-off.

The dual modes are not ``L^2``-orthonormal, so a frontend that encodes ``P`` with them is
not an isometry. Neither a monotone nor a positive ``f_0`` is needed: the construction
holds for the two-stream equilibrium ``f_0 = v^2 M(v)`` as for the Maxwellian ``M``.

# Throws

- `ArgumentError` when the splines of `x_modes`, of `even_v` or of `odd_v` do not share one
  basis, when `even_v` and `odd_v` differ in length, when an `x` mode is constant,
  when ``D_x`` is singular because a ``\cos mx`` lacks its ``\sin mx``, and when ``K_v`` is
  numerically singular, ``\mathrm{cond}(K_v) > 1/(100\,\mathrm{eps}(T))``. Each message
  names its condition.
"""
function dual_parity_modes end

@doc raw"""
    klimontovich_deposition(basis_x, basis_v, weights) -> KlimontovichDepositionLayer

Return the deposition of weighted particles onto the tensor-product spline space
`basis_x ⊗ basis_v`, a periodic basis in ``x`` and a clamped one in ``v``. Called with the
particle positions `x` and velocities `v`, the layer returns the coefficient matrix

```math
\tilde c = M_x^{-1}\, c\, M_v^{-1}, \qquad
c_{ab} = \sum_{i} w_i\, \Lambda^x_a(x_i)\, \Lambda^v_b(v_i),
```

the coefficients of the ``L^2`` projection ``\Pi_h \mathbf{J}`` of the Klimontovich density
``\mathbf{J} = \sum_i w_i \delta(x - x_i) \delta(v - v_i)``, not the moments ``c``. A
particle outside the velocity domain adds nothing.

Defined in a package extension: load `GeometricMachineLearning` to make it available.

Because ``\Pi_h`` is self-adjoint, ``\langle \Pi_h \mathbf{J}, g\rangle =
\langle \mathbf{J}, \Pi_h g \rangle`` for every ``g``. The Galerkin-truncated Lie–Poisson
bracket on ``\tilde c`` matches the particle bracket of two smooth observables up to the
truncation error of the spline space, and it does not satisfy the Jacobi identity at a fixed
mesh.

# Throws

- `DimensionMismatch` when `x`, `v` and `weights` differ in length.
"""
function klimontovich_deposition end

@doc raw"""
    klimontovich_bridge_loss(deposition, spline_layer; ε = eps(T)) -> KlimontovichBridgeLoss

Return the loss that ties a particle frontend to the spline frontend `spline_layer` of
[`build_spline_parity_factors`](@ref) through the deposition `deposition` of
[`klimontovich_deposition`](@ref). Called with the particle latent state `z` (a `NamedTuple`
`(q, p)`), the particle positions `x`, velocities `v` and the parameters `ps` of
`spline_layer`, it returns

```math
\frac{\lVert z - \mathcal{E}_{\mathrm{spline}}(\tilde c) \rVert_2^2}
     {\lVert \mathcal{E}_{\mathrm{spline}}(\tilde c) \rVert_2^2 + \varepsilon},
```

with ``\tilde c`` the deposited coefficients, not the moments.

Defined in a package extension: load `GeometricMachineLearning` to make it available.
"""
function klimontovich_bridge_loss end

@doc raw"""
    kruskal_oberman_free_energy(df0, grid_v, wavenumbers) -> KruskalObermanFreeEnergy

Return the Kruskal–Oberman free energy of the Vlasov–Poisson equations linearised at a
monotone equilibrium ``f_0(v)``, on the spatial Fourier modes ``\cos kx``, ``\sin kx`` of
`wavenumbers` and the symmetric velocity grid `grid_v`.

Defined in a package extension: load `GeometricMachineLearning` to make it available.

A state is the matrix `C` of the values of ``\delta f`` at the nodes `grid_v`, one row per
spatial mode, in the order ``\cos k_1 x, \sin k_1 x, \cos k_2 x, \dots``, each mode normalised
in ``L^2(0, 2\pi)``. The velocity space is nodal, with the grid width ``\Delta v`` as the
weight of each node, so the Galerkin matrix of a multiplication is diagonal. Called with `C`,
the free energy returns `(value, gradient, rhs)`: the value

```math
H_{\mathrm{KO}} = \frac12 \sum_q \Delta v\, w_0(v_q) \lVert C_{:, q} \rVert^2
               + \frac12 \sum_i \frac{\rho_i^2}{k_i^2}, \qquad
w_0 = -\frac{v}{f_0'(v)}, \quad \rho_i = \sum_q \Delta v\, C_{iq},
```

its gradient with respect to `C`, and the Hamiltonian vector field

```math
\partial_t \begin{pmatrix} \delta f^+ \\ \delta f^- \end{pmatrix}
= \begin{pmatrix} 0 & \mathcal{A} \\ \mathcal{A} & 0 \end{pmatrix}
  \begin{pmatrix} \delta H / \delta f^+ \\ \delta H / \delta f^- \end{pmatrix},
\qquad \mathcal{A} = f_0'(v)\, \partial_x ,
```

whose Galerkin matrix is ``\mathrm{i}k K_v`` on the mode ``k``, with
``K_v = \mathrm{diag}(\Delta v\, f_0'(v_q))``. The operator is skew, so the free energy is
conserved, and the vector field is the linearised Vlasov–Poisson right-hand side
``-v\, \partial_x \delta f + f_0'(v)\, \partial_x \Phi``, ``-\partial_x^2 \Phi = \rho``.

For a non-monotone ``f_0``, such as the two-stream equilibrium, ``w_0`` is not defined; the
canonical coordinates of that case are those of [`dual_parity_modes`](@ref).

# Throws

- `ArgumentError` when a wavenumber is zero (the constant mode has no conjugate), when
  `grid_v` is not symmetric or has an odd number of nodes, and when ``f_0'(v) \ge 0`` at a
  node ``v > 0``, which names [`dual_parity_modes`](@ref).
- `DimensionMismatch` when `C` does not have one row per mode and one column per node.
"""
function kruskal_oberman_free_energy end
