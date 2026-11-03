```@meta
CurrentModule = VlasovMethods
```

# Frontends of a multi-modal symplectic autoencoder

A symplectic autoencoder for the Vlasov equation can learn from several representations of one
distribution function: a tensor-product spline, its values on a phase-space grid, a low-rank
factorisation of those values, and a particle ensemble. Each representation needs a frontend that
maps it to the coordinates ``(q, p)`` of a shared core. The functions on this page build the
Vlasov-specific frontends. They are a package extension, and load with
`GeometricMachineLearning`:

```julia
using VlasovMethods, GeometricMachineLearning
```

## Parity layers

The frontends split the velocity axis into its even and its odd part, ``f^+`` and ``f^-``, and
project both onto ``L^2``-orthonormal bases ``\Phi_x \otimes \Phi_v``. The coordinates
``q`` are the projections of ``f^+``, the coordinates ``p`` those of ``f^-``.

  - [`build_spline_parity_factors`](@ref) returns the layer of a spline. Its metric factors
    ``(R_x, R_v)`` make the decoded modes orthonormal in the mass matrices. ``R_v`` is a
    square-root factor that commutes with the velocity reflection, so that an even latent state
    decodes to an even function.
  - [`build_grid_parity_factors`](@ref) returns the layer of a uniform grid, with the factors
    ``\sqrt{\Delta x}`` and ``\sqrt{\Delta v}`` of the cell weights.
  - [`tie_continuum_frontends!`](@ref) sets the bases of the grid layer and of a
    `LowRankStiefelLayer` from those of the spline layer. A function given as a spline, as grid
    values or as a low-rank factorisation then maps to the same latent state, up to the
    quadrature error of the grid.

The velocity basis needs an even number of functions. A clamped spline basis of degree ``p`` on
``N`` cells has ``N + p`` functions, so ``N + p`` must be even.

## Canonical coordinates of the linearised bracket

The ``L^2``-orthonormal coordinates are not canonical for the Vlasov bracket. Linearised at a
homogeneous equilibrium ``f_0(v)``, the bracket of an even and an odd coordinate is
``-D_x \otimes K_v``, with ``D_x`` the derivative in the spatial modes and ``K_v`` the overlap
of the velocity modes weighted with ``f_0'``. [`dual_parity_modes`](@ref) returns the odd modes
``h = e B^{-1}``, ``B = -D_x \otimes K_v``, for which the bracket is the identity. The
construction needs an invertible ``B``: no constant spatial mode, every ``\cos mx`` with its
``\sin mx``, as many even as odd velocity modes, and an invertible ``K_v``. It does not need a
monotone ``f_0``, so it serves the two-stream equilibrium too.

## The particle bridge

[`klimontovich_deposition`](@ref) deposits weighted particles on a tensor-product spline space
and returns the coefficients ``\tilde c = M_x^{-1} c\, M_v^{-1}`` of the ``L^2`` projection of
the Klimontovich density, not the moments ``c``. The spline frontend acts on coefficients, so
[`klimontovich_bridge_loss`](@ref) compares the latent state of a particle frontend with the
spline frontend applied to ``\tilde c``.

The Galerkin-truncated bracket on ``\tilde c`` matches the particle bracket of two smooth
observables up to the truncation error of the spline space, which falls with the mesh width. It
is not a Poisson bracket: it fails the Jacobi identity at every fixed mesh, by a defect that
decreases under refinement.

## The Kruskal–Oberman free energy

For a monotone equilibrium, ``f_0'(v) < 0`` for ``v > 0``,
[`kruskal_oberman_free_energy`](@ref) returns the free energy of the linearised Vlasov–Poisson
equations on spatial Fourier modes and a symmetric velocity grid, with its gradient and its
Hamiltonian vector field. The vector field is that of the skew operator
``\bigl(\begin{smallmatrix} 0 & \mathcal{A} \\ \mathcal{A} & 0 \end{smallmatrix}\bigr)``,
``\mathcal{A} = f_0'(v)\, \partial_x``, and it equals the linearised Vlasov–Poisson right-hand
side, so the free energy is conserved. For a non-monotone ``f_0`` the free energy is not defined,
and the canonical coordinates are those of [`dual_parity_modes`](@ref).
