```@meta
CurrentModule = VlasovMethods
```

# Grid-based reduced tensors

The grid Vlasov–Poisson model of [Vlasov–Poisson](vlasov_poisson.md#Grids) is
``\dot{f} = -[f, h]`` with the Arakawa bracket on an ``n_x \times n_v`` grid. The bracket is
bilinear, so it is a tensor of three indices over the grid nodes:

```math
[f, h]_I = \sum_{J, K} T_{IJK} \, f_J \, h_K .
```

The Hamiltonian splits into the kinetic part ``v_j^2 / 2``, which does not depend on ``f``, and
the potential ``\phi_i``, which depends on ``x`` only. A Galerkin reduction approximates
``f \approx \sum_j P_j a_j`` with a reduced basis ``P``, and projects the equation onto the test
basis. The tensors on this page are the reduced forms of ``T`` that this reduction needs. They
are lazy: an entry is summed when it is read, over the stencil of the bracket only, and nothing of
the full size ``(n_x n_v)^3`` is ever formed.

| type | entry | contracts the third index with |
|:--|:--|:--|
| [`PotentialReducedTensor`](@ref)`(T, Pi, Pj, Pk)` | ``\sum_{I,J,K} T_{IJK} P^i_{Ii} P^j_{Jj} P^k_{Kk}`` | a basis `Pk` of grid functions |
| [`VelocityReducedMatrix`](@ref)`(T, Pi, Pj, v²)` | ``\sum_{I,J,K} T_{IJK} P^i_{Ii} P^j_{Jj} \, v_{j(K)}^2 / 2`` | the kinetic energy |
| [`FullyReducedTensor`](@ref)`(T, Pi, Pj, Px)` | ``\sum_{I,J,K} T_{IJK} P^i_{Ii} P^j_{Jj} P^x_{i(K) k}`` | a basis `Px` of functions of ``x`` |

The grid nodes are counted in the order of `vec(dist.values)` of a
[`GridDistribution`](@ref), the ``x`` index fastest. [`bracket_tensor`](@ref) builds ``T`` from
the `Arakawa` bracket of GeometricBrackets:

```@example tensors
using GeometricBrackets: Arakawa
using VlasovMethods

nx, nv = 4, 5
T = bracket_tensor(Arakawa(nx, nv, 0.5, 0.5))
size(T)
```

Two reduced bases of three and two columns give the reduced tensors:

```@example tensors
P = [cos(π * (i + 2j) / (nx * nv)) for i in 1:(nx * nv), j in 1:3]
Px = [cos(2π * (i - 1) * k / nx) for i in 1:nx, k in 0:1]
v = range(-1.0, 1.0; length = nv)

(potential = size(PotentialReducedTensor(T, P, P, P)),
 velocity = size(VelocityReducedMatrix(T, P, P, collect(v) .^ 2)),
 fully = size(FullyReducedTensor(T, P, P, Px)))
```

The projections need one row per grid node, or one per ``x``-node for `Px`, and the element type
of the tensor. The constructors check both, and an ``x`` or ``v`` direction of fewer than three
nodes, which the stencil of the bracket would visit twice.

## Moments on the grid

[`velocity_moments`](@ref) computes the density, the momentum and the kinetic energy density of a
`GridDistribution` at each ``x``-node, by the rectangle rule in ``v``. The stencils behind it,
such as [`_apply_∫dv!`](@ref), act on the state vector `vec(dist.values)` and write into a
preallocated vector.

## Collision tensors

The grid discretisation of a collision operator ``C[f]`` that is at most quadratic in ``f`` is a
tensor of the same kind. The internal types `QuadraticCollisions`, `CollisionTensor` and
`ReducedCollisionTensor` hold its stencil, the lazy tensor on the grid and its reduction. The
internal functions [`_get_MC̃_quadratic`](@ref) and [`_get_MC̃_cubic`](@ref) assemble the
Galerkin forms of ``\rho^2 C[f]``. The cubic tensor conserves the mass, the momentum and the
energy; the quadratic one conserves the mass and the energy only. Both are valid while ``f``
vanishes at the ends of the velocity grid, because the velocity stencil wraps around. None of
these is exported, and no model of the package uses them; the tests check their conservation.
