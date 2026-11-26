```@meta
CurrentModule = VlasovMethods
```

# Distribution functions

The distribution function ``f(x, v)`` of a species is the state of every model in the package.
VlasovMethods stores it in one of five forms. Each form is a subtype of
`DistributionFunction{DT, XD, VD}`, with the element type `DT`, the number `XD` of position
dimensions and the number `VD` of velocity dimensions.

| type | stores | used by |
|:--|:--|:--|
| [`ParticleDistribution`](@ref) | positions, velocities and weights of ``N`` particles | every particle model |
| [`SplineDistribution`](@ref) | the coefficients of a spline in ``v`` | the collision operators, the entropy |
| [`GridDistribution`](@ref) | the values at the nodes of a 1D1V phase-space grid | the grid Vlasov–Poisson model, the reduced tensors |
| [`TensorSplineDistribution`](@ref) | the coefficients of a tensor product spline in ``x`` and ``v`` | the spline Vlasov–Poisson model |
| [`LowRankDistribution`](@ref) | the factors ``X``, ``S`` and ``V`` of a 1D1V phase-space function of low rank | the low-rank methods |

The collision operators use two of them at once: the particles carry the state, and a spline
holds the projection of the particles, from which the operator computes ``f_s'/f_s`` or
``\log |f_s|``. [Projections](projections.md) describes the conversions between the forms.

A spline or a grid distribution is a callable: `dist(z...)` evaluates it at the phase-space point
`z`, the `XD` position components first and the `VD` velocity components after them.

## Particles

A [`ParticleDistribution`](@ref) is the Klimontovich measure

```math
f_p(x, v) = \sum_{\alpha = 1}^{N} w_\alpha \, \delta(x - x_\alpha) \, \delta(v - v_\alpha)
```

of ``N`` particles. `ParticleDistribution(xdim, vdim, npart)` allocates the particles of `xdim`
position and `vdim` velocity components, all zero:

```@example particles
using VlasovMethods

dist = ParticleDistribution(1, 1, 5)
dist.particles.x .= [0.1 0.3 0.5 0.7 0.9]
dist.particles.v .= [-1.0 -0.5 0.0 0.5 1.0]
dist.particles.w .= 1 / 5
dist.particles.z
```

The particles are a `ParticleList` of
[ParticleMethods](https://github.com/JuliaPlasma/ParticleMethods.jl). It stores one matrix with
one column per particle, and gives four views of it:

  - `particles.x`, the `xdim × N` positions,
  - `particles.v`, the `vdim × N` velocities,
  - `particles.z`, the `(xdim + vdim) × N` phase-space coordinates, positions first,
  - `particles.w`, the `1 × N` weights.

A write to one view changes the others. The splitting methods advance `particles.z`, and the
collision operators advance `particles.v`.

## Splines in velocity

A [`SplineDistribution`](@ref) holds a B-spline in the velocity,

```math
f_s(v) = \sum_I \hat{f}_I \, \varphi_I(v) ,
```

with one basis function per index ``I`` in one velocity dimension and a tensor-product basis in
more. It has no position dependence. The constructor takes the dimensions, the number of
breakpoints, the spline order, the velocity domain, the width of an optional oversized end cell,
and the boundary condition:

```@example splines
using VlasovMethods

sdist = SplineDistribution(1, 1, 41, 4, (-10.0, 10.0))
```

The order is ``k = p + 1``, the convention of
[SimpleSplines](https://github.com/JuliaGNI/SimpleSplines.jl): order `4` is the cubic basis that
the collision manuscripts use. A clamped basis of order ``k`` on ``n`` breakpoints has
``n + k - 2`` basis functions, so the basis above has 43 of them:

```@example splines
nbasis(sdist.basis)
```

The fields of a spline distribution are the `spline`, the `basis`, the `quadrature` and the
`coefficients`. The spline shares its coefficient array, so a projection that writes the
coefficients changes the function that `sdist.spline` evaluates, without a rebuild:

```@example splines
project_Maxwellian(sdist)
(sdist.spline(0.0), 1 / sqrt(2π))
```

In two or three velocity dimensions, the basis is a `TensorProductBasis`, and the coefficients are
an array of `VD` dimensions. The mass matrix is the Kronecker product of the one-dimensional mass
matrices. It is never formed: a solve is one one-dimensional solve per axis.

```@example splines
sdist2 = SplineDistribution(1, 2, 9, 3, (-4.0, 4.0))
size(sdist2.coefficients)
```

### Particles that leave the domain

The velocity domain of a spline is bounded, and the velocities of the particles are not. A
particle outside the domain contributes nothing to the projection, so its weight would be lost
from ``f_s``. The projection therefore throws a `DomainError` for such a particle. Two remedies
exist: a larger domain, or oversized end cells. The argument `length_big_cell > 0` appends one
cell of that width at each end of the domain:

```@example splines
wide = SplineDistribution(1, 1, 41, 4, (-10.0, 10.0), 5.0)
domain(wide.basis)
```

### The boundary condition

The boundary condition decides which polynomials the basis reproduces. The moments
``\int v^k f_s \, dv`` of the projection of a particle distribution agree with the particle
moments ``\sum_\alpha w_\alpha v_\alpha^k`` only where ``v^k`` is in the span of the basis.
`polynomial_reproduction` reports the largest such degree:

```@example splines
[bc => polynomial_reproduction(SplineDistribution(1, 1, 41, 4, (-10.0, 10.0), 0, bc).basis)
 for bc in (Free(), Periodic(), Dirichlet())]
```

`Free()`, the default, is the unconstrained clamped basis. It is the only one of the three that
reproduces ``1``, ``v`` and ``v^2``, so it is the only one on which the mass, the momentum and the
energy of ``f_s`` are those of the particles. [`check_conservation_basis`](@ref) throws where a
basis does not reproduce a given degree.

## Grids in phase space

A [`GridDistribution`](@ref) holds the values ``f_{ij} = f(x_i, v_j)`` at the nodes of an
``n_x \times n_v`` grid. The ``x`` axis is periodic, so its right end is not stored. The ``v``
axis is bounded, and both of its ends are nodes:

```@example grid
using VlasovMethods

grid = GridDistribution(8, 65, (0.0, 2π), (-6.0, 6.0))
grid.values .= (1 .+ 0.1 .* cos.(grid.x)) .* exp.(-grid.v' .^ 2 ./ 2) ./ sqrt(2π)
(step(grid.x), step(grid.v), grid(π / 4, 0.0))
```

`grid(x, v)` interpolates bilinearly. It is periodic in ``x`` and zero outside the velocity range.
[`velocity_moments`](@ref) integrates over ``v`` at each ``x``-node and returns the density, the
momentum and the kinetic energy density:

```@example grid
velocity_moments(grid).density
```

## The entropy

The collision operators dissipate the entropy of the spline projection,

```math
S_h = \int_\Omega f_s \log |f_s| \, dv .
```

A [`CollisionEntropy`](@ref) holds the spline distribution of a collision operator, and calling it
evaluates ``S_h`` on the Gauß–Legendre nodes of the basis:

```@example splines
entropy = CollisionEntropy(sdist)
entropy()
```

The projected unit Maxwellian has the entropy ``-\tfrac{1}{2}(1 + \log 2\pi) \approx -1.4189``.
The collision operators take a `CollisionEntropy` as their second argument, and use its spline
distribution for every projection.
