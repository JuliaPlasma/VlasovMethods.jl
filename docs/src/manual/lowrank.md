```@meta
CurrentModule = VlasovMethods
```

# Dynamical low rank

A dynamical low-rank method evolves a 1D1V distribution function in the factored form

```math
f(x, v) \approx \sum_{a, b = 1}^{r} X_a(x) \, S_{ab} \, V_b(v) ,
```

with ``r`` x factors ``X_a``, ``r`` v factors ``V_b`` and an ``r \times r`` coefficient matrix
``S``. The state costs ``(n_x + n_v) r + r^2`` numbers in place of ``n_x n_v``.

## The factor space

The factors live in a factor space, a subtype of
[`VlasovMethods.AbstractFactorSpace`](@ref). The space gives the lengths of the factors,
`nx, nv = size(space)`, and the weighted inner products of the x and the v factors with their
orthonormalisations: `VlasovMethods.inner_x`, `VlasovMethods.inner_v`,
`VlasovMethods.orthonormalize_x` and `VlasovMethods.orthonormalize_v`. These four functions and
the abstract type are documented, but not exported.

[`GridFactorSpace`](@ref) is the factor space of the nodes of a
[`GridDistribution`](@ref): a factor is a vector of its values at the nodes, and the inner
products are the rectangle rules of the grid,
``\langle a, b \rangle_x = h_x \sum_i a_i b_i`` and ``\langle a, b \rangle_v = h_v \sum_j a_j b_j``;
the rule in ``v`` is the rule of `velocity_moments` of a grid distribution. The
orthonormalisations return `Q` and `R` of `A = Q R`, with ``h_x Q^\mathsf{T} Q = I`` or
``h_v Q^\mathsf{T} Q = I``.

## The low-rank distribution

[`LowRankDistribution`](@ref) holds the factor space, the factors `X`, `S` and `V`, and the
number `nfixed` of leading columns of `V` that span fixed velocity modes. It factorises the node
values of a grid distribution at a rank `r`. With `fixed = [1 v v²]`, the first three columns of
`V` span 1, ``v`` and ``v^2``, and the rest of `V` comes from the weighted truncated SVD of the
part of ``f`` outside them. The factorisation is then the weighted orthogonal projection of ``f``
onto the span of `V`, which keeps the density, the momentum and the energy of ``f``:

```@example lowrank
using VlasovMethods

nx, nv = 32, 64
dist = GridDistribution(zeros(nx, nv), (0.0, 4π), (-10.0, 10.0))
dist.values .= [exp(-(v - 0.5sin(x / 2))^2 / 2) / sqrt(2π) for x in dist.x, v in dist.v]

v = collect(dist.v)
d = LowRankDistribution(dist, 5; fixed = hcat(one.(v), v, v .^ 2))
m, md = velocity_moments(dist), velocity_moments(d)
(rank = size(d.S, 1), nfixed = d.nfixed,
 moment_errors = map((a, b) -> maximum(abs, a - b) / maximum(abs, a), m, md))
```

Without fixed modes the factorisation is the weighted truncated SVD of rank `r`, and the moments
are not kept. Where the part of ``f`` outside the fixed modes has a numerical rank below
`r - nfixed`, the remaining columns of `V` come from a Householder QR, with no random numbers, so
a factorisation is the same for every state of the global random number generator.

`velocity_moments(d)` and `evaluate(d, x, v)` work on the factors and never form the
``n_x \times n_v`` array; `GridDistribution(d)` forms it, for tests and diagnostics:

```@example lowrank
f = GridDistribution(d)
maximum(abs, f.values - dist.values) / maximum(abs, dist.values)
```

A step of a method replaces the factors of the mutable `LowRankDistribution`, also by factors
of another rank. `LowRankDistribution(space, X, S, V; nfixed)` wraps given factors on any factor
space, and checks their sizes but not their orthonormality.

## The factored fields

A dynamical low-rank step evolves the factors through the projections of the full-order
right-hand side ``F(f) = [h, f]`` of the grid model, with ``h = \phi \otimes 1 + 1 \otimes v^2/2``:
the K-field ``\dot{K} = h_v F(K V^\mathsf{T}) V``, the L-field
``\dot{L} = h_x F(X L^\mathsf{T})^\mathsf{T} X`` and the S-field
``\dot{S} = h_x h_v X^\mathsf{T} F(X S V^\mathsf{T}) V``. The workspace
[`VlasovMethods.LowRankFields1D1V`](@ref) holds the Arakawa bracket, the Green's matrix and the
kinetic energy ``v^2/2`` of the grid, and the functions `VlasovMethods._kfield!`,
`VlasovMethods._lfield!`, `VlasovMethods._sfield!` and `VlasovMethods._potential!` evaluate the
K-, L- and S-fields and the potential on the factors, without forming ``f``. The Hamiltonian
``h`` has rank 2, so each term of Arakawa's stencil projects onto the other factor through
``r \times r`` matrices of shifted columns, and a field costs ``O((n_x + n_v) r^2)`` beyond the
product with the Green's matrix. Each field also takes a given, frozen potential ``\phi``. The
workspace and the functions are documented, but not exported.

## The projector splitting

[`ProjectorSplitting`](@ref) is the projector-splitting integrator of Lubich and Oseledets, a
reference method: it splits the tangent projector of the low-rank manifold into a K-substep with
`V` fixed, a backward S-substep, and an L-substep with `X` fixed, and orthonormalises `K` and `L`
by a QR in the weighted products of the factor space. `ProjectorSplitting(:lie)` is the
first-order composition; `ProjectorSplitting(:strang)` is the second-order composition of
Einkemmer and Lubich, which runs its L-substeps at a frozen potential and predicts the midpoint
potential by a half step. Each substep is one RK4 step. [`lowrank_step!`](@ref) advances a
`LowRankDistribution` by one step in place, with the fields of the grid as the right-hand side,
and sets `nfixed = 0`, because the splitting does not keep fixed velocity modes:

```@example lowrank
using PoissonSolvers: FFTWBasis, Potential

model = VlasovPoisson(dist, Potential(FFTWBasis((0.0, 4π), nx)))
w = VlasovMethods.LowRankFields1D1V(d.space, model.work)
Y = LowRankDistribution(dist, 5)
method = ProjectorSplitting(:strang)
for n in 0:9
    lowrank_step!(Y, method, w, 0.05n, 0.05)
end
(rank = size(Y.S, 1), nfixed = Y.nfixed,
 orthonormality = maximum(abs, VlasovMethods.inner_x(Y.space, Y.X, Y.X) - one(Y.S)))
```

The step takes any field with methods of the internal functions `_lowrank_kfield!`,
`_lowrank_sfield!` and `_lowrank_lfield!`, and for `:strang` also `_lowrank_potential`, which the
docstring of [`lowrank_step!`](@ref) lists.
