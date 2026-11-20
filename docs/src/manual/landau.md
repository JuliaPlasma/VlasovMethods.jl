```@meta
CurrentModule = VlasovMethods
```

# Landau collisions

The Landau operator describes the collisions of charged particles. For one species in two velocity
dimensions it reads

```math
\frac{\partial f}{\partial t} = \nu \, \nabla_v \cdot \int U(v - v') \,
  \bigl( f(v') \nabla_v f(v) - f(v) \nabla_{v'} f(v') \bigr) \, dv' ,
```

with the kernel

```math
U(u) = |u|^{\gamma + 2} \left( \mathbb{I} - \frac{u \otimes u}{|u|^2} \right) .
```

The exponent ``\gamma = -3`` is the Coulomb case, and ``\gamma = 0`` gives the Maxwell molecules.
[`kernel`](@ref) evaluates ``U``. The operator conserves the mass, the momentum and the energy,
and it decreases the entropy ``\int f \log f \, dv``.

The MP-L manuscript writes the operator as a metric bracket and discretises it on particles. The
package implements both of its methods. Both project the particles onto a spline ``f_s`` of two
velocity dimensions, and both use the gradient of the discrete entropy
``S_h = \int f_s \log |f_s| \, dv``.

| model | method | vector field | time integrator |
|:--|:--|:--|:--|
| [`PairwiseLandau`](@ref) | main text | ``\nu \, \mathbb{G}(v) \nabla S_h(v)``, summed over particle pairs | [`discrete_gradient_step!`](@ref) |
| [`Landau`](@ref) | appendix | least-squares solve with the spline matrix ``\mathbb{L}`` | [`Picard_iterate_Landau_nls!`](@ref) |

## The particle-pair method

[`PairwiseLandau`](@ref) is the gradient system

```math
\dot{v}_\gamma = - \nu \sum_{\alpha \ne \gamma} w_\alpha \, U(v_\gamma - v_\alpha)
  \bigl( \psi_\gamma - \psi_\alpha \bigr) , \qquad
\psi_\alpha = \sum_k \mathbb{L}_k \nabla \varphi_k(v_\alpha) ,
```

where ``\mathbb{L}_k`` is the ``L^2`` projection of ``1 + \log |f_s|`` onto the spline basis
([`compute_J!`](@ref)). The sum over the particle pairs makes the matrix ``\mathbb{G}`` of
[`apply_G!`](@ref) symmetric and negative semi-definite, and its null space holds the momentum and
the energy. The cost of one evaluation is ``O(N^2)`` kernel evaluations.

[`discrete_gradient_step!`](@ref) takes one step of the discrete gradient method

```math
\frac{v^{n+1} - v^n}{\Delta t} = \nu \, \mathbb{G}\!\left( \frac{v^n + v^{n+1}}{2} \right)
  \bar{\nabla} S_h(v^n, v^{n+1}) ,
```

solved by a fixed-point iteration to the absolute residual `abstol`. The step writes the new
velocities into the particles, and returns the number of iterations and the final residual. It
conserves the momentum to rounding, the energy to the tolerance of the iteration, and the discrete
gradient gives the discrete H-theorem. Two discrete gradients are implemented: the
[`GonzalezGradient`](@ref), the default and the choice of the manuscript, and the
[`AVFGradient`](@ref), the average of the gradient along the segment.

The example below relaxes the BKW profile of 200 particles over 20 steps, on quadratic splines of
8 cells per axis:

```@example pairwise
using VlasovMethods

dist = initialize!(ParticleDistribution(1, 2, 200), BKWDistribution(0.2))

# the velocity domain covers the sample with a margin
halfwidth = 1.5 * maximum(abs, dist.particles.v)
sdist = SplineDistribution(1, 2, 9, 3, (-halfwidth, halfwidth))

model = PairwiseLandau(dist, CollisionEntropy(sdist); γ = -3, ν = 1.0)

w = vec(dist.particles.w)
entropy = CollisionEntropy(sdist)
function diagnostics(v)
    projection(v, dist, sdist)
    return (momentum = v * w, energy = sum(w .* vec(sum(abs2, v; dims = 1))) / 2,
        entropy = entropy())
end

before = diagnostics(copy(dist.particles.v))
iterations = [discrete_gradient_step!(model, 1e-4; abstol = 1e-10).iterations for k in 1:20]
after = diagnostics(copy(dist.particles.v))

(iterations = iterations,
 momentum_drift = maximum(abs, after.momentum .- before.momentum),
 energy_drift = abs(after.energy - before.energy) / before.energy,
 entropy = (before.entropy, after.entropy))
```

The momentum and the energy drift at the level of rounding, and the entropy decreases.

### Errors

The step refuses what the scheme cannot handle. Two particles of equal velocities throw a
`DomainError`, because the kernel of a particle pair at ``u = 0`` is not defined; for
``\gamma < -2`` so do two velocities whose distance is too small for the kernel to be finite. A particle that leaves the
velocity domain of the spline throws the `DomainError` of the projection. A fixed-point iteration
that does not reach `abstol` in `max_iters` iterations throws an `ErrorException`. In each case,
the particles keep the velocities of the last accepted step. The docstring of
[`discrete_gradient_step!`](@ref) lists the residual floor that bounds a useful `abstol`.

## The spline method

[`Landau`](@ref) is the method of the appendix of MP-L. It assembles the spline matrix

```math
\mathbb{L}_{ij} = \frac{1}{2} \int\!\!\int
  \bigl( \nabla \varphi_i(v) - \nabla \varphi_i(v') \bigr) \cdot
  |f_s(v)| \, U(v - v') \, |f_s(v')| \,
  \bigl( \nabla \varphi_j(v) - \nabla \varphi_j(v') \bigr) \, dv \, dv'
```

of [`compute_L!`](@ref) on the tensor-product quadrature grid, at a cost of ``O(Q^2)`` kernel
evaluations for ``Q`` quadrature points, and the matrices
``\mathbb{K}_{k,\alpha} = w_\alpha \nabla \varphi_k(v_\alpha)`` of [`compute_K!`](@ref). The
velocities follow from the least-squares problem ``\mathbb{K} \dot{v} = \mathbb{L} \, \mathbb{L}_k``.
The appendix needs ``\mathbb{K} \mathbb{K}^{+} = \mathbb{I}``, which requires at least as many
particle velocity components as basis functions less one, ``dN \ge M - 1``, and a particle in the
support of every basis function.

[`Picard_iterate_Landau_nls!`](@ref) takes one implicit-midpoint step, solved by a Picard
iteration. The implicit midpoint rule is the appendix's choice. It is not a discrete gradient, so
no discrete H-theorem holds for the spline method: the entropy of a step need not decrease. The
discrete H-theorem holds for the particle-pair method with its discrete gradient step. The caller
keeps the velocities of the two previous steps and the two previous vector fields:

```julia
landau = Landau(dist, CollisionEntropy(sdist); γ = -3, ν = 1.0)

v = copy(dist.particles.v)
v_prev_2 = copy(v)
rhs = zeros(size(v)..., 2)
for k in 1:nsteps
    VlasovMethods.Picard_iterate_Landau_nls!(landau, Δt, k, k * Δt, v, v_prev_2, rhs;
        abstol = 1e-10)
    v_prev_2 .= v
    v .= dist.particles.v
end
```

!!! warning "Known issue of the spline method"
    The implementation solves the least-squares problem for each velocity component alone, and
    its vector field has the opposite sign of the appendix's ``-\mathbb{K}^{+} \mathbb{L} \,
    \mathbb{L}_k``. Under it, the temperature anisotropy grows at the rate at which it should
    decay, and the energy drifts by ``10^{-8}`` to ``10^{-7}`` in the runs of
    [Particle methods for the Landau operator](@ref). K96 of `KNOWN_ISSUES.md` records the measurements.
    The particle-pair method has neither defect.

## The reference cases

The three two-dimensional initial conditions of [Initial conditions](initial_conditions.md) are
the reference cases of the operator: the [`AnisotropicBiMaxwellian`](@ref), whose temperatures
isotropise at the closed rate ``8`` for ``\gamma = 0``, the [`DoubleMaxwellian2D`](@ref), and the
[`BKWDistribution`](@ref), an exact solution for ``\gamma = 0``. The script
`scripts/reference_forms.jl` derives their closed forms, and
[Particle methods for the Landau operator](@ref) runs them under both methods.
