```@meta
CurrentModule = VlasovMethods
```

# Reduced-order models

A reduced-order model replaces the ``2 N`` particle coordinates of a Vlasov–Poisson run by a few
coefficients in a reduced basis. The basis is learned from snapshots of full-order runs at several
values of a parameter, and the reduced model then runs at new values of that parameter. This page
describes the building blocks; [Symplectic model reduction](../scripts/reduced_basis.md) runs the
experiments of the paper with them.

The workflow has four stages.

 1. **Sweep.** March the full-order model at every sample of a parameter space, and save the
    states into [`Snapshots`](@ref).
 2. **Bases.** Compute a symplectic basis of the particle snapshots with [`particle_basis`](@ref),
    and a POD basis of the field snapshots with [`field_basis`](@ref).
 3. **Reduced march.** March the reduced model with a reduced electric field, either the full
    field of the lifted particles, [`ReducedElectricField`](@ref), or its hyper-reduction by the
    discrete empirical interpolation method, [`DEIMElectricField`](@ref).
 4. **Diagnostics.** Compare the full-order and the reduced runs with the functions of
    [Diagnostics](#Diagnostics).

The bases and the DEIM indices come from
[ReducedBasisMethods](https://github.com/JuliaRCM/ReducedBasisMethods.jl), and the parameter
spaces from [ReducedComplexityModeling](https://github.com/JuliaRCM/ReducedComplexityModeling.jl).

## The electric field

The marches on this page take the field as an [`ElectricField`](@ref). Its interface is two
functions: `update!(f, x, w, t)` deposits the positions `x` with the weights `w` and solves for
the potential, and `efield!(f, e, x)` writes the field at the positions `x` into `e`. A field
holds its own buffers, so neither call allocates.

| field | evaluates |
|:--|:--|
| [`PoissonField`](@ref)`(potential)` | the full-order field of a `Potential` of PoissonSolvers |
| [`ScaledField`](@ref)`(field, χ)` | the field divided by ``\chi^2``, for the rescaled time of a parameter sweep |
| [`ReducedElectricField`](@ref)`(field, Ψ, X)` | lift by ``\Psi``, solve, evaluate at every particle, project by ``\Psi^{\mathsf T}`` |
| [`DEIMElectricField`](@ref)`(field, Ψ, Ψₑ, p, X)` | the same, evaluated at the ``k_e`` DEIM points `p` of the field basis ``\Psi_e`` only |

[`energy`](@ref) returns the field energy ``\frac{1}{2} \int (\partial_x \phi)^2 \, dx`` of a
field after its last `update!`.

### The parameter scale

The bump-on-tail sweep varies the wave number ``\kappa``. The particles are drawn once, at the
reference wave number ``\kappa_0``, and a sample at ``\kappa`` marches them on the same domain
with the scale ``\chi = \kappa / \kappa_0``. The scale enters the march as a time-step scale,
``\Delta t \, \chi``, and the field of the march as ``E / \chi^2``. Every march reads ``\chi``
from its parameters, `params.χ`, and a `ScaledField` must be built with the same value; a march
throws an `ArgumentError` where the two differ.

## Sweeps and snapshots

[`Snapshots`](@ref) holds the positions `X`, the velocities `V`, the accelerations `A` and the
potentials `Φ` of every saved state of every sample, and the field energy `W`, the kinetic energy
`K` and the momentum `M`. [`IntegratorParameters`](@ref) holds the time step and the counts of
the sweep. [`integrate_vp`](@ref) marches the full-order model of one sample with the leapfrog
method, and saves its states:

```julia
using PoissonSolvers
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods

κ₀ = 0.3
pspace = ParameterSpace(Parameter(:χ, collect(range(0.1, 0.5; length = 10) ./ κ₀)))

particles = initialize!(ParticleDistribution(1, 1, 10_000), BumpOnTail(κ = κ₀)).particles
poisson = Potential(PeriodicBasisSpline((0.0, 2π / κ₀), 4, 16))   # cubic: the order is 4

dt, nt = 0.1, 250
nh = length(PoissonSolvers.coefficients(poisson))
IP = IntegratorParameters(dt, nt, nt + 1, nh, length(particles), length(pspace))
SS = Snapshots(Float64, 1, length(particles), nh, nt + 1, pspace)

for p in 1:length(pspace)
    params = pspace(p)
    field = ScaledField(PoissonField(poisson), params.χ)
    integrate_vp(particles, params, field, SS, IP, p)
end
```

The march does not write the accelerations `SS.A`, which the field basis reduces. The field of each
sample, deposited from the saved positions and evaluated at them, gives them, as in the script
`scripts/bump_on_tail_training.jl`:

```julia
using VlasovMethods: efield!, update!

w = vec(particles.w)
a = zeros(size(SS.A, 2))
for p in 1:length(pspace), ts in axes(SS.A, 3)
    field = ScaledField(PoissonField(poisson), pspace(p).χ)
    x = @view SS.X[1, :, ts, p]
    update!(field, x, w, 0.0)
    SS.A[1, :, ts, p] .= efield!(field, a, x)
end
```

`h5save(io, SS)` and `Snapshots(io)` write and read the snapshots in HDF5.

## Reduced bases

[`particle_basis`](@ref) stacks the positions over the velocities, one column per saved state and
sample, and passes the snapshot matrix to the `ReducedBasis` of ReducedBasisMethods with the given
algorithm. A cotangent lift, `CotangentLiftSVD()` or `CotangentLiftEVD()`, gives a symplectic
basis. [`field_basis`](@ref) computes the POD basis of the accelerations.

```julia
using ReducedBasisMethods: CotangentLiftEVD, PODEVD, get_DEIM_indices

Φ = particle_basis(CotangentLiftEVD(), SS; tolerance = 1e-9).Φ
Ψₑ = field_basis(PODEVD(), SS; tolerance = 1e-5).Φ
p = get_DEIM_indices(Ψₑ)
```

## The reduced march

[`reduced_integrate_vp`](@ref) marches the reduced coordinates ``\Phi^{\mathsf T} x`` and
``\Phi^{\mathsf T} v`` with the same leapfrog method as the full-order march, and saves the lifted
states into its own snapshots. It takes a [`ReducedIntegratorCache`](@ref) of the basis size,
which one march after another can share:

```julia
IC = ReducedIntegratorCache(IP, size(Φ, 2))
X₀ = Φ' * vec(particles.x)
ROM = Snapshots(Float64, 1, length(particles), nh, nt + 1, pspace)

for s in 1:length(pspace)
    params = pspace(s)
    field = DEIMElectricField(ScaledField(PoissonField(poisson), params.χ), Φ, Ψₑ, p, X₀)
    reduced_integrate_vp(particles, Φ, params, field, ROM, IP, IC, s)
end
```

Both marches are inferred and allocate nothing inside the time loop.

## Diagnostics

The diagnostics take arrays rather than snapshots, so they serve any data of the same shape.

| function | returns |
|:--|:--|
| [`trajectory_errors`](@ref) | the root-mean-square and the maximum errors of the positions and the velocities |
| [`potential_error`](@ref) | the error of the potential coefficients |
| [`field_energy_error`](@ref) | the relative error of the field energy |
| [`energy_drift`](@ref) | the relative drift of an energy history |
| [`growth_rate`](@ref) | the growth rate of a field energy history, fitted to its maxima |
| [`best_and_worst`](@ref) | the samples whose growth rates agree best and worst |
| [`saved_indices`](@ref) | the final, three-quarter and half-way saves of a sweep |

## The analytic quadratic model

The quadratic potential ``H(q, p) = \frac{1}{2} |p|^2 + \frac{1}{2} \beta^2 |q|^2`` has a closed
solution, which tests a reduction without the error of a field solve.
[`AnalyticHamiltonian`](@ref)`(β)` is its Hamiltonian and its gradient, in the form that the
reduced ODE and PODE problems of ReducedBasisMethods take. [`lift_solution`](@ref) lifts a reduced
PODE solution back to the positions and the velocities. `scripts/analytic_quadratic.jl` runs the
experiment.

## The Hamiltonian neural network

[`train_hnn`](@ref) trains a Hamiltonian neural network on the reduced field of a sweep: its
gradient is fitted to the reduced accelerations. It is a package extension, which loads with
`GeometricOptimizers` and `AbstractNeuralNetworks`. `scripts/bump_on_tail_train_hnn.jl` trains it
on the coefficients of the bump-on-tail sweep.
