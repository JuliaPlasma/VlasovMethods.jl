```@meta
CurrentModule = VlasovMethods
```

# VlasovMethods.jl

VlasovMethods implements numerical methods for the Vlasov equation of a collisionless plasma, and
for the collision operators that are added to it. It contains particle methods, spline and grid
discretisations of the distribution function, structure-preserving time integrators, and
reduced-order models that are built from the snapshots of particle runs.

The package serves as the reference implementation of several papers. Each numerical experiment of
these papers has a driver script in `scripts/`, and the section [Scripts](scripts/index.md) of
this manual describes how to run them.

## What the package contains

The models are the equations that the package discretises. Each model holds a distribution
function and the data of its equation, and a method advances it in time.

| model | equation | distribution | time integration |
|:--|:--|:--|:--|
| [`VlasovPoisson`](@ref) | electrostatic Vlasov–Poisson, 1D1V | particles | [`SplittingMethod`](@ref), Strang splitting |
| [`VlasovPoisson`](@ref) | electrostatic Vlasov–Poisson, 1D1V | grid | [`GeometricIntegrator`](@ref), any ODE method |
| [`VlasovWithExternalFields`](@ref) | Vlasov equation in prescribed fields, 3D3V | particles | [`SplittingMethod`](@ref), drift and Boris push |
| [`LenardBernstein`](@ref) | Lenard–Bernstein collision operator | particles and spline | [`GeometricIntegrator`](@ref) |
| [`ConservativeLenardBernstein`](@ref) | conservative Lenard–Bernstein operator | particles and spline | [`GeometricIntegrator`](@ref) |
| [`RescaledConservativeLenardBernstein`](@ref) | the same, rescaled in time | particles and spline | [`GeometricIntegrator`](@ref) |
| [`MetriplecticLenardBernstein`](@ref) | metriplectic Lenard–Bernstein operator | particles and spline | [`Picard_iterate_over_particles!`](@ref) |
| [`Landau`](@ref) | Landau operator, spline method, 2V | particles and spline | [`Picard_iterate_Landau_nls!`](@ref) |
| [`PairwiseLandau`](@ref) | Landau operator, particle-pair method, 2V | particles and spline | [`discrete_gradient_step!`](@ref) |

Around the models, the package provides

  - three representations of the distribution function, and the projections between them;
  - initial conditions, from the Maxwellian to the bump-on-tail and the BKW solution;
  - the reduced tensors of a grid discretisation, for a Galerkin reduction of the grid model;
  - snapshots, reduced bases, reduced and hyper-reduced electric fields and the reduced march of
    the particle Vlasov–Poisson model;
  - the frontends of a multi-modal symplectic autoencoder;
  - plotting functions for every figure of the papers.

## How this manual is organised

The sections introduce the package from the bottom up.

  - [**Getting started**](manual/getting_started.md) runs a complete Vlasov–Poisson simulation
    and a complete collision simulation.
  - **Distribution functions**: [Representations](manual/distributions.md),
    [Initial conditions](manual/initial_conditions.md) and [Projections](manual/projections.md)
    describe how a distribution function is stored, how it is initialised, and how the
    representations are converted into each other.
  - **Collisionless models**: [Running a model](manual/methods.md) describes the methods, the time
    integrators and the files that a run writes; [Vlasov–Poisson](manual/vlasov_poisson.md) and
    [External electromagnetic fields](manual/external_fields.md) describe the models.
  - **Collision operators**: [Lenard–Bernstein collisions](manual/lenard_bernstein.md) and
    [Landau collisions](manual/landau.md).
  - **Model reduction**: [Grid-based reduced tensors](manual/grid_tensors.md),
    [Reduced-order models](manual/reduced_order.md) and
    [Multi-modal autoencoder frontends](manual/autoencoder.md).
  - [**Plotting**](manual/plotting.md) lists the figures that the package draws.
  - [**Scripts**](scripts/index.md) gives, for each paper, the command of every experiment, the
    parameters, and the results of a run.
  - The [**Library**](library/index.md) holds the docstring of every documented function and
    type, ordered by the sections above.

## Installation

VlasovMethods is registered in the General registry. Install it with

```julia
using Pkg
Pkg.add("VlasovMethods")
```

The plotting functions are defined in a package extension. They become available when
`CairoMakie` is loaded. The adapter for the fields of `ElectromagneticFields`, the trainer of the
Hamiltonian neural network and the autoencoder frontends are package extensions, too; each page
names the packages that load its extension.

## The papers

The reduced-order models reproduce the experiments of

  - T. M. Tyranowski and M. Kraus, *Symplectic model reduction methods for the Vlasov equation*,
    Contributions to Plasma Physics 63, e202200046 (2023),
    [doi:10.1002/ctpp.202200046](https://doi.org/10.1002/ctpp.202200046),
    [arXiv:1910.06026](https://arxiv.org/abs/1910.06026).

The collision operators implement the particle methods of three manuscripts. This manual names them
by their abbreviations:

  - **PD-LB**, the particle discretisation of the Lenard–Bernstein operator;
  - **MP-LB**, the metriplectic particle methods for the Lenard–Bernstein operator;
  - **MP-L**, the metriplectic particle methods for the Landau operator.

Where a docstring names an equation label such as `eq:velocity_ode`, it is the label of that
equation in the manuscript of its operator.
