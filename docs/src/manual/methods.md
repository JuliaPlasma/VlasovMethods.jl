```@meta
CurrentModule = VlasovMethods
```

# Running a model

A simulation in VlasovMethods has three parts.

 1. A **distribution function** holds the state: the particles, the spline coefficients or the
    grid values. [Distribution functions](distributions.md) describes the three types.
 2. A **model** holds the distribution and the data of its equation: the potential of a
    Vlasov–Poisson problem, the external field, or the spline basis and the collision frequency
    of a collision operator.
 3. A **method** advances the model in time. It is built from the model, the time span and the
    time step.

The methods are thin wrappers around
[GeometricIntegrators](https://github.com/JuliaGNI/GeometricIntegrators.jl). A method holds the
`GeometricEquations` problem of its model and the integrator that solves it. [`run!`](@ref)
advances the integrator over the whole time span, writes every state into an HDF5 file, and
copies the final state back into the distribution of the model.

## The two method types

| method | built by | problem | state | models |
|:--|:--|:--|:--|:--|
| [`SplittingMethod`](@ref) | `SplittingMethod(model, tspan, tstep)` | `SODEProblem` | `dist.particles.z`, positions and velocities | particle [`VlasovPoisson`](@ref), [`VlasovWithExternalFields`](@ref) |
| [`GeometricIntegrator`](@ref) | `GeometricIntegrator(model, tspan, tstep[, method])` | `ODEProblem` | `dist.particles.v[1, :]` or `vec(dist.values)` | [`LenardBernstein`](@ref), [`ConservativeLenardBernstein`](@ref), [`RescaledConservativeLenardBernstein`](@ref), grid [`VlasovPoisson`](@ref) |

A `SplittingMethod` composes the exact flows of the parts of a split vector field with the Strang
composition of GeometricIntegrators. A `GeometricIntegrator` takes the vector field of the model as
one `ODEProblem` and solves it with any ODE method of GeometricIntegrators: the default is
`RK438()` for [`LenardBernstein`](@ref) and `ImplicitMidpoint()` for the conservative collision
operators and the grid Vlasov–Poisson model.

!!! note "Two functions of the same name"
    `GeometricIntegrators` exports a function `GeometricIntegrator` as well. Where both packages
    are loaded with `using`, write `VlasovMethods.GeometricIntegrator(model, tspan, tstep)`.

The models of the metriplectic and the Landau operators are not advanced by `run!`. Each has a
step function, which the caller calls once per time step:

| model | step function | scheme |
|:--|:--|:--|
| [`MetriplecticLenardBernstein`](@ref) | [`Picard_iterate_over_particles!`](@ref) | implicit midpoint rule, Picard iteration |
| [`Landau`](@ref) | [`Picard_iterate_Landau_nls!`](@ref) | implicit midpoint rule, Picard iteration |
| [`PairwiseLandau`](@ref) | [`discrete_gradient_step!`](@ref) | discrete gradient method, fixed-point iteration |

The pages of these models show the loop.

## The output of `run!`

`run!(method, h5file)` creates the file `h5file`, or overwrites it. It writes one state per time
of the range `t₀:Δt:t₁`, the initial state included.

  - A `SplittingMethod` writes the dataset `z` of the size `(nd, npart, ntime)`. `z[:, :, s]` is
    the state `dist.particles.z` at the time slice `s`: the positions in the first rows and the
    velocities in the last rows.
  - A `GeometricIntegrator` writes the dataset `z` of the size `(n, ntime)`, one column per time,
    and the dataset `t` of the times.

The run holds one state in memory at a time, so the length of a run is limited by the disk, not by
the memory. After the last step, `run!` copies the final state into the distribution of the model
and returns that distribution.

A time span that is not a multiple of the step ends at the last time of the range `t₀:Δt:t₁`,
before `t₁`.

## Ensembles

[`sodeensemble`](@ref) builds the training data of a reduced model: one Vlasov–Poisson problem per
member of a list of particle distributions, as a `GeometricEquations.SODEEnsemble` that
`GeometricIntegrators.integrate` solves in one call.
[Vlasov–Poisson](vlasov_poisson.md#Ensembles-of-runs) shows how to use it.
