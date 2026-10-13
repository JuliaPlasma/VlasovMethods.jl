module VlasovMethods

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using HDF5
using LinearAlgebra
using MultiIndexArrays: multiindex, linearindex, GridTensor, MultiIndexAxis,
                        MultiIndexLazyArray
using Parameters
using ParticleMethods
using PoissonSolvers
using ProgressMeter
using Random
using ReducedBasisMethods: ReducedTensor
using Sobol
using SpecialFunctions
using StaticArrays

using SimpleSplines
import SimpleSplines: basis, coefficients, derivative, evaluate, mass_matrix, mass_operator

import GeometricEquations
import GeometricEquations: ntime

# `import` and not `using`: `SimpleSolvers` exports generic names (`Options`, `Newton`,
# `Picard`, `solve`), and every use here is qualified as `SimpleSolvers.…`.
import SimpleSolvers

using GeometricBrackets: Arakawa, poisson_apply

# `import A.B` binds only `B`, so the bare module name needs an import of its own for the call
# sites that qualify names as `GeometricIntegrators.…`.
import GeometricIntegrators
import GeometricIntegrators.Integrators

# The `Extrapolators` submodule was flattened in GeometricIntegrators 0.18: the names live in
# GeometricIntegratorsBase and are re-exported at the top level, so they are imported from
# there rather than through a submodule path.
import GeometricIntegrators: extrapolate!, HermiteExtrapolation, MidpointExtrapolation

# utilities

include("cache.jl")

# abstract types

include("distributions/distribution.jl")
include("projections/projection.jl")
include("methods/method.jl")
include("models/model.jl")
include("examples/example.jl")
include("sampling/sampling.jl")
include("entropies/entropy.jl")

export initialize!

# The spline machinery lives in SimpleSplines: the tensor-product basis, its mass matrix, the
# L² projection and the local evaluation of a basis and its gradient at a point, with per-axis
# degrees, domains and boundary conditions, and with the Kronecker structure of the mass
# operator used rather than assembled. These names are re-exported so that a caller working
# with a `SplineDistribution` does not have to depend on SimpleSplines directly.
export Spline, derivative
export BSplineBasis, PeriodicBSplineBasis, RecombinedBSplineBasis, TensorProductBasis
export UniformMesh, GradedMesh, RandomMesh, GeneralMesh
export Free, Periodic, Dirichlet, Neumann, Natural, Robin, Constraint
export polynomial_reproduction
# `..` builds a domain and every mesh constructor takes one, so it is re-exported for the same
# reason SimpleSplines re-exports it: `UniformMesh(40, -10 .. 10)` failing with `.. not
# defined` is a poor first experience.
export ..
export ncells, nbasis, degree, order, breakpoints, meshwidth, domain
export evaluate, evaluate_all, evaluate_all!, basis_index, local_width
export quadrature_nodes, quadrature_weights, basis_values
export l2_projection, l2_projection!, mass_operator, mass_matrix, mass_solve!

# distribution functions

include("distributions/maxwellian.jl")
include("distributions/particle_distribution.jl")
include("distributions/spline_distribution.jl")
include("distributions/grid_distribution.jl")

export ParticleDistribution
export SplineDistribution
export GridDistribution, velocity_moments
export check_conservation_basis, project_function, project_Maxwellian

# plotting

include("plotting.jl")

export plot_distribution, plot_particles

# entropy models

include("entropies/collision_entropy.jl")

export CollisionEntropy

# projections

include("projections/potential.jl")
include("projections/distribution.jl")
include("projections/density.jl")

export projection

# numerical methods

include("methods/splitting.jl")
include("methods/geometric_integrator.jl")
include("methods/Landau_solver.jl")

export run!
export SplittingMethod
export GeometricIntegrator
export Picard_iterate_over_particles

# Vlasov models

include("models/vlasov_model.jl")

include("models/collision_operator.jl")
include("models/landau.jl")
include("models/lenard_bernstein.jl")
include("models/lenard_bernstein_conservative.jl")
include("models/lenard_bernstein_metriplectic.jl")
include("models/rescaled_lenard_bernstein_conservative.jl")
include("models/vlasov_poisson.jl")

# The `ElectromagneticFields` adapter lives in `ext/VlasovMethodsElectromagneticFieldsExt.jl`
# and loads when both packages are loaded. The field interface is defined by the included model.
include("models/vlasov_with_external_fields.jl")

export VlasovPoisson
export VlasovWithExternalFields, UniformField, electric_field, magnetic_field
export sodeensemble
export LenardBernstein
export ConservativeLenardBernstein
export RescaledConservativeLenardBernstein
export MetriplecticLenardBernstein
export Landau

# Example Problems

include("examples/bumpontail.jl")
include("examples/normal.jl")
include("examples/uniform.jl")
include("examples/shiftednormalv.jl")
include("examples/shifteduniform.jl")
include("examples/doublemaxwellian.jl")
include("examples/bump.jl")
include("examples/sum_maxwellian.jl")

export BumpOnTail, NormalDistribution, UniformDistribution, ShiftedNormalV,
       ShiftedUniformDistribution, DoubleMaxwellian, Bump, SumMaxwellian

# grid-based phase-space discretisation

include("gridbased/moments.jl")

export _apply_∫dv!

include("gridbased/reduced_tensors.jl")

export bracket_tensor, PotentialReducedTensor, VelocityReducedMatrix, FullyReducedTensor

include("gridbased/collisions.jl")

# The grid 1D1V Vlasov–Poisson scheme adds methods whose signatures name `VlasovPoisson` and
# `GeometricIntegrator`, so it is included after both.
include("gridbased/vlasov_poisson_1d1v.jl")

# The particle-based reduction in `src/particles/` is not included: each of the four files
# names a binding that neither the module nor its dependencies define —
# `PoissonSolverPBSplines`, `PBSpline`, `ElectricField`, `ParameterSpace` — so including one
# breaks the load.

end
