module VlasovMethods

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Distances
using HDF5
using LinearAlgebra
using LinearSolve
using NonlinearSolve
using OffsetArrays
using Parameters
using ParticleMethods
# using Plots
using PoissonSolvers
using ProgressMeter
using QuadratureRules
using Random
using Sobol
using SimpleSolvers
using SpecialFunctions
using StaticArrays
using Trapz

# import DifferentialEquations
# import NaNMath

import Base: Callable

using SimpleSplines
import SimpleSplines: basis, coefficients, derivative, evaluate, mass_matrix, mass_operator

import GeometricEquations
import GeometricEquations: ntime

# `import A.B` binds only `B`, so the bare name has to be imported in its own right: four
# call sites qualify names as `GeometricIntegrators.…` and were throwing `UndefVarError`.
import GeometricIntegrators
import GeometricIntegrators.Integrators

# The `Extrapolators` submodule was flattened in GeometricIntegrators 0.18 — the names now
# live in GeometricIntegratorsBase and are re-exported at the top level. `import
# GeometricIntegrators.Extrapolators` therefore only warns rather than failing, and leaves
# the binding undefined, so the eight call sites throw at run time instead of at load.
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

# The spline machinery lives in SimpleSplines. Everything the four files that used to sit in
# `src/splines/` provided — the tensor-product basis, its mass matrix, the L² projection, the
# local evaluation of a basis and its gradient at a point — is there, with per-axis degrees,
# domains and boundary conditions, and with the Kronecker structure of the mass operator used
# rather than assembled.
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
export evaluate_all, evaluate_all!, basis_index, local_width
export quadrature_nodes, quadrature_weights, basis_values
export l2_projection, l2_projection!, mass_operator, mass_matrix, mass_solve!

# distribution functions

include("distributions/maxwellian.jl")
include("distributions/particle_distribution.jl")
include("distributions/spline_distribution.jl")

export ParticleDistribution
export SplineDistribution
export check_conservation_basis, project_function, project_Maxwellian

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
include("methods/diffeq_integrator.jl")
include("methods/geometric_integrator.jl")
include("methods/Landau_solver.jl")

export run!
export run
export SplittingMethod
export DiffEqIntegrator
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

export VlasovPoisson
export LenardBernstein
export ConservativeLenardBernstein
export RescaledConservativeLenardBernstein
export MetriplecticLenardBernstein
export Landau

# Example Problems

include("examples/bumpontail.jl")
include("examples/normal.jl")
include("examples/uniform.jl")
include("examples/twostream.jl")
include("examples/shiftednormalv.jl")
include("examples/shifteduniform.jl")
include("examples/doublemaxwellian.jl")
include("examples/bump.jl")
include("examples/sum_maxwellian.jl")

export BumpOnTail, NormalDistribution, UniformDistribution, ShiftedNormalV,
       ShiftedUniformDistribution, DoubleMaxwellian, Bump, SumMaxwellian

# include("electric_field.jl")

# export ElectricField, PoissonField, ExternalField
# export ScaledField, ScaledPoissonField, ScaledExternalField

# include("vlasov_poisson.jl")

# export VPIntegratorParameters, VPIntegratorCache, integrate_vp!

# include("sampling.jl")

# export draw_g_accept_reject, draw_g_importance_sampling, weight_f

# include("visualisation.jl")

# export plot_particles, plot_distribution

end
