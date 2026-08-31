module VlasovMethods

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Distances
using FastGaussQuadrature
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
using QuadGK
using Random
using Sobol
using SimpleSolvers
using SpecialFunctions
using StaticArrays
using Trapz

# import DifferentialEquations
# import NaNMath

import Base: Callable

import BSplineKit
import BSplineKit: AbstractBSplineBasis, BSplineBasis, PeriodicBSplineBasis,
                   RecombinedBSplineBasis
import BSplineKit: BSplineOrder
import BSplineKit.BSplines

import GeometricEquations
import GeometricEquations: ntime
import GeometricIntegrators.Integrators
import GeometricIntegrators.Extrapolators

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
include("splines/nd_spline.jl")

export initialize!

# splines

include("splines/spline_nd.jl")

export SplineND, L2projection!

# include("splines/2d_spline.jl")
include("splines/2d_spline_new.jl")
export TwoDSpline
export evaluate, evaluate_first_derivative

# distribution functions

include("distributions/maxwellian.jl")
include("distributions/particle_distribution.jl")
include("distributions/spline_distribution.jl")

export ParticleDistribution
export SplineDistribution

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
