module VlasovMethods

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using BSplineKit
using BSplineKit.BSplines
using HDF5
using OffsetArrays
using Parameters
using ParticleMethods
using PoissonSolvers
using ProgressMeter
using Random
using Sobol
using SpecialFunctions
using LinearAlgebra
using QuadGK
using HCubature
using Plots
using SimpleSolvers

import GeometricEquations
import GeometricEquations: ntime
import GeometricIntegrators.Integrators
import GeometricIntegrators.Extrapolators
# import DifferentialEquations


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
# include("splines/2d_spline.jl")
include("splines/2d_spline_new.jl")
export TwoDSpline
export evaluate, evaluate_first_derivative

# distribution functions

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

include("models/vlasov_poisson.jl")
include("models/lenard_bernstein.jl")
include("models/conservative_lb.jl")
include("models/landau.jl")

export VlasovPoisson
export LenardBernstein
export ConservativeLenardBernstein
export Landau
export LB_rhs
export CLB_rhs
export compute_matrices, compute_U, Landau_rhs!, stuff, compute_L

# Example Problems

include("examples/bumpontail.jl")
include("examples/normal.jl")
include("examples/uniform.jl")
include("examples/twostream.jl")
include("examples/shiftednormalv.jl")
include("examples/shifteduniform.jl")
include("examples/doublemaxwellian.jl")
include("examples/bump.jl")


export BumpOnTail, NormalDistribution, UniformDistribution, ShiftedNormalV, ShiftedUniformDistribution, DoubleMaxwellian, Bump


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
