module VlasovMethods

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using BSplineKit
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
using Random

import GeometricEquations
import GeometricEquations: ntime
import GeometricIntegrators.Integrators
import DifferentialEquations


# abstract types

include("distributions/distribution.jl")
include("projections/projection.jl")
include("methods/method.jl")
include("models/model.jl")
include("examples/example.jl")
include("sampling/sampling.jl")
include("entropies/entropy.jl")


export initialize!

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

export run!
export run
export SplittingMethod
export DiffEqIntegrator
export GeometricIntegrator


# Vlasov models

include("models/vlasov_poisson.jl")
include("models/lenard_bernstein.jl")
include("models/conservative_lb.jl")

export VlasovPoisson
export LenardBernstein
export ConservativeLenardBernstein
export LB_rhs
export CLB_rhs

# Example Problems

include("examples/bumpontail.jl")
include("examples/normal.jl")
include("examples/uniform.jl")
include("examples/twostream.jl")
include("examples/shiftednormalv.jl")
include("examples/shifteduniform.jl")
include("examples/doublemaxwellian.jl")


export BumpOnTail, NormalDistribution, UniformDistribution, ShiftedNormalV, ShiftedUniformDistribution, DoubleMaxwellian


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
