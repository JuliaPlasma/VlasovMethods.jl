module VlasovMethods

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

export projection

# numerical methods

include("methods/splitting.jl")
include("methods/diffeq_integrator.jl")

export run!
export run
export SplittingMethod
export DiffEqIntegrator


# Vlasov models

include("models/vlasov_poisson.jl")
include("models/lenard_bernstein.jl")

export VlasovPoisson
export LenardBernstein
export LB_rhs

# Example Problems

include("examples/bumpontail.jl")
include("examples/normal.jl")
include("examples/twostream.jl")
include("examples/shiftednormalv.jl")


export BumpOnTail, NormalDistribution, ShiftedNormalV


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
