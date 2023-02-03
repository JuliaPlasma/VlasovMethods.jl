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

import GeometricEquations
import GeometricEquations: ntime
import GeometricIntegrators.Integrators


# abstract types

include("distributions/distribution.jl")
include("projections/projection.jl")
include("methods/method.jl")
include("models/model.jl")
include("examples/example.jl")
include("sampling/sampling.jl")


export initialize!

# distribution functions

include("distributions/particle_distribution.jl")

export ParticleDistribution


# projections

include("projections/potential.jl")


# numerical methods

include("methods/splitting.jl")

export run!
export SplittingMethod


# Vlasov models

include("models/vlasov_poisson.jl")

export VlasovPoisson


# Example Problems

include("examples/bumpontail.jl")
include("examples/normal.jl")
include("examples/twostream.jl")

export BumpOnTail, NormalDistribution




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
