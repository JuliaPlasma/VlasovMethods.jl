
abstract type VlasovExample end

"""
    initialize!(dist::ParticleDistribution, example::VlasovExample[, sampling])

Draw a particle sample of the initial condition `example` into `dist`, and return `dist`.

Every example writes the velocities and the weights of the particles. Some write the positions
as well, and some leave them as they are; the docstring of each example says which.
`sampling` selects the sampling method of the examples that sample by rejection: `BumpOnTail`
takes `AcceptRejectSampling()`, its default, or `ImportanceSampling()`, and `Bump` takes
`AcceptRejectSampling()` only. The other examples ignore it.
"""
function initialize! end
