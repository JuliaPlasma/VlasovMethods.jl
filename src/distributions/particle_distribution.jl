
struct ParticleDistribution{XD, VD, PT} <: DistributionFunction{XD,VD}
    particles::PT

    function ParticleDistribution(xdim, vdim, particles::PT) where {PT <: ParticleList}
        new{xdim, vdim, PT}(particles)
    end
end

function ParticleDistribution(xdim, vdim, npart::Int; kwargs...)
    z = zeros(xdim+vdim+1, npart)
    vars = ( 
            x = 1:xdim,
            v = xdim+1:xdim+vdim,
            z = 1:xdim+vdim,
            w = xdim+vdim+1:xdim+vdim+1,
           )
    list = ParticleList(z; variables = vars, kwargs...)
    ParticleDistribution(xdim, vdim, list)
end


Base.size(dist::ParticleDistribution) = size(dist.particles)
