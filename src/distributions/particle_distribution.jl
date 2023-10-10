
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
xdim(::ParticleDistribution{XD,VD}) where {XD,VD} = XD
vdim(::ParticleDistribution{XD,VD}) where {XD,VD} = VD




# function ParticleDistribution(h5::H5DataStore, path::AbstractString = "/")
#     g = h5[path]
#     name = _name(g)

#     minimum = read(g["minimum"])
#     maximum = read(g["maximum"])
#     samples = read(g["samples"])

#     Parameter(name, minimum, maximum, samples)
# end

# function h5save(h5::H5DataStore, dist::ParticleDistribution; path::AbstractString = "/")
#     n = length(dist.particles)
#     g = _create_group(h5, path)

#     attributes(g)["XD"] = xdim(dist)
#     attributes(g)["VD"] = vdim(dist)

#     z = create_dataset(h5, "z", eltype(dist.particles), ((n, 1), (n, -1)), chunk=(n,1))
#     t = create_dataset(h5, "t", eltype(dist.particles), ((1,), (-1,)), chunk=(1,))

    
#     t[1] = 0
# end

# function h5load(::Type{Parameter}, h5::H5DataStore; path::AbstractString = "/")
#     Parameter(h5, path)
# end

