@doc raw"""
    ParticleDistribution(xdim, vdim, npart; kwargs...)
    ParticleDistribution(xdim, vdim, particles::ParticleList)

The distribution function as the Klimontovich measure of ``N`` particles,

```math
f_p(x, v) = \sum_{\alpha = 1}^{N} w_\alpha \, \delta(x - x_\alpha) \, \delta(v - v_\alpha) ,
```

with `xdim` position and `vdim` velocity components per particle. It is a
`DistributionFunction{DT, xdim, vdim}`, beside `SplineDistribution` and `GridDistribution`.

The field `particles` is a `ParticleList` of ParticleMethods with the variables `x`, the
`xdim × N` positions, `v`, the `vdim × N` velocities, `z`, the `(xdim + vdim) × N` phase-space
coordinates, positions first, and `w`, the `1 × N` weights. All four are views of one
`(xdim + vdim + 1) × N` matrix, so a write to one of them changes the others.

The first form allocates `npart` particles of `Float64` zeros and passes `kwargs` to the
`ParticleList` constructor. The second form wraps an existing list, whose element type is the
element type of the distribution.

```julia
dist = ParticleDistribution(1, 1, 1000)
dist.particles.w .= 1 / 1000
```
"""
struct ParticleDistribution{DT, XD, VD, PT <: ParticleList{DT}} <:
       DistributionFunction{DT, XD, VD}
    particles::PT

    function ParticleDistribution(xdim, vdim, particles::PT) where {PT <: ParticleList}
        new{eltype(particles), xdim, vdim, PT}(particles)
    end
end

function ParticleDistribution(xdim, vdim, npart::Int; kwargs...)
    z = zeros(xdim+vdim+1, npart)
    vars = (
        x = 1:xdim,
        v = (xdim + 1):(xdim + vdim),
        z = 1:(xdim + vdim),
        w = (xdim + vdim + 1):(xdim + vdim + 1)
    )
    list = ParticleList(z; variables = vars, kwargs...)
    ParticleDistribution(xdim, vdim, list)
end

Base.eltype(::ParticleDistribution{T}) where {T} = T
Base.size(dist::ParticleDistribution) = size(dist.particles)

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
