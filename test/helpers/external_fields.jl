# The fixtures the external-field test files share. It is included by
# `test/models/vlasov_with_external_fields.jl`, `test/integration/electromagnetic_fields_ext.jl`
# and `test/quality/jet.jl`, and is listed in no `runtests.jl`.

using ParticleMethods: ParticleList
using StaticArrays
using VlasovMethods

# A `ParticleDistribution{T,3,3}` from a matrix of element type `T`, as
# `test/integration/vlasov_poisson.jl` builds the `Float32` case: the two-argument constructor
# always fills `Float64` zeros.
function particle_distribution(::Type{T}, npart) where {T}
    z = zeros(T, 7, npart)
    vars = (x = 1:3, v = 4:6, z = 1:6, w = 7:7)
    return ParticleDistribution(3, 3, ParticleList(z; variables = vars))
end

# One particle at `x₀` with velocity `v₀`.
function one_particle(::Type{T}, x₀, v₀) where {T}
    dist = particle_distribution(T, 1)
    dist.particles.x .= SVector{3, T}(x₀)
    dist.particles.v .= SVector{3, T}(v₀)
    return dist
end

# `@allocated` through a fixed-arity barrier whose arguments have concrete types, so that a Julia
# 1.11 closure or splat boxing is not counted. It calls the function once before it measures.
# Coverage instrumentation allocates on every counted line, so the assertions that use it do not
# run under it rather than being weakened.
function allocations_s(f, z, t, z̄, t̄, params)
    f(z, t, z̄, t̄, params)
    return @allocated f(z, t, z̄, t̄, params)
end

# A generic barrier for the control: it shows that the machinery sees an allocation, so a zero
# from the fixed-arity barriers above is not vacuous.
function allocations(f, a)
    f(a)
    return @allocated f(a)
end
