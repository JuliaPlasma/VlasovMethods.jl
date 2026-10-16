# The field of check 1: the analytic quadratic model's acceleration `-β²x`, as a test-only
# `ElectricField`. It is wrapped as `ReducedElectricField(field, Φ, X)` with the cotangent-lift `Φ`,
# which is the one basis that both the positions and the velocities use.
#
# `update!` is a no-op, because the model's field is a function of the position alone and there is
# no potential to solve for; `coefficients` returns an empty `Vector{T}`, which `Snapshots` allows at
# `nh = 0`; and `energy` is zero, the field being external to the particles.

using VlasovMethods: ElectricField
# `import` and not `using`: the methods below extend the field interface.
import VlasovMethods: update!, efield!, energy, coefficients

struct QuadraticField{T} <: ElectricField
    β::T
end

update!(::QuadraticField, x, w, t) = nothing

function efield!(f::QuadraticField, e, x)
    β² = f.β^2
    @inbounds for i in eachindex(e)
        e[i] = -β² * x[i]
    end
    return e
end

energy(::QuadraticField{T}) where {T} = zero(T)

coefficients(::QuadraticField{T}) where {T} = T[]
