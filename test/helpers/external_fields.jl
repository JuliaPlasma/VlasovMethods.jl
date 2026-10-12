# The fixtures and checks the external-field test files share. It is included by
# `test/models/vlasov_with_external_fields.jl`, `test/integration/electromagnetic_fields_ext.jl`,
# `test/methods/splitting.jl` and `test/quality/jet.jl`, and is listed in no `runtests.jl`.

using LinearAlgebra: ×, ⋅
using ParticleMethods: ParticleList
using StaticArrays
using Test
using VlasovMethods
using VlasovMethods: s_drift!, s_lorentz!, v_drift!, v_lorentz!

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

# `½|v|²` of the first particle stays within the rounding of `N` steps of `½|v₀|²` at every time
# point of the solution. The Boris push is a rotation where `E = 0`, which keeps `|v|` exactly in
# exact arithmetic; an explicit kick grows the norm by a factor per step instead and leaves this
# bound far behind.
function test_energy(sol, v₀)
    T = eltype(v₀)
    N = lastindex(sol.q)
    E₀ = v₀ ⋅ v₀ / 2
    drift = zero(T)
    for n in 0:N
        v = sol.q[n]
        drift = max(drift, abs((v[4, 1]^2 + v[5, 1]^2 + v[6, 1]^2) / 2 - E₀))
    end
    @test drift ≤ N * eps(T) * E₀
end

# `s_drift!` and `s_lorentz!` of a model in `field`, with three particles, `q = 2` and `m = 3`,
# are type-stable and allocate nothing at the model's concrete argument types. A runtime dispatch
# is what the allocation figures measure only indirectly.
function test_push_does_not_allocate(::Type{T}, field) where {T}
    model = VlasovWithExternalFields(
        particle_distribution(T, 3), field; charge = T(2), mass = T(3))
    method = SplittingMethod(model, (zero(T), T(1)), T(0.1))
    params = method.equation.parameters
    z̄ = copy(model.distribution.particles.z)
    z = similar(z̄)
    t̄ = zero(T)
    @test @inferred(s_drift!(z, t̄ + T(0.1), z̄, t̄, params)) === nothing
    @test @inferred(s_lorentz!(z, t̄ + T(0.1), z̄, t̄, params)) === nothing
    if Base.JLOptions().code_coverage == 0
        @test allocations_s(s_drift!, z, t̄ + T(0.1), z̄, t̄, params) == 0
        @test allocations_s(s_lorentz!, z, t̄ + T(0.1), z̄, t̄, params) == 0
    end
    # a control: the barrier machinery sees an allocation, so the zeros above are not vacuous
    @test allocations(collect, z) > 0
end

# The right-hand side of each of the two maps of a model in `field`, against values written out
# here: two particles at `x` with velocity `v`, `q = 2`, `m = 3`, at `t = 0.4`. The orbit, energy
# and allocation checks cannot see a zero `ẋ` in the drift, or a `qm` dropped from the Lorentz
# acceleration: both leave them exactly as they are. The two controls require a nonzero `E` and a
# nonzero acceleration at `x`, so that a dropped `E` or `q/m` shows up in the acceleration.
function test_vector_fields(::Type{T}, field, x, v) where {T}
    q, m, t = T(2), T(3), T(0.4)
    dist = particle_distribution(T, 2)
    dist.particles.x .= x
    dist.particles.v .= v
    model = VlasovWithExternalFields(dist, field; charge = q, mass = m)
    method = SplittingMethod(model, (zero(T), T(1)), T(0.1))
    params = method.equation.parameters
    z = copy(model.distribution.particles.z)
    ż = similar(z)
    a = (q / m) * (electric_field(field, t, x) + v × magnetic_field(field, t, x))

    @test !all(iszero, a)
    @test !all(iszero, electric_field(field, t, x))

    @test @inferred(v_drift!(ż, t, z, params)) === nothing
    @test @inferred(v_lorentz!(ż, t, z, params)) === nothing

    v_drift!(ż, t, z, params)
    for i in axes(z, 2)
        @test ż[1:3, i] == v
        @test all(iszero, ż[4:6, i])
    end

    v_lorentz!(ż, t, z, params)
    for i in axes(z, 2)
        @test all(iszero, ż[1:3, i])
        @test ż[4:6, i] == collect(a)
    end
end
