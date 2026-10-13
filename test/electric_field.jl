using LinearAlgebra
using ParticleMethods: ParticleList
using PoissonSolvers
using Random
using SimpleSplines: SplineQuadrature, stiffness_matrix
using Test
using VlasovMethods
using VlasovMethods: ElectricField, PoissonField, coefficients, efield!, energy,
                     projection!,
                     update!
using VlasovMethods: _electric_field

# `@allocated` through a fixed-arity barrier whose arguments have concrete types, so that a Julia
# 1.11 closure or splat boxing is not counted. Each call runs once before it is measured.
function _alloc_update(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    (f(a, b, c, d); @allocated f(a, b, c, d))
end
function _alloc_efield(f::F, a::A, b::B, c::C) where {F, A, B, C}
    (f(a, b, c); @allocated f(a, b, c))
end

# The particle distribution of `x` and `w`, for the deposit of `projection!`: the same two
# columns the field's own deposit takes as vectors.
function particle_distribution(::Type{T}, x::AbstractVector, w::AbstractVector) where {T}
    n = length(x)
    z = zeros(T, 3, n)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    d = ParticleDistribution(1, 1, ParticleList(z; variables = vars))
    d.particles.x .= x'
    d.particles.w .= w'
    return d
end

# The field of `n` seeded positions in `(0, 1)` with equal weights, on `ncells` periodic spline
# cells of degree `order`.
function poisson_field(::Type{T}; n = 40, ncells = 16, order = 3) where {T}
    Random.seed!(1234)
    q = T(0.2) .+ T(0.6) .* rand(T, n)
    w = fill(T(1) / n, n)
    potential = Potential(PeriodicBasisSpline((zero(T), one(T)), order, ncells))
    return PoissonField(potential), q, w
end

@testset "ElectricField and PoissonField, $T" for T in (Float64, Float32)
    @testset "the interface" begin
        field, q, w = poisson_field(T)
        @test field isa ElectricField
        @test field.potential isa Potential
        # the buffers are of the element type of the potential's right-hand side
        @test eltype(field.work) === T
        @test eltype(field.nodes) === T
        @test eltype(field.weights) === T

        # `f(e, x, w, t)` is `update!` followed by `efield!`, and `f(e, x)` is `efield!` alone
        # on the potential that deposits
        e1 = similar(q)
        @test field(e1, q, w, zero(T)) === e1
        e2 = similar(q)
        @test field(e2, q) === e2
        @test e1 == e2
        @test field(q, w, zero(T)) == e1

        # `update!` ignores `t`: the field is time-independent
        @test update!(field, q, w, T(3.5)) === nothing
    end

    @testset "the deposit matches projection!" begin
        field, q, w = poisson_field(T)
        update!(field, q, w, zero(T))

        # the same deposit through the matrix `_deposit!` of `projection!`
        p = field.potential
        other = Potential(basis(p))
        projection!(other, particle_distribution(T, q, w))
        PoissonSolvers.update!(other)

        @test PoissonSolvers.rhs(p) == PoissonSolvers.rhs(other)
        @test coefficients(field) == coefficients(other)
        @test coefficients(field) === PoissonSolvers.coefficients(p)
    end

    @testset "the field matches _electric_field" begin
        field, q, w = poisson_field(T)
        update!(field, q, w, zero(T))
        e = similar(q)
        @test efield!(field, e, q) === e
        work = similar(field.work)
        @test e == [_electric_field(field.potential, work, x) for x in q]
    end

    @testset "the field energy matches the stiffness matrix" begin
        field, q, w = poisson_field(T)
        update!(field, q, w, zero(T))
        b = basis(field.potential)
        S = stiffness_matrix(SplineQuadrature(b))
        c = coefficients(field)
        # `∫(∂ₓϕ)² = ϕᵀSϕ`, assembled two ways: the quadrature the constructor tabulated and the
        # stiffness matrix of the basis
        @test energy(field) ≈ dot(c, S, c) / 2 rtol = 10sqrt(eps(T))
    end

    @testset "update! and efield! are inferred and allocate nothing" begin
        field, q, w = poisson_field(T; n = 200)
        e = similar(q)
        t = zero(T)

        @test @inferred(update!(field, q, w, t)) === nothing
        @test @inferred(efield!(field, e, q)) === e
        @test @inferred(field(e, q, w, t)) === e

        # coverage instrumentation allocates on every counted line, so under it the assertion
        # does not run rather than being weakened
        if Base.JLOptions().code_coverage == 0
            @test _alloc_update(update!, field, q, w, t) == 0
            @test _alloc_efield(efield!, field, e, q) == 0
            @test _alloc_update(field, e, q, w, t) == 0
        end
        # the controls: the barrier sees an allocation when there is one
        @test _alloc_update((f, x, w, t) -> copy(x), field, q, w, t) > 0
        @test _alloc_efield((f, e, x) -> similar(e), field, e, q) > 0
    end
end
