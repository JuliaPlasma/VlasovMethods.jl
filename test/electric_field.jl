using LinearAlgebra
using ParticleMethods: ParticleList
using PoissonSolvers
using SimpleSplines: SplineQuadrature, stiffness_matrix
using Test
using VlasovMethods
using VlasovMethods: ElectricField, PoissonField, coefficients, efield!, energy,
                     projection!,
                     update!
using VlasovMethods: _electric_field

include("helpers/allocation_barriers.jl")
include("helpers/poisson_field.jl")

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

# The field of the positions that reach the last cell of the domain, on a `DirichletBasisSpline`:
# a recombined basis leaves padding in the block of basis values at the two ends, which the
# positions in `(0.2, 0.8)` that `poisson_field` draws never reach.
function edge_poisson_field(::Type{T}) where {T}
    q = T.([0.05, 0.5, 0.93, 0.97, 0.995, 0.9999])
    w = fill(T(1) / length(q), length(q))
    return PoissonField(Potential(DirichletBasisSpline((zero(T), one(T)), 3, 16))), q, w
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

        # the `e` that `f(x, w, t)` allocates is of the element type of `x`: the comparison
        # above crosses types and would not notice a buffer of another type
        @test eltype(field(q, w, zero(T))) === T

        # `update!` ignores `t`: the field is time-independent
        @test update!(field, q, w, T(3.5)) === nothing
    end

    # A recombined basis has no single width: a function near an end spans the union of two parent
    # supports, so `local_width` is the largest block and a cell with fewer nonzero functions
    # leaves the tail of the block zero. Those padding entries carry indices past `nbasis`, so a
    # deposit that added them instead of skipping them raises a `BoundsError`. The positions of
    # `edge_poisson_field` reach the last cell, which is where that happens.
    @testset "the deposit matches projection!, $(nameof(fixture))" for fixture in (
        poisson_field, edge_poisson_field)
        field, q, w = fixture(T)
        p = field.potential
        fixture === edge_poisson_field && @test basis(p) isa RecombinedBSplineBasis
        @test update!(field, q, w, zero(T)) === nothing

        # the same deposit through the matrix `_deposit!` of `projection!`
        other = Potential(basis(p))
        projection!(other, particle_distribution(T, q, w))
        PoissonSolvers.update!(other)

        @test PoissonSolvers.rhs(p) == PoissonSolvers.rhs(other)
        @test coefficients(field) == coefficients(other)
        @test coefficients(field) === PoissonSolvers.coefficients(p)
    end

    @testset "the deposit is the weighted basis values" begin
        # The independent expectation of the deposit: `evaluate` one basis index at a time, which
        # is not the block evaluation the deposit loop reads. The weights are unequal and neither
        # is the particle count's inverse, so a deposit that dropped them, or took one for a
        # mean, fails here; an index that is off by one puts the value beside its own.
        for basis_constructor in (PeriodicBasisSpline, DirichletBasisSpline)
            b = basis_constructor((zero(T), one(T)), 3, 16)
            field = PoissonField(Potential(b))
            q = T[0.3, 0.72]
            w = T[2.5, T(1) / 3]

            update!(field, q, w, zero(T))

            expected = zeros(T, nbasis(b))
            for i in eachindex(q)
                expected .+= w[i] .* [evaluate(b, j, q[i]) for j in 1:nbasis(b)]
            end
            @test PoissonSolvers.rhs(field.potential) == expected
        end
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

        # `energy` is the diagnostic a snapshot march saves at every state, so it is inferred and
        # allocates nothing like the two step functions
        @test @inferred(energy(field)) isa T

        # coverage instrumentation allocates on every counted line, so under it the assertion
        # does not run rather than being weakened
        if Base.JLOptions().code_coverage == 0
            @test _alloc_update(update!, field, q, w, t) == 0
            @test _alloc_efield(efield!, field, e, q) == 0
            @test _alloc_update(field, e, q, w, t) == 0
            @test _alloc_energy(energy, field) == 0
        end
        # the controls: the barrier sees an allocation when there is one. `similar` rather than
        # `copy`: on Julia 1.11 a `copy` whose result is dead is elided and reads zero, while
        # `similar` always allocates.
        @test _alloc_update((f, x, w, t) -> similar(x), field, q, w, t) > 0
        @test _alloc_efield((f, e, x) -> similar(e), field, e, q) > 0
        @test _alloc_energy(f -> similar(f.nodes), field) > 0
    end
end
