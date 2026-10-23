using LinearAlgebra
using ParticleMethods: ParticleList
using PoissonSolvers
using SimpleSplines: SplineQuadrature, stiffness_matrix
using Test
using VlasovMethods
using VlasovMethods: ElectricField, PoissonField, coefficients, efield!, energy,
                     projection!,
                     update!, _with_scale
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

# The phrase of the exception a call raises, or a fixed word when it raises nothing: the guard
# of `ScaledField` is the only check on a `χ` of the wrong sign, so its phrase is the evidence.
function refusal(f)
    try
        f()
        return "no exception"
    catch thrown
        return sprint(showerror, thrown)
    end
end

@testset "check 1: ScaledField, $T" for T in (Float64, Float32)
    field, q, w = poisson_field(T)
    t = zero(T)

    # `χ = 1` is the inner field exactly
    one_ = ScaledField(field, one(T))
    @test one_ isa ElectricField
    @test one_.field === field
    e = similar(q)
    e′ = similar(q)
    update!(field, q, w, t)
    efield!(field, e′, q)
    update!(one_, q, w, t)
    @test efield!(one_, e, q) === e
    @test e == e′
    @test energy(one_) == energy(field)
    @test coefficients(one_) == coefficients(field)

    # `χ = 2` divides the field and the energy by four, exactly, and an integer `χ` is accepted
    two = ScaledField(field, 2)
    update!(two, q, w, t)
    e2 = similar(q)
    @test efield!(two, e2, q) === e2
    @test e2 == e′ ./ 4
    @test energy(two) == energy(field) / 4
    @test eltype(e2) === T

    # a `Float64` `χ` on a `Float32` field stays `Float32`
    if T === Float32
        half = ScaledField(field, 0.5)
        update!(half, q, w, t)
        @test eltype(efield!(half, similar(q), q)) === Float32
        @test eltype(energy(half)) === Float32
    end

    # `update!` forwards: the scaled field's coefficients are those of the same deposit into a
    # plain field of the same basis, and a drop of the forward leaves them at zero
    a = ScaledField(
        PoissonField(Potential(PeriodicBasisSpline((zero(T), one(T)), 3, 16))), 2)
    b = PoissonField(Potential(PeriodicBasisSpline((zero(T), one(T)), 3, 16)))
    update!(a, q, w, t)
    update!(b, q, w, t)
    @test coefficients(a) == coefficients(b)
    @test any(!iszero, coefficients(a))
    @test coefficients(ScaledField(field, 2)) == coefficients(field)

    # a `χ` that is neither finite nor positive is refused, and its phrase says which value
    for χ in (0, 0.0, -1, -1.0, T(NaN), T(Inf), T(-Inf))
        @test_throws ArgumentError ScaledField(field, χ)
        @test occursin("χ must be finite and positive",
            refusal(() -> ScaledField(field, χ)))
    end

    @testset "the field and the energy allocate nothing" begin
        sc = ScaledField(field, T(2))
        es = similar(q)
        @test @inferred(efield!(sc, es, q)) === es
        @test @inferred(energy(sc)) isa T

        # `update!` runs once per step of a march and `coefficients` once per save, so both are
        # inferred and allocate nothing on the scaled field as well
        @test @inferred(update!(sc, q, w, t)) === nothing
        @test @inferred(coefficients(sc)) isa typeof(coefficients(field))

        # coverage instrumentation allocates on every counted line, so under it the assertion
        # does not run rather than being weakened
        if Base.JLOptions().code_coverage == 0
            @test _alloc_efield(efield!, sc, es, q) == 0
            @test _alloc_energy(energy, sc) == 0
            # `_alloc_energy` is the one-argument barrier, which measures any one-argument read
            @test _alloc_update(update!, sc, q, w, t) == 0
            @test _alloc_energy(coefficients, sc) == 0
        end
        @test _alloc_efield((f, e, x) -> similar(e), sc, es, q) > 0
        @test _alloc_energy(f -> similar(f.field.nodes), sc) > 0
        @test _alloc_update((f, x, w, t) -> similar(x), sc, q, w, t) > 0
    end
end

@testset "_with_scale refuses a stored χ that the march would discard, $T" for T in (Float64, Float32)
    # `_with_scale` is the one place a march reads the scale of its field: a `ScaledField` is
    # rebuilt at the parameters' `χ`, and a stored `χ` that differs would leave the step scale and
    # the field scale out of step. The phrase is asserted because the constructor's own refusal of
    # a non-finite or non-positive `χ` raises the same exception type.
    field, q, _ = poisson_field(T)
    χ = T(2)

    @test _with_scale(field, χ) === field

    scaled = _with_scale(ScaledField(field, χ), χ)
    @test scaled isa ScaledField
    @test scaled.χ == χ
    @test efield!(scaled, similar(q), q) == efield!(field, similar(q), q) ./ 4

    for stored in (T(0.375), T(0.25))
        @test_throws ArgumentError _with_scale(ScaledField(field, stored), χ)
        msg = refusal(() -> _with_scale(ScaledField(field, stored), χ))
        @test occursin("a ScaledField at χ", msg)
        @test occursin("ScaledField at χ = $stored", msg)
        @test occursin("params.χ = $χ", msg)
    end

    # the comparison is exact: a stored `χ` one ulp away is another scale
    close = nextfloat(χ)
    @test_throws ArgumentError _with_scale(ScaledField(field, close), χ)
    @test occursin("a ScaledField at χ", refusal(() -> _with_scale(ScaledField(field, close), χ)))
end

@testset "the energy of a ScaledField carries χ², the potential of its field" begin
    # The scaled field's acceleration `E/χ²` is the force of the energy `W/χ²` per unit weight:
    # `wᵢ aᵢ = -∂W_χ/∂xᵢ`. That identity, not the square of the divided field, fixes the divisor,
    # so the energy divides by `χ²` and not by `χ⁴`. Here `χ = 3`, which no power of two hides.
    field, q, w = poisson_field(Float64)
    χ = 3.0
    sc = ScaledField(field, χ)

    function energy_at(x)
        update!(sc, x, w, 0.0)
        return energy(sc)
    end

    update!(sc, q, w, 0.0)
    force = w .* efield!(sc, similar(q), q)

    h = 1e-6
    gradient = map(eachindex(q)) do i
        qp, qm = copy(q), copy(q)
        qp[i] += h
        qm[i] -= h
        (energy_at(qp) - energy_at(qm)) / 2h
    end
    scale = maximum(abs, force)
    @test scale > 0
    @test maximum(abs, force .+ gradient) ≤ 1e-8 * scale

    # the controls: the energy of the divided field, `W/χ⁴`, and the unscaled energy are not the
    # potential of that force
    @test maximum(abs, force .+ gradient ./ χ^2) > 0.1 * scale
    @test maximum(abs, force .+ gradient .* χ^2) > 0.1 * scale
end
