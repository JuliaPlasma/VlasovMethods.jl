using LinearAlgebra
using PoissonSolvers
using Random
using ReducedBasisMethods: get_DEIM_indices, get_DEIM_interpolation_matrix
using Test
# `import` and not `using`: the `CountingField` stub below extends `update!` and `efield!`, which
# a `using` import does not allow.
import VlasovMethods: DEIMElectricField, ReducedElectricField, efield!, energy, update!
using VlasovMethods

include("../helpers/allocation_barriers.jl")
include("../helpers/poisson_field.jl")

# A seeded random orthonormal `n × k` basis, from the QR of a Gaussian matrix.
function orthogonal(::Type{T}, n, k; seed = 42) where {T}
    Random.seed!(seed)
    return Matrix(qr(randn(T, n, k)).Q)
end

# A full-order field that records how many positions each `efield!` is asked for: the
# hyper-reduction of the DEIM path is that the field is evaluated at the interpolation points
# rather than at every particle, which is what the timing script measures.
struct CountingField{FT} <: ElectricField
    field::FT
    counts::Vector{Int}
end

update!(f::CountingField, x, w, t) = update!(f.field, x, w, t)
function efield!(f::CountingField, e, x)
    push!(f.counts, length(x))
    return efield!(f.field, e, x)
end

# The buffers and the bases of a reduced field follow the element type of its input. A
# `Float64` lift under a `Float32` input shows up in `f.y` first, and a `Float64` projection in
# `f.Ψₑ`: `mul!` rounds either into a `Float32` buffer instead of failing, and every comparison
# of values would still pass.
function test_field_types(f, ::Type{T}) where {T}
    @test eltype(f.x) === T
    @test eltype(f.y) === T
    @test eltype(f.e) === T
    @test eltype(f.Ψₓ) === T
    @test eltype(f.Ψₚ) === T
    @test eltype(f.Ψₑ) === T
end

@testset "ReducedElectricField and DEIMElectricField, $T" for T in (Float64, Float32)
    # the tolerance of the reduced field against the full-order one: the `√eps(T)` of a
    # `Float64` run, and the `Float32` accuracy of the Poisson solve, `1e-5`
    local tol = T === Float64 ? sqrt(eps(T)) : T(1e-5)

    @testset "check 1: the three-argument constructor takes a Matrix basis" begin
        # the field `Ψₚ` is typed `ΨₚT`, the parameter the struct declares, so a `Matrix` basis
        # is accepted: a field typed `ΨₑT` rejects it, because the constructor passes `Ψ` as
        # `Ψₚ` and `Ψ'` as `Ψₑ`.
        field, q, _ = poisson_field(T; n = 20)
        Ψ = orthogonal(T, 20, 4)
        X = Ψ' * q

        rfield = ReducedElectricField(field, Ψ, X)

        # `Ψ` is the position basis and `Ψ'` the projection of the field onto it
        @test rfield isa ReducedElectricField
        @test rfield.Ψₓ === Ψ
        @test rfield.Ψₚ === Ψ
        @test rfield.Ψₑ == Ψ'
        @test rfield.field === field

        # the buffers follow the element type of the input
        test_field_types(rfield, T)
        @test size(rfield.x) == size(q)
        @test size(rfield.y) == size(q)
        @test size(rfield.e) == size(q)
    end

    @testset "check 2: the reduced field matches the full-order field" begin
        # at full basis size the two fields must agree exactly, which is the only setting where
        # they can. `Ψ` is a seeded random orthogonal matrix, not the identity: the identity
        # hides a transpose.
        field, q, w = poisson_field(T; n = 1000)
        Ψ = orthogonal(T, 1000, 1000)

        update!(field, q, w, zero(T))
        E_full = similar(q)
        efield!(field, E_full, q)

        # the reduced coordinates of the positions, lifted back by the basis
        X = Ψ' * q
        rfield = ReducedElectricField(field, Ψ, X)
        E_red = similar(X)
        rfield(E_red, X, w, zero(T))

        # the reduced field lifted back onto the particles equals the full-order field
        @test norm(Ψ * E_red .- E_full) / norm(E_full) ≤ tol
    end

    @testset "check 3: the DEIM selection" begin
        # the orthonormal input below: column 2 has its largest magnitude at a
        # negative entry, so a selection by `argmax` picks row 1 again and leaves `ΠₑᵀΨₑ`
        # singular. The selection by `max|r|` picks row 2.
        field, q, w = poisson_field(T; n = 3)
        Ψₑ = T[1 0; 0 -2/√5; 0 -1/√5]
        p = get_DEIM_indices(Ψₑ)

        kₑ = size(Ψₑ, 2)
        @test length(unique(p)) == kₑ

        Ψₓ = orthogonal(T, 3, 2; seed = 5)
        X = Ψₓ' * q
        dfield = DEIMElectricField(field, Ψₓ, Ψₑ, p, X)
        test_field_types(dfield, T)

        E = similar(X)
        @test dfield(E, X, w, zero(T)) === E
        @test all(isfinite, E)
    end

    @testset "check 4: the DEIM error bound" begin
        # the bound of the paper `deim.jl` cites, with the field sample `e` and `‖Ψₓᵀ‖₂ = 1`
        field, q, w = poisson_field(T; n = 1000)
        k = 20
        Ψₓ = orthogonal(T, 1000, k)
        @test norm(Ψₓ' * Ψₓ - Matrix{T}(I, k, k)) ≤ sqrt(eps(T))

        X = Ψₓ' * q
        rfield = ReducedElectricField(field, Ψₓ, X)
        E_red = similar(X)
        rfield(E_red, X, w, zero(T))

        # at `m = N` the field basis is square and the DEIM interpolation is a permutation, so
        # the two paths agree to round-off
        Ψₑ = orthogonal(T, 1000, 1000; seed = 7)
        dfield = DEIMElectricField(field, Ψₓ, Ψₑ, get_DEIM_indices(Ψₑ), X)
        test_field_types(dfield, T)
        E_deim = similar(X)
        dfield(E_deim, X, w, zero(T))
        @test norm(E_deim .- E_red) / norm(E_red) ≤ tol

        # at `m = 1` the DEIM approximation is at its coarsest, and the measured error is
        # bounded by the paper's bound
        φ = orthogonal(T, 1000, 1; seed = 11)
        p = get_DEIM_indices(φ)
        dfield = DEIMElectricField(field, Ψₓ, φ, p, X)
        test_field_types(dfield, T)
        dfield(E_deim, X, w, zero(T))

        # the field sample: the field at every lifted position, on the deposit both paths made
        x = Ψₓ * X
        e = similar(x)
        efield!(field, e, x)

        # the rank-one path is the interpolation formula itself: the field sample at the lifted
        # positions, reduced through the one interpolation point `p` and lifted back by `Ψₓ`.
        # The expectation is computed here from the sample, not read back from the field, so a
        # projection that returns zero fails it. The sample is asserted nonzero first, because
        # a zero field would make every comparison vacuous.
        @test norm(e) > zero(T)
        expected = Ψₓ' * φ * (φ[p, :] \ e[p])
        @test norm(expected) > zero(T)
        @test norm(E_deim .- expected) ≤ 10 * sqrt(eps(T)) * norm(expected)

        bound = opnorm(inv(φ[p, :])) * norm(e - φ * (φ' * e)) + 10 * eps(T) * norm(e)
        @test norm(E_deim .- E_red) ≤ bound
    end

    @testset "check 5: the DEIM field is evaluated at the interpolation points" begin
        # the field of the DEIM path is evaluated at the `kₑ` interpolation points, not at every
        # particle: `scripts/electric_field_timing.jl` times what this counts
        field, q, _ = poisson_field(T; n = 60)
        n, k, kₑ = 60, 4, 4
        Ψₓ = orthogonal(T, n, k; seed = 1)
        Ψₑ = orthogonal(T, n, kₑ; seed = 2)
        p = get_DEIM_indices(Ψₑ)
        X = Ψₓ' * q

        counting = CountingField(field, Int[])

        rfield = ReducedElectricField(counting, Ψₓ, X)
        E = similar(X)
        efield!(rfield, E, X)
        @test counting.counts == [n]

        empty!(counting.counts)
        dfield = DEIMElectricField(counting, Ψₓ, Ψₑ, p, X)
        test_field_types(dfield, T)
        efield!(dfield, E, X)
        @test counting.counts == [kₑ]
    end

    @testset "check 6: the index and the matrix forms agree bit for bit" begin
        # the matrix form reads its row indices from `Πₑ` and calls the index form, so the two
        # fields are the same: `isequal` compares the bits, including the sign of a zero
        field, q, w = poisson_field(T; n = 200)
        n, k, kₑ = 200, 20, 8
        Ψₓ = orthogonal(T, n, k; seed = 3)
        Ψₑ = orthogonal(T, n, kₑ; seed = 4)
        X = Ψₓ' * q

        ifield = DEIMElectricField(field, Ψₓ, Ψₑ, get_DEIM_indices(Ψₑ), X)
        mfield = DEIMElectricField(field, Ψₓ, Ψₑ, get_DEIM_interpolation_matrix(Ψₑ), X)
        test_field_types(ifield, T)
        test_field_types(mfield, T)

        update!(ifield, X, w, zero(T))
        E_index = efield!(ifield, similar(X), X)
        update!(mfield, X, w, zero(T))
        E_matrix = efield!(mfield, similar(X), X)

        @test norm(E_index) > zero(T)
        @test isequal(E_index, E_matrix)
    end

    @testset "the diagnostics are the full-order field's" begin
        field, q, w = poisson_field(T; n = 20)
        Ψ = orthogonal(T, 20, 4)
        X = Ψ' * q

        update!(field, q, w, zero(T))
        rfield = ReducedElectricField(field, Ψ, X)
        @test energy(rfield) == energy(field)
        @test coefficients(rfield) === coefficients(field)
    end

    @testset "update! and efield! are inferred and allocate nothing" begin
        field, q, w = poisson_field(T; n = 200)
        n, k, kₑ = 200, 20, 8
        Ψₓ = orthogonal(T, n, k; seed = 3)
        Ψₑ = orthogonal(T, n, kₑ; seed = 4)
        X = Ψₓ' * q

        rfield = ReducedElectricField(field, Ψₓ, X)
        dfield = DEIMElectricField(field, Ψₓ, Ψₑ, get_DEIM_indices(Ψₑ), X)
        test_field_types(rfield, T)
        test_field_types(dfield, T)
        E = similar(X)
        t = zero(T)

        for f in (rfield, dfield)
            @test @inferred(update!(f, X, w, t)) === nothing
            @test @inferred(efield!(f, E, X)) === E
            @test @inferred(f(E, X, w, t)) === E
        end

        # coverage instrumentation allocates on every counted line, so under it the assertion
        # does not run rather than being weakened
        if Base.JLOptions().code_coverage == 0
            for f in (rfield, dfield)
                @test _alloc_update(update!, f, X, w, t) == 0
                @test _alloc_efield(efield!, f, E, X) == 0
            end
        end
        # the controls: the barrier sees an allocation when there is one. `similar` rather than
        # `copy`: on Julia 1.11 a `copy` whose result is dead is elided and reads zero, while
        # `similar` always allocates.
        @test _alloc_update((f, x, w, t) -> similar(x), rfield, X, w, t) > 0
        @test _alloc_efield((f, e, x) -> similar(e), rfield, E, X) > 0
    end
end
