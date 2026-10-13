using GeometricBrackets: Arakawa
using LinearAlgebra: qr
using MultiIndexArrays: GridTensor, MultiIndexAxis, MultiIndexLazyArray
using Random
using Test
using VlasovMethods

# Seeded orthonormal projections, converted to `DT`.
function orthonormal(rng::AbstractRNG, ::Type{DT}, n::Int, m::Int) where {DT}
    return DT.(Matrix(qr(randn(rng, n, m)).Q))
end

# The 1D1V bracket case on the `nx × nv` grid: the `x`-grid is periodic and the `v`-grid bounded,
# so neither spacing is a power of two, and `v²` over the asymmetric velocity range is not
# symmetric in `n₂`. `Pi` and `Pj` are drawn with different seeds and have different numbers of
# columns, so a swap of the two index slots cannot pass a comparison against the reference.
function fixture(::Type{DT}, nx::Int, nv::Int, mi = 3, mj = 2, mk = 4) where {DT}
    local dist = GridDistribution(zeros(DT, nx, nv), (0, 2π), (-2.3, 3.1))
    local N = nx * nv
    return (; nx, nv, N, gt = bracket_tensor(Arakawa(nx, nv, step(dist.x), step(dist.v))),
        v² = collect(dist.v .^ 2), Pi = orthonormal(Xoshiro(1), DT, N, mi),
        Pj = orthonormal(Xoshiro(2), DT, N, mj), Pk = orthonormal(Xoshiro(3), DT, N, mk),
        Px = orthonormal(Xoshiro(4), DT, nx, mk))
end

# The width-`w` periodic neighbourhood test of `J` around `I`, written from the wrapped offsets by
# hand and never through `stencil_indices`, so that a defect in the stencil cannot reach the
# reference and the code under test together.
function in_neighbourhood(I::CartesianIndex, J::CartesianIndex, w::Int, g::Tuple)
    all(mod.(Tuple(J - I) .+ w, g) .<= 2w)
end

# The width-`w` `GridTensor` whose nonzero entries lie in the width-`w` neighbourhood of the first
# index and nowhere else: the contract of the locality interface, on a table that is not Arakawa's.
function neighbourhood_tensor(rng::AbstractRNG, ::Type{DT}, g::Tuple, w::Int) where {DT}
    local n = prod(g)
    local L = LinearIndices(g)
    local C = rand(rng, DT, n, n, n) .- DT(1/2)
    local ax = MultiIndexAxis(g...)
    return GridTensor(
        MultiIndexLazyArray(DT,
            (I, J, K) -> in_neighbourhood(I, J, w, g) && in_neighbourhood(I, K, w, g) ?
                         C[L[I], L[J], L[K]] : zero(DT), ax, ax, ax),
        w)
end

# The hand projection of the materialised tensor: each slice of the third index is contracted with
# `Πiᵀ` and `Πj`. It calls no `ReducedTensor`, no `stencil_indices` and no `axes(parent(…))`, so
# a wrong derivation of a grid index in the code cannot be matched by the same error here.
function hand_projection(tensor, Pi, Pj)
    local A = Array(tensor)
    local Rd = zeros(eltype(A), size(Pi, 2), size(Pj, 2), size(A, 3))
    for K in axes(A, 3)
        Rd[:, :, K] .= Pi' * A[:, :, K] * Pj
    end
    return Rd
end

# The third index of `Rd` contracted with the columns of `W`, or with the vector `W`.
function contract(Rd, W::AbstractVecOrMat)
    reshape(reshape(Rd, :, size(Rd, 3)) * W, size(Rd, 1), size(Rd, 2), size(W)[2:end]...)
end

# `√eps(T)` of the largest reference entry, after asserting that the reference is not zero. The
# tolerance is the round-off of the two sums, and the entries of a projection are far below `1`
# on these grids, so the bound is relative to the reference rather than to one.
function matches(x, ref)
    maximum(abs, ref) > 0 &&
        maximum(abs, Array(x) .- ref) ≤ sqrt(eps(eltype(ref))) * maximum(abs, ref)
end

# The allocations of the second call, as in `collisions.jl`: the arguments are concrete, and the
# `where` clause forces specialisation on the function argument.
_alloc_call2(f::F, a::A, b::B) where {F, A, B} = (f(a, b); @allocated f(a, b))

@testset "Reduced grid tensors" begin
    @testset "bracket_tensor is the bracket on its own grid" begin
        # a non-square grid, so that the extents swapped or one extent used twice fail; the
        # reference reads the `Arakawa` on its own nodes, not through a `MultiIndexAxis`
        local a = Arakawa(5, 4, 0.5, 0.5)
        local gt = bracket_tensor(a)
        local C = vec(CartesianIndices((5, 4)))
        @test gt isa GridTensor{Float64, 3}
        @test size(gt) == (20, 20, 20)
        @test size(axes(parent(gt), 1)) == (5, 4)
        @test Array(gt) == [a(I, J, K) for I in C, J in C, K in C]
        @test eltype(bracket_tensor(Arakawa(5, 4, 0.5f0, 0.5f0))) == Float32
    end

    @testset "the three are the hand projection" begin
        for DT in (Float32, Float64), g in ((8, 8), (16, 16), (8, 6))

            local (; nx, nv, gt, v², Pi, Pj, Pk, Px) = fixture(DT, g...)
            local Rd = hand_projection(gt, Pi, Pj)

            local prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
            @test size(prt) == (3, 2, 4)
            @test matches(prt, contract(Rd, Pk))

            # `v²` has one entry per velocity node and is halved, so a `v⁴` over the same grid, a
            # `v²` of `N` entries and a missing `½` each fail here
            local vrm = VelocityReducedMatrix(gt, Pi, Pj, v²)
            @test size(vrm) == (3, 2)
            @test matches(vrm, contract(Rd, repeat(v² ./ 2, inner = nx)))

            # `Px` has one row per `x` node and is read at the grid's `x` index
            local frt = FullyReducedTensor(gt, Pi, Pj, Px)
            @test size(frt) == (3, 2, 4)
            @test matches(frt, contract(Rd, repeat(Px, nv, 1)))
        end
    end

    @testset "the contract tensor is read through the tensor's own stencil" begin
        for DT in (Float32, Float64), w in (1, 2)

            local (; nx, nv, v², Pi, Pj, Pk, Px) = fixture(DT, 8, 6)
            local gt = neighbourhood_tensor(Xoshiro(1), DT, (nx, nv), w)
            local Rd = hand_projection(gt, Pi, Pj)

            # the width is the tensor's, and the summation is the nested stencil the contract
            # allows: a width hard-coded to 1, or `m, n ∈ S(K)`, fails one of these
            @test matches(PotentialReducedTensor(gt, Pi, Pj, Pk), contract(Rd, Pk))
            @test matches(VelocityReducedMatrix(gt, Pi, Pj, v²),
                contract(Rd, repeat(v² ./ 2, inner = nx)))
            @test matches(FullyReducedTensor(gt, Pi, Pj, Px), contract(Rd, repeat(Px, nv, 1)))
        end
    end

    @testset "a grid of one dimension is only a potential" begin
        for DT in (Float32, Float64)
            # a width-1 tensor on a one-dimensional grid still satisfies the contract, and the
            # third index of a `PotentialReducedTensor` is one per grid node
            local gt = neighbourhood_tensor(Xoshiro(1), DT, (7,), 1)
            local Pi, Pj, Pk = (orthonormal(Xoshiro(s), DT, 7, m)
            for (s, m) in ((1, 3), (2, 2), (3, 4)))
            local prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
            @test size(prt) == (3, 2, 4)
            @test matches(prt, contract(hand_projection(gt, Pi, Pj), Pk))

            # a velocity index and an `x` index are one per node of a two-dimensional grid
            @test_throws ArgumentError VelocityReducedMatrix(gt, Pi, Pj, ones(DT, 7))
            @test_throws ArgumentError FullyReducedTensor(gt, Pi, Pj, Pk)
        end
    end

    @testset "the constructors refuse what they cannot project" begin
        for DT in (Float32, Float64)
            local OT = DT == Float32 ? Float64 : Float32
            local (; nx, nv, N, gt, v², Pi, Pj, Pk, Px) = fixture(DT, 8, 6)

            # the valid input, on the non-square grid that shows a guard against the other extent
            @test PotentialReducedTensor(gt, Pi, Pj, Pk) isa PotentialReducedTensor{DT}
            @test VelocityReducedMatrix(gt, Pi, Pj, v²) isa SubArray{DT, 2}
            @test FullyReducedTensor(gt, Pi, Pj, Px) isa PotentialReducedTensor{DT}

            # a guard against the wrong extent refuses the valid input, and a deleted guard lets
            # each of these construct
            @test_throws DimensionMismatch PotentialReducedTensor(
                gt, Pi, Pj, randn(Xoshiro(5), DT, N - 1, 4))
            @test_throws DimensionMismatch VelocityReducedMatrix(
                gt, Pi, Pj, randn(Xoshiro(5), DT, nx))
            @test_throws DimensionMismatch FullyReducedTensor(
                gt, Pi, Pj, randn(Xoshiro(5), DT, nv, 4))

            # the element type of every projection is the tensor's, not a promoted one
            @test_throws ArgumentError PotentialReducedTensor(gt, Pi, Pj, OT.(Pk))
            @test_throws ArgumentError VelocityReducedMatrix(gt, Pi, Pj, OT.(v²))
            @test_throws ArgumentError FullyReducedTensor(gt, Pi, Pj, OT.(Px))

            # the checks of `ReducedTensor` are the ones that run, and the constructors keep no
            # copy
            local long = orthonormal(Xoshiro(6), DT, N + 1, 3)
            @test_throws DimensionMismatch PotentialReducedTensor(gt, long, Pj, Pk)
            @test_throws DimensionMismatch VelocityReducedMatrix(gt, long, Pj, v²)
            @test_throws DimensionMismatch FullyReducedTensor(gt, long, Pj, Px)
            @test_throws ArgumentError PotentialReducedTensor(gt, OT.(Pi), Pj, Pk)
            @test_throws ArgumentError VelocityReducedMatrix(gt, OT.(Pi), Pj, v²)
            @test_throws ArgumentError FullyReducedTensor(gt, OT.(Pi), Pj, Px)

            # a direction shorter than `2w + 1` makes the stencil name a node twice, which
            # `ReducedTensor` refuses; the 2 × 5 grid is too short for width 1 in its first direction
            local short = neighbourhood_tensor(Xoshiro(7), DT, (2, 5), 1)
            local Qi, Qj, Qk = (orthonormal(Xoshiro(s), DT, 10, m)
            for (s, m) in ((1, 3), (2, 2), (3, 4)))
            @test_throws ArgumentError PotentialReducedTensor(short, Qi, Qj, Qk)
            @test_throws ArgumentError VelocityReducedMatrix(short, Qi, Qj, randn(Xoshiro(8), DT, 5))
            @test_throws ArgumentError FullyReducedTensor(short, Qi, Qj, Qk[1:2, :])
        end
    end

    @testset "the entries, their bounds, their type and their index forms" begin
        for DT in (Float32, Float64)
            local (; gt, v², Pi, Pj, Pk, Px) = fixture(DT, 8, 6)
            local prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
            local vrm = VelocityReducedMatrix(gt, Pi, Pj, v²)
            local frt = FullyReducedTensor(gt, Pi, Pj, Px)

            # a bare `@test_throws BoundsError` cannot fail: `ReducedTensor` and the projections
            # throw one as well, so the error has to name the array that was asked for
            for (x, I) in ((prt, (0, 1, 1)), (prt, (4, 1, 1)), (prt, (1, 1, 5)), (
                frt, (0, 1, 1)),
                (frt, (4, 1, 1)), (frt, (1, 1, 5)), (vrm, (0, 1)), (vrm, (1, 3)), (
                vrm, (1, 0)),
                (vrm, (4, 1)))
                local err = try
                    x[I...]
                    nothing
                catch e
                    e
                end
                @test err isa BoundsError && err.a === x
            end

            # the accumulator and the `½` stay in `DT`, which `eltype(Array(vrm))` cannot see:
            # `setindex!` converts a promoted entry back on the way in
            @test @inferred(prt[2, 1, 3]) isa DT
            @test @inferred(vrm[2, 1]) isa DT
            @test @inferred(frt[2, 1, 3]) isa DT

            # another integer type and a linear index go through `AbstractArray`'s fallbacks
            @test prt[Int32(2), Int32(1), Int32(3)] == prt[2, 1, 3]
            @test vrm[Int32(2), Int32(1)] == vrm[2, 1]
            @test frt[Int32(2), Int32(1), Int32(3)] == frt[2, 1, 3]
            @test prt[7] == Array(prt)[7]
            @test vrm[5] == Array(vrm)[5]
            @test frt[7] == Array(frt)[7]

            # past the rank, `size` and `axes` are `AbstractArray`'s and not the type's own
            @test size(prt, 4) == 1
            @test axes(prt, 4) == Base.OneTo(1)
            @test size(vrm, 3) == 1
            @test axes(vrm, 3) == Base.OneTo(1)
            @test size(frt, 4) == 1
            @test axes(frt, 4) == Base.OneTo(1)
        end
    end

    @testset "the entries allocate nothing" begin
        # `getindex` is what a reduced right-hand side evaluates, and it is called once per node
        # of the output
        for DT in (Float32, Float64), g in ((8, 8), (16, 16)), m in (2, 4)
            local (; gt, v², Pi, Pj, Pk, Px) = fixture(DT, g..., m, m, m)
            for x in (PotentialReducedTensor(gt, Pi, Pj, Pk),
                VelocityReducedMatrix(gt, Pi, Pj, v²),
                FullyReducedTensor(gt, Pi, Pj, Px))
                local out = Array{DT}(undef, size(x))
                local I = ntuple(_ -> 2, ndims(x))
                # coverage instrumentation allocates on every counted line, so under it the
                # assertions do not run rather than being weakened
                if Base.JLOptions().code_coverage == 0
                    @test _alloc_call2((r, o) -> copyto!(o, r), x, out) == 0
                    @test _alloc_call2((r, i) -> r[i...], x, I) == 0
                end
                # the controls: the barrier sees an allocation when there is one
                @test _alloc_call2((r, o) -> collect(r), x, out) > 0
                @test _alloc_call2((r, i) -> collect(r), x, I) > 0
            end
        end
    end
end
