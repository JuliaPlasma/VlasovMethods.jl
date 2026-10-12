using GeometricBrackets: Arakawa
using LinearAlgebra: Matrix, qr
using MultiIndexArrays: GridTensor, MultiIndexAxis, MultiIndexLazyArray
using Random
using Test
using VlasovMethods

# The 1D1V grid of the bracket case: the `x`-grid is periodic and the `v`-grid bounded, so neither
# spacing is a power of two, and `v²` over the asymmetric velocity range is not symmetric in `n₂`.
function grid(::Type{DT}, nx::Int, nv::Int) where {DT}
    GridDistribution(zeros(DT, nx, nv), (0, 2π), (-2.3, 3.1))
end

# The `GridTensor` of the Arakawa bracket on the `nx × nv` grid. `MultiIndexLazyArray` refuses a
# callable struct, so the coefficient function is a closure over the `Arakawa`.
function arakawa_tensor(::Type{DT}, nx::Int, nv::Int, hx, hv) where {DT}
    local a = Arakawa(nx, nv, hx, hv)
    local ax = MultiIndexAxis(1:nx, 1:nv)
    return GridTensor(MultiIndexLazyArray(DT, (I, J, K) -> a(I, J, K), ax, ax, ax), 1)
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

# Seeded orthonormal projections, converted to `DT`. `Pi` and `Pj` are drawn with different seeds,
# so a swap of the two index slots cannot pass a comparison against the reference.
function orthonormal(rng::AbstractRNG, ::Type{DT}, n::Int, m::Int) where {DT}
    return DT.(Matrix(qr(randn(rng, n, m)).Q)[:, 1:m])
end

# The hand projection of the materialised tensor: each slice of the third index is contracted with
# `Πiᵀ` and `Πj`. It calls no `ReducedTensor`, no `stencil_indices` and no `axes(parent(…))`, so a
# wrong derivation of a grid index in the code cannot be matched by the same error here.
function hand_projection(A::AbstractArray{DT, 3}, Pi, Pj) where {DT}
    local Rd = zeros(DT, size(Pi, 2), size(Pj, 2), size(A, 3))
    for K in axes(A, 3)
        Rd[:, :, K] .= Pi' * A[:, :, K] * Pj
    end
    return Rd
end

function contract(Rd::AbstractArray{DT, 3}, W::AbstractMatrix) where {DT}
    reshape(reshape(Rd, :, size(Rd, 3)) * W, size(Rd, 1), size(Rd, 2), size(W, 2))
end

function contract(Rd::AbstractArray{DT, 3}, w::AbstractVector) where {DT}
    reshape(reshape(Rd, :, size(Rd, 3)) * w, size(Rd, 1), size(Rd, 2))
end

# `√eps(T)` of the largest reference entry. The tolerance is the round-off of the two sums, and
# the entries of a projection are far below `1` on these grids, so the bound is relative to the
# reference rather than to one.
bound(::Type{DT}, R) where {DT} = sqrt(eps(DT)) * maximum(abs, R)

@testset "Reduced grid tensors" begin
    @testset "the three types are the hand projection" begin
        for DT in (Float32, Float64), g in ((8, 8), (16, 16), (8, 6))

            local nx, nv = g
            local N = nx * nv
            local dist = grid(DT, nx, nv)
            local gt = arakawa_tensor(DT, nx, nv, step(dist.x), step(dist.v))
            local v² = collect(dist.v .^ 2)

            # 3, 2 and 4 columns keep `m ≤ 4` and make `Πᵢ ≠ Πⱼ`
            local Pi = orthonormal(Xoshiro(1), DT, N, 3)
            local Pj = orthonormal(Xoshiro(2), DT, N, 2)
            local Pk = orthonormal(Xoshiro(3), DT, N, 4)
            local Px = orthonormal(Xoshiro(4), DT, nx, 4)
            @test Pi != Pj

            local A = Array(gt)
            local Rd = hand_projection(A, Pi, Pj)

            local prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
            local refp = contract(Rd, Pk)
            @test size(prt) == (3, 2, 4)
            @test maximum(abs, refp) > 0
            @test maximum(abs, Array(prt) .- refp) ≤ bound(DT, refp)

            # `v²` has one entry per velocity node and the `½` is inside, so a `v⁴` over the same
            # grid, a `v²` of `N` entries and a missing `½` each fail here
            local vrm = VelocityReducedMatrix(gt, Pi, Pj, v²)
            local refv = contract(Rd, repeat(v² ./ 2, inner = nx))
            @test size(vrm) == (3, 2)
            @test maximum(abs, refv) > 0
            @test maximum(abs, Array(vrm) .- refv) ≤ bound(DT, refv)

            # `Px` has one row per `x` node and is read at the grid's `x` index
            local frt = FullyReducedTensor(gt, Pi, Pj, Px)
            local reff = contract(Rd, repeat(Px, nv, 1))
            @test size(frt) == (3, 2, 4)
            @test maximum(abs, reff) > 0
            @test maximum(abs, Array(frt) .- reff) ≤ bound(DT, reff)
        end
    end

    @testset "the contract tensor is read through the tensor's own stencil" begin
        for DT in (Float32, Float64), w in (1, 2)

            local nx, nv = 8, 6
            local N = nx * nv
            local dist = grid(DT, nx, nv)
            local v² = collect(dist.v .^ 2)
            local gt = neighbourhood_tensor(Xoshiro(1), DT, (nx, nv), w)
            local Pi = orthonormal(Xoshiro(1), DT, N, 3)
            local Pj = orthonormal(Xoshiro(2), DT, N, 2)
            local Pk = orthonormal(Xoshiro(3), DT, N, 4)
            local Px = orthonormal(Xoshiro(4), DT, nx, 4)
            local Rd = hand_projection(Array(gt), Pi, Pj)

            local refp = contract(Rd, Pk)
            local refv = contract(Rd, repeat(v² ./ 2, inner = nx))
            local reff = contract(Rd, repeat(Px, nv, 1))
            @test maximum(abs, refp) > 0
            @test maximum(abs, refv) > 0
            @test maximum(abs, reff) > 0

            # the width is the tensor's, and the summation is the nested stencil the contract
            # allows: a width hard-coded to 1, or `m, n ∈ S(K)`, fails one of these
            @test maximum(abs, Array(PotentialReducedTensor(gt, Pi, Pj, Pk)) .- refp) ≤
                  bound(DT, refp)
            @test maximum(abs, Array(VelocityReducedMatrix(gt, Pi, Pj, v²)) .- refv) ≤
                  bound(DT, refv)
            @test maximum(abs, Array(FullyReducedTensor(gt, Pi, Pj, Px)) .- reff) ≤
                  bound(DT, reff)
        end
    end

    @testset "a grid of one dimension is only a potential" begin
        local DT = Float64
        # a width-1 tensor on a one-dimensional grid still satisfies the contract, so the dense
        # reference is what the stencil sum has to reproduce
        local gt = neighbourhood_tensor(Xoshiro(1), DT, (7,), 1)
        local Pi = orthonormal(Xoshiro(1), DT, 7, 3)
        local Pj = orthonormal(Xoshiro(2), DT, 7, 2)
        local Pk = orthonormal(Xoshiro(3), DT, 7, 4)
        local Px = orthonormal(Xoshiro(4), DT, 7, 4)

        # the third index of a `PotentialReducedTensor` is one per grid node, so a
        # one-dimensional grid serves it
        local prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
        local ref = contract(hand_projection(Array(gt), Pi, Pj), Pk)
        @test size(prt) == (3, 2, 4)
        @test maximum(abs, ref) > 0
        @test maximum(abs, Array(prt) .- ref) ≤ bound(DT, ref)

        # a velocity index and an `x` index are one per node of a two-dimensional grid
        @test_throws ArgumentError VelocityReducedMatrix(gt, Pi, Pj, ones(DT, 7))
        @test_throws ArgumentError FullyReducedTensor(gt, Pi, Pj, Px)
    end

    @testset "the constructor refuses what it cannot project" begin
        local DT = Float32
        local nx, nv = 8, 6
        local N = nx * nv
        local dist = grid(DT, nx, nv)
        local gt = arakawa_tensor(DT, nx, nv, step(dist.x), step(dist.v))
        local v² = collect(dist.v .^ 2)
        local Pi = orthonormal(Xoshiro(1), DT, N, 3)
        local Pj = orthonormal(Xoshiro(2), DT, N, 2)
        local Pk = orthonormal(Xoshiro(3), DT, N, 4)
        local Px = orthonormal(Xoshiro(4), DT, nx, 4)

        # the valid input, on the non-square grid that shows a guard against the other extent
        @test PotentialReducedTensor(gt, Pi, Pj, Pk) isa PotentialReducedTensor
        @test VelocityReducedMatrix(gt, Pi, Pj, v²) isa VelocityReducedMatrix
        @test FullyReducedTensor(gt, Pi, Pj, Px) isa FullyReducedTensor

        # a guard against the wrong extent refuses the valid input, and a deleted guard lets each
        # of these construct
        @test_throws DimensionMismatch PotentialReducedTensor(
            gt, Pi, Pj, randn(Xoshiro(5), DT, N - 1, 4))
        @test_throws DimensionMismatch VelocityReducedMatrix(
            gt, Pi, Pj, randn(Xoshiro(5), DT, nx))
        @test_throws DimensionMismatch FullyReducedTensor(
            gt, Pi, Pj, randn(Xoshiro(5), DT, nv, 4))

        # the element type of every projection is the tensor's, not a promoted one
        @test_throws ArgumentError PotentialReducedTensor(gt, Pi, Pj, Float64.(Pk))
        @test_throws ArgumentError VelocityReducedMatrix(gt, Pi, Pj, Float64.(v²))
        @test_throws ArgumentError FullyReducedTensor(gt, Pi, Pj, Float64.(Px))

        # the checks of `ReducedTensor` are the ones that run, and the constructors keep no copy
        local long = orthonormal(Xoshiro(6), DT, N + 1, 3)
        @test_throws DimensionMismatch PotentialReducedTensor(gt, long, Pj, Pk)
        @test_throws DimensionMismatch VelocityReducedMatrix(gt, long, Pj, v²)
        @test_throws DimensionMismatch FullyReducedTensor(gt, long, Pj, Px)
        @test_throws ArgumentError PotentialReducedTensor(gt, Float64.(Pi), Pj, Pk)
        @test_throws ArgumentError VelocityReducedMatrix(gt, Float64.(Pi), Pj, v²)
        @test_throws ArgumentError FullyReducedTensor(gt, Float64.(Pi), Pj, Px)

        # a direction shorter than `2w + 1` makes the stencil name a node twice, which
        # `ReducedTensor` refuses; the 2 × 5 grid is long enough for width 1 per direction
        local short = neighbourhood_tensor(Xoshiro(7), DT, (2, 5), 1)
        local Qi = orthonormal(Xoshiro(1), DT, 10, 3)
        local Qj = orthonormal(Xoshiro(2), DT, 10, 2)
        local Qk = orthonormal(Xoshiro(3), DT, 10, 4)
        local Qx = orthonormal(Xoshiro(4), DT, 2, 2)
        local w² = randn(Xoshiro(8), DT, 5)
        @test_throws ArgumentError PotentialReducedTensor(short, Qi, Qj, Qk)
        @test_throws ArgumentError VelocityReducedMatrix(short, Qi, Qj, w²)
        @test_throws ArgumentError FullyReducedTensor(short, Qi, Qj, Qx)
    end

    @testset "an index off the array throws its own BoundsError" begin
        for DT in (Float32, Float64)
            local nx, nv = 8, 6
            local N = nx * nv
            local dist = grid(DT, nx, nv)
            local gt = arakawa_tensor(DT, nx, nv, step(dist.x), step(dist.v))
            local v² = collect(dist.v .^ 2)
            local Pi = orthonormal(Xoshiro(1), DT, N, 3)
            local Pj = orthonormal(Xoshiro(2), DT, N, 2)
            local Pk = orthonormal(Xoshiro(3), DT, N, 4)
            local Px = orthonormal(Xoshiro(4), DT, nx, 4)
            local prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
            local vrm = VelocityReducedMatrix(gt, Pi, Pj, v²)
            local frt = FullyReducedTensor(gt, Pi, Pj, Px)

            # a bare `@test_throws BoundsError` cannot fail: `ReducedTensor` and the projections
            # throw one as well, so the error has to name the array that was asked for
            for (x, i, j, k) in ((prt, 0, 1, 1), (prt, 4, 1, 1), (prt, 1, 1, 5),
                (frt, 0, 1, 1), (frt, 4, 1, 1), (frt, 1, 1, 5))
                local err = try
                    x[i, j, k]
                    nothing
                catch e
                    e
                end
                @test err isa BoundsError && err.a === x
            end
            for (x, i, j) in ((vrm, 0, 1), (vrm, 1, 3), (vrm, 1, 0), (vrm, 4, 1))
                local err = try
                    x[i, j]
                    nothing
                catch e
                    e
                end
                @test err isa BoundsError && err.a === x
            end
        end
    end

    @testset "the entries, their type and their index forms" begin
        for DT in (Float32, Float64)
            local nx, nv = 8, 6
            local N = nx * nv
            local dist = grid(DT, nx, nv)
            local gt = arakawa_tensor(DT, nx, nv, step(dist.x), step(dist.v))
            local v² = collect(dist.v .^ 2)
            local Pi = orthonormal(Xoshiro(1), DT, N, 3)
            local Pj = orthonormal(Xoshiro(2), DT, N, 2)
            local Pk = orthonormal(Xoshiro(3), DT, N, 4)
            local Px = orthonormal(Xoshiro(4), DT, nx, 4)
            local prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
            local vrm = VelocityReducedMatrix(gt, Pi, Pj, v²)
            local frt = FullyReducedTensor(gt, Pi, Pj, Px)

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
        # of the output. A fresh process with `--check-bounds=auto`, because `julia_cmd()` carries
        # the parent's bounds setting, which the later flag overrides.
        local prog = """
            # A test runner can pass on a `JULIA_LOAD_PATH` without the active project, and
            # `--project` alone then does not make this package loadable.
            pushfirst!(LOAD_PATH, "@")
            using VlasovMethods
            using GeometricBrackets: Arakawa
            using LinearAlgebra: Matrix, qr
            using MultiIndexArrays: GridTensor, MultiIndexAxis, MultiIndexLazyArray
            using Random

            # The allocations of the second call, behind a barrier that specialises on every
            # argument, so that only the call itself is measured.
            allocations(f::F, a::A, b::B) where {F, A, B} = (f(a, b); @allocated f(a, b))
            entry(f::F, a::A, i::I) where {F, A, I} = (f(a, i); @allocated f(a, i))

            Base.JLOptions().check_bounds == 0 || error("this program measures allocations")

            for DT in (Float32, Float64), g in ((8, 8), (16, 16)), m in (2, 4)
                nx, nv = g
                N = nx * nv
                dist = GridDistribution(zeros(DT, nx, nv), (0, 2π), (-2.3, 3.1))
                a = Arakawa(nx, nv, step(dist.x), step(dist.v))
                ax = MultiIndexAxis(1:nx, 1:nv)
                gt = GridTensor(MultiIndexLazyArray(DT, (I, J, K) -> a(I, J, K), ax, ax, ax), 1)
                v² = collect(dist.v .^ 2)
                Pi = DT.(Matrix(qr(randn(Xoshiro(1), N, m)).Q)[:, 1:m])
                Pj = DT.(Matrix(qr(randn(Xoshiro(2), N, m)).Q)[:, 1:m])
                Pk = DT.(Matrix(qr(randn(Xoshiro(3), N, m)).Q)[:, 1:m])
                Px = DT.(Matrix(qr(randn(Xoshiro(4), nx, m)).Q)[:, 1:m])

                prt = PotentialReducedTensor(gt, Pi, Pj, Pk)
                vrm = VelocityReducedMatrix(gt, Pi, Pj, v²)
                frt = FullyReducedTensor(gt, Pi, Pj, Px)

                for (label, x) in (("PotentialReducedTensor", prt), ("VelocityReducedMatrix", vrm),
                        ("FullyReducedTensor", frt))
                    out = Array{DT}(undef, size(x))
                    n = allocations((r, o) -> copyto!(o, r), x, out)
                    # Coverage instrumentation allocates on every counted line, so under it the
                    # assertion does not run rather than being weakened.
                    if Base.JLOptions().code_coverage == 0
                        n == 0 || error("\$label materialised \$n bytes at \$g, m = \$m")
                    end
                    # the barrier sees an allocation when there is one, so the zeros are not
                    # vacuous
                    allocations((r, o) -> collect(r), x, out) > 0 ||
                        error("the barrier sees no allocation")
                end

                for (x, inds) in ((prt, (2, 1, 1)), (vrm, (2, 1)), (frt, (2, 1, 1)))
                    n = entry((r, i) -> r[i...], x, inds)
                    if Base.JLOptions().code_coverage == 0
                        n == 0 || error("one entry of \$x allocated \$n bytes at \$g, m = \$m")
                    end
                end
            end
            """
        local cmd = `$(Base.julia_cmd()) --startup-file=no --check-bounds=auto --project=$(Base.active_project()) -e $prog`
        @test success(cmd)
    end
end
