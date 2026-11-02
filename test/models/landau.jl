using LinearAlgebra
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, TensorProductBasis,
                     TensorProductQuadrature
using Test
using VlasovMethods

const VM = VlasovMethods

# `n1²` particles of uniform weight at the cell centres of an `n1 × n1` grid on the velocity
# domain `[-2, 2]²`, so that the L² projection stays positive (KNOWN_ISSUES.md, K8), and a cubic
# clamped basis of 8 × 8 cells: M = 121 basis functions, Q = 1600 quadrature points.
function landau_setup(::Type{T}, n1) where {T}
    npart = n1^2
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pdist = ParticleDistribution(1, 2, ParticleList(zeros(T, 4, npart); variables = vars))
    g = [T(-2) + (T(i) - T(1) / 2) * T(4) / n1 for i in 1:n1]
    pdist.particles.v[1, :] .= repeat(g; outer = n1)
    pdist.particles.v[2, :] .= repeat(g; inner = n1)
    pdist.particles.w[1, :] .= one(T) / npart
    axis = BSplineBasis(UniformMesh(8, T(-2) .. T(2)), 3, Free())
    B = TensorProductBasis((axis, axis))
    sdist = SplineDistribution(1, 2, B, TensorProductQuadrature(B), zeros(T, size(B)...))
    return Landau(pdist, CollisionEntropy(sdist)), pdist
end

# The bytes of one call after a warm call, through barriers of fixed arity whose arguments have
# concrete types.
function bytes3(f::F, a::A, b::B, c::C) where {F, A, B, C}
    f(a, b, c)
    return @allocated f(a, b, c)
end
function bytes4(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    f(a, b, c, d)
    return @allocated f(a, b, c, d)
end
function bytes5(f::F, a::A, b::B, c::C, d::D, e::E) where {F, A, B, C, D, E}
    f(a, b, c, d, e)
    return @allocated f(a, b, c, d, e)
end
function bytes6(f::F, a::A, b::B, c::C, d::D, e::E, g::G) where {F, A, B, C, D, E, G}
    f(a, b, c, d, e, g)
    return @allocated f(a, b, c, d, e, g)
end

# The controls: calls of each arity that allocate.
copy_first(a, b, c, d, e) = copy(a)
copy_first(a, b, c, d, e, g) = copy(a)
# A least-squares solve of the vector field, on the four-argument barrier of the vector field.
solve(K, LJ, _, _) = K \ LJ

@testset "Landau buffers, $T, Np = $(n1^2)" for T in (Float32, Float64), n1 in (16, 24)

    landau, pdist = landau_setup(T, n1)
    cache = landau.cache[T]
    sdist = cache.sdist
    @test length(sdist) == 121
    @test size(cache.D[1]) == (121, 1600)
    v = copy(pdist.particles.v)
    v̇ = similar(v)

    vf = bytes4(VM.collisional_vectorfield!, v̇, v, nothing, landau)
    # the two solves on the matrices that the evaluation leaves in the cache
    ceiling = bytes4(solve, cache.K1, cache.LJ, nothing, nothing) +
              bytes4(solve, cache.K2, cache.LJ, nothing, nothing)

    @test @inferred(VM.collisional_vectorfield!(v̇, v, nothing, landau)) === nothing
    @test @inferred(VM.compute_L!(cache.L, sdist, landau)) === cache.L
    @test @inferred(VM.compute_J!(cache.J, sdist, landau)) === cache.J
    @test @inferred(VM.compute_K!(cache.K1, cache.K2, v, sdist, landau)) ==
          (cache.K1, cache.K2)
    @test @inferred(VM._projection!(
        cache.rhs, cache.vals, cache.mwork, v, pdist, sdist)) ===
          sdist.spline

    # the vector field is the minimum-norm solution of the two least-squares problems
    @test v̇[1, :] == (cache.K1 \ cache.LJ) .* landau.ν
    @test v̇[2, :] == (cache.K2 \ cache.LJ) .* landau.ν

    # coverage instrumentation allocates on every counted line, so under it the assertions
    # do not run rather than being weakened
    if Base.JLOptions().code_coverage == 0
        @test bytes3(VM.compute_L!, cache.L, sdist, landau) == 0
        @test bytes3(VM.compute_J!, cache.J, sdist, landau) == 0
        @test bytes5(VM.compute_K!, cache.K1, cache.K2, v, sdist, landau) == 0
        @test bytes6(
            VM._projection!, cache.rhs, cache.vals, cache.mwork, v, pdist, sdist) ==
              0
        @test 0 < vf ≤ ceiling
    end
    # the controls: each barrier sees an allocation
    @test bytes3(projection, v, pdist, sdist) > 0
    @test ceiling > 0
    @test bytes5(copy_first, v, v, v, v, v) > 0
    @test bytes6(copy_first, v, v, v, v, v, v) > 0
end

@testset "the in-place 2D projection is the allocating one" begin
    landau, pdist = landau_setup(Float64, 16)
    cache = landau.cache[Float64]
    v = copy(pdist.particles.v)
    v[1, 1:7] .+= 0.1
    projection(v, pdist, cache.sdist)
    expected = copy(cache.sdist.coefficients)
    fill!(cache.sdist.coefficients, 0)
    VM._projection!(cache.rhs, cache.vals, cache.mwork, v, pdist, cache.sdist)
    @test cache.sdist.coefficients == expected
    # the load array is cleared on each call, so a second call gives the same coefficients
    VM._projection!(cache.rhs, cache.vals, cache.mwork, v, pdist, cache.sdist)
    @test cache.sdist.coefficients == expected
end

@testset "compute_L! reads the gradient tabulations stored in the cache" begin
    landau, pdist = landau_setup(Float64, 16)
    cache = landau.cache[Float64]
    sdist = cache.sdist
    projection(pdist.particles.v, pdist, sdist)

    # `LandauCache` stores the tabulations of its own quadrature grid
    @test cache.D == map(Matrix, VM.gradient_tabulations(sdist))

    L = VM.compute_L!(zeros(121, 121), sdist, landau)
    # `L` is quadratic in the tabulations, and a factor of 2 is exact in floating point, so
    # doubling the stored tables gives exactly 4L when `compute_L!` reads them
    foreach(D -> D .*= 2, cache.D)
    @test VM.compute_L!(zeros(121, 121), sdist, landau) == 4 .* L
    foreach(D -> D ./= 2, cache.D)
    @test VM.compute_L!(zeros(121, 121), sdist, landau) == L
end
