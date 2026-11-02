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
    @test isapprox(v̇[1, :], (cache.K1 \ cache.LJ) .* landau.ν; rtol = 8 * eps(T))
    @test isapprox(v̇[2, :], (cache.K2 \ cache.LJ) .* landau.ν; rtol = 8 * eps(T))

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

@testset "compute_L! rescales the values of X in place and keeps its sparsity pattern" begin
    landau, pdist = landau_setup(Float64, 16)
    cache = landau.cache[Float64]
    sdist = cache.sdist
    projection(pdist.particles.v, pdist, sdist)
    D = map(copy, cache.D)

    VM.compute_L!(zeros(121, 121), sdist, landau)
    for c in 1:2
        # the same structure as `D`, and the values of `D` scaled by the weighted density
        @test cache.X[c].colptr == D[c].colptr
        @test cache.X[c].rowval == D[c].rowval
        @test cache.D[c] == D[c]
        @test cache.X[c] == D[c] * Diagonal(cache.s)
    end

    # a second call with a different density rescales the values again, on the same pattern
    colptr = map(X -> copy(X.colptr), cache.X)
    rowval = map(X -> copy(X.rowval), cache.X)
    sdist.coefficients .*= 2
    VM.compute_L!(zeros(121, 121), sdist, landau)
    for c in 1:2
        @test cache.X[c].colptr == colptr[c]
        @test cache.X[c].rowval == rowval[c]
        @test cache.X[c] == D[c] * Diagonal(cache.s)
    end
end

@testset "a Landau on a one-dimensional velocity space does not build" begin
    # `LandauCache` builds the buffers of a two-dimensional velocity space, and no method of
    # `compute_L!` takes another, so the constructor raises a `MethodError`
    pdist1 = ParticleDistribution(1, 1, 10)
    sdist1 = SplineDistribution(1, 1, 21, 4, (-6.0, 6.0), 0, Free())
    @test_throws MethodError Landau(pdist1, CollisionEntropy(sdist1))
end

@testset "compute_L! refuses a spline distribution on another quadrature" begin
    pd = ParticleDistribution(1, 2, 100)
    g = v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π) + 0.05
    sd1 = SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free())
    sd2 = SplineDistribution(1, 2, 4, 2, (-2.5, 2.5), 0, Free())
    project_function(g, sd1)
    project_function(g, sd2)
    landau1 = Landau(pd, CollisionEntropy(sd1))
    M = length(sd1)
    @test length(sd2) == M

    # the same sizes on another domain: the stored tabulations are not those of `sd2`'s grid
    err = try
        VM.compute_L!(zeros(M, M), sd2, landau1)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test startswith(err.msg, "compute_L!: `sdist` does not share the quadrature")

    # `similar` keeps the quadrature, so a copy of the cache's distribution is accepted
    sd1b = similar(Float64, sd1)
    sd1b.coefficients .= sd1.coefficients
    @test VM.compute_L!(zeros(M, M), sd1b, landau1) ==
          VM.compute_L!(zeros(M, M), sd1, landau1)

    # the same basis with another quadrature: sharing the basis is not enough
    sd3 = SplineDistribution(1, 2, sd1.basis, TensorProductQuadrature(sd1.basis; nq = 4),
        copy(sd1.coefficients))
    @test_throws ArgumentError VM.compute_L!(zeros(M, M), sd3, landau1)
end

@testset "compute_J! refuses, and compute_K! takes, a spline distribution on another quadrature" begin
    pd = ParticleDistribution(1, 2, 100)
    pd.particles.v .= [2 * sin(i * α) for i in 1:2, α in 1:100]
    pd.particles.w .= 1 / 100
    g = v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π) + 0.05
    sd1 = SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free())
    project_function(g, sd1)
    landau1 = Landau(pd, CollisionEntropy(sd1))
    M1 = length(sd1)
    N = size(pd.particles.v, 2)
    refuses_J(sd) =
        try
            VM.compute_J!(zeros(length(sd)), sd, landau1)
            nothing
        catch e
            e isa ArgumentError &&
                startswith(e.msg, "compute_J!: `sdist` does not share the quadrature")
        end

    # the same sizes on another domain, and another number of cells
    foreign = (SplineDistribution(1, 2, 4, 2, (-2.5, 2.5), 0, Free()),
        SplineDistribution(1, 2, 6, 2, (-3.0, 3.0), 0, Free()))
    for sd in foreign
        project_function(g, sd)
        @test refuses_J(sd) === true
    end

    # `similar` keeps the quadrature, so a copy of the cache's distribution is accepted
    sd1b = similar(Float64, sd1)
    sd1b.coefficients .= sd1.coefficients
    @test VM.compute_J!(zeros(M1), sd1b, landau1) == VM.compute_J!(zeros(M1), sd1, landau1)

    # the same basis with another quadrature: sharing the basis is not enough
    sd3 = SplineDistribution(1, 2, sd1.basis, TensorProductQuadrature(sd1.basis; nq = 4),
        copy(sd1.coefficients))
    @test refuses_J(sd3) === true

    # `compute_K!` reads no quadrature, so it takes any `sdist` of the cache's basis order and
    # gives the `K` of a `Landau` built on that distribution
    for sd in (foreign..., sd3)
        M = length(sd)
        K = VM.compute_K!(zeros(M, N), zeros(M, N), pd.particles.v, sd, landau1)
        landau = Landau(pd, CollisionEntropy(sd))
        @test K == VM.compute_K!(zeros(M, N), zeros(M, N), pd.particles.v, sd, landau)
        @test any(!iszero, K[1]) && any(!iszero, K[2])
    end
end
