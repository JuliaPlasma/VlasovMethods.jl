using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, SplineQuadrature, nbasis
using Test
using VlasovMethods

const VM = VlasovMethods

function setup(::Type{T}) where {T}
    N = 64
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, N); variables = vars))
    pdist.particles.v[1, :] .= range(T(-1.5), T(1.5); length = N)
    pdist.particles.w[1, :] .= one(T) / N
    axis = BSplineBasis(UniformMesh(17, T(-2) .. T(2)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))
    return pdist, sdist
end

@testset "CacheDict" begin
    # `CacheDict(p)` stores its parent under the element type of `p`, so the parent is found under
    # that type, and a lookup of it builds no second cache
    @testset "$(nameof(M)) keeps its parent under its element type, $T" for M in (
            LenardBernstein, ConservativeLenardBernstein,
            RescaledConservativeLenardBernstein, MetriplecticLenardBernstein),
        T in (Float32, Float64)

        pdist, sdist = setup(T)
        model = M(pdist, CollisionEntropy(sdist))
        p = parent(model.cache)
        @test eltype(p) === T
        @test eltype(typeof(p)) === T
        @test model.cache[T] === p
        @test length(model.cache.caches) == 1

        # another element type builds one cache of its own, once
        other = T === Float64 ? Float32 : Float64
        c = model.cache[other]
        @test eltype(c) === other
        @test model.cache[other] === c
        @test length(model.cache.caches) == 2
    end
end
