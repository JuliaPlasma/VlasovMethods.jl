using ParticleMethods: ParticleList
using Random
using Statistics
using Test
using VlasovMethods
using VlasovMethods: ImportanceSampling

# The mean-of-a-sample assertions below use `atol = 3.5/sqrt(np)`. For a uniform distribution of
# width w the standard error of the mean is w/sqrt(12 np), so on the widest domain here that
# tolerance is about 3.0 standard errors — each such assertion therefore fails about 0.27 % of
# the time, and there are 120 of them across the uniform and shifted-uniform testsets. That is
# roughly a one-in-four chance that some assertion fails on any given run.
#
# Seeding fixes the stream, so the outcome is deterministic and a failure is reproducible rather
# than intermittent. The seed is set here rather than left to the caller because this file runs
# as part of the suite, after other testsets have already consumed from the global RNG.
#
# The 3σ tolerance is deliberate. Widening it to 5/sqrt(np) would make the tests robust rather
# than merely reproducible, but that is a change to what they assert.
Random.seed!(0x5c1f2a08)

npart_list = [10, 100, 1000, 10000, 100000]
domainx = (-5.0, 5.0)

@testset "1D Particle Distribution Tests - Normal" begin
    for np in npart_list
        dist = initialize!(ParticleDistribution(1, 1, np), NormalDistribution())

        @test eltype(dist) == Float64
        @test size(dist) == (3, np)
        @test length(dist.particles.v) == np

        @test mean(dist.particles.v) ≈ 0.0 atol = 3.5/sqrt(np)
        @test sum(dist.particles.w) ≈ 1.0 atol = 1e-8
    end
end

@testset "1D Particle Distribution Tests - Shifted Normal" begin
    for np in npart_list
        for shift in -2.0:1.0:+2.0
            dist = initialize!(ParticleDistribution(1, 1, np), ShiftedNormalV(domainx, shift))

            @test eltype(dist) == Float64
            @test size(dist) == (3, np)
            @test length(dist.particles.v) == np

            @test mean(dist.particles.v) ≈ shift atol = 3.5/sqrt(np)
            @test sum(dist.particles.w) ≈ 1.0 atol = 1e-8
        end
    end
end

@testset "1D Particle Distribution Tests - Double Maxwellian" begin
    for np in npart_list
        for shift in -2.0:1.0:+2.0
            dist = initialize!(ParticleDistribution(1, 1, np), DoubleMaxwellian(shift = shift))

            @test eltype(dist) == Float64
            @test size(dist) == (3, np)
            @test length(dist.particles.v) == np

            @test mean(dist.particles.v) ≈ 0.0 atol = 3.5/sqrt(np)
            @test sum(dist.particles.w) ≈ 1.0 atol = 1e-8
        end
    end
end

@testset "1D Particle Distribution Tests - Uniform" begin
    domain_list = [(-2.0, 2.0), (-1.0, 1.0), (-2.0, 0.0), (0.0, +2.0)]

    for np in npart_list
        for domainv in domain_list
            dist = initialize!(ParticleDistribution(1, 1, np), UniformDistribution(domainx, domainv))

            @test eltype(dist) == Float64
            @test size(dist) == (3, np)
            @test length(dist.particles.v) == np

            @test mean(dist.particles.v) ≈ (domainv[1] + domainv[2])/2 atol = 3.5/sqrt(np)
            @test sum(dist.particles.w) ≈ 1.0 atol = 1e-8
        end
    end
end

@testset "1D Particle Distribution Tests - Shifted Uniform" begin
    domain_list = [(-2.0, 2.0), (-1.0, 1.0), (-2.0, 0.0), (0.0, +2.0)]

    for np in npart_list
        for shift in -2.0:1.0:+2.0
            for domainv in domain_list
                dist = initialize!(ParticleDistribution(1, 1, np),
                    ShiftedUniformDistribution(domainx, domainv, shift))

                @test eltype(dist) == Float64
                @test size(dist) == (3, np)
                @test length(dist.particles.v) == np

                @test mean(dist.particles.v) ≈ (domainv[1] + domainv[2])/2 + shift atol = 3.5/sqrt(np)
                @test sum(dist.particles.w) ≈ 1.0 atol = 1e-8
            end
        end
    end
end

# The 1D1V distribution of `np` particles of the element type `T`, in the layout of
# `ParticleDistribution(1, 1, np)`, which is always `Float64`.
function particle_distribution(::Type{T}, np) where {T}
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    return ParticleDistribution(1, 1, ParticleList(zeros(T, 3, np); variables = vars))
end

@testset "1D Particle Distribution Tests - Bump on Tail" begin
    params = BumpOnTail()
    L = 2π / params.κ

    # The accept–reject draw, on a `Float64` and on a `Float32` distribution: every particle is
    # drawn inside the periodic domain with the equal weight `L / np`, and about the share `α` of
    # them is in the tail at `v₀ = 4.5`. In the bulk, `v > 3` is a share of about `1e-3`.
    for dist in (ParticleDistribution(1, 1, 1000), particle_distribution(Float32, 1000))
        Random.seed!(1234)
        @test initialize!(dist, params) === dist
        x, v = vec(dist.particles.x), vec(dist.particles.v)
        @test all(!iszero, x)
        @test all(0 .≤ x .≤ L)
        @test sum(dist.particles.w) ≈ L
        @test count(>(3), v) / length(v) ≈ params.α atol = 0.03
    end

    # The importance-sampling draw takes every uniform proposal and weights it by the `x`-marginal
    # `1 - ε cos(κx)`, whose integral over the domain is `L`.
    Random.seed!(1234)
    dist = ParticleDistribution(1, 1, 1000)
    @test initialize!(dist, params, ImportanceSampling()) === dist
    x, w = vec(dist.particles.x), vec(dist.particles.w)
    @test all(0 .≤ x .≤ L)
    @test w ≈ (1 .- params.ε .* cos.(params.κ .* x)) .* L ./ 1000
    @test sum(w) ≈ L rtol = 1e-2
    @test count(>(3), vec(dist.particles.v)) / 1000 ≈ params.α atol = 0.03
end
