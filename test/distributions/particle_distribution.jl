using Random
using Statistics
using Test
using VlasovMethods

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

@testset "1D Particle Distribution Tests - Bump on Tail" begin
    # The accept–reject draw of `src/examples/bumpontail.jl`: the positions are drawn on the one
    # period `[0, 2π/κ]` and every particle is written, so no particle keeps the zero the
    # constructor gave it, and every weight is the uniform `2π/κ/N`. The seed is set here: the
    # draw consumes the global RNG, and the tests above have already drawn from it.
    κ = 0.3
    Random.seed!(1234)
    dist = initialize!(ParticleDistribution(1, 1, 1000), BumpOnTail(κ = κ))

    x = vec(dist.particles.x)
    @test length(x) == 1000
    @test all(!iszero, x)
    @test all(0 .≤ x .≤ 2π / κ)
    @test sum(dist.particles.w) ≈ 2π / κ atol = 1e-8
    # the draw is the one the seed gives, so a second draw of the same seed is the same cloud
    Random.seed!(1234)
    again = initialize!(ParticleDistribution(1, 1, 1000), BumpOnTail(κ = κ))
    @test vec(again.particles.x) == x

    # the importance-sampling draw of the same file: the second method, which writes every
    # particle at its proposal and weights it by `fₓ` instead of the uniform `2π/κ/N`
    params = BumpOnTail(κ = κ)
    Random.seed!(1234)
    weighted = initialize!(ParticleDistribution(1, 1, 1000), params,
        VlasovMethods.ImportanceSampling())

    xw = vec(weighted.particles.x)
    @test all(!iszero, xw)
    @test all(0 .≤ xw .≤ 2π / κ)
    # the weight is `fₓ(x)·2π/κ/N`, and another grouping of that product differs by a few ulps
    @test vec(weighted.particles.w) ≈
          VlasovMethods.fₓ.(xw, Ref(params)) .* (2π / κ / 1000) rtol = 4eps(Float64)
    # the bump takes `α = 0.1` of the particles to `v₀ = 4.5`, so the mean of `v` is near
    # `α v₀ = 0.45`, and a tail branch that never fires leaves it at zero
    @test mean(vec(weighted.particles.v)) > 0.2
end
