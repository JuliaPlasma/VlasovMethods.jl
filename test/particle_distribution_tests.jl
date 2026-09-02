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
# than intermittent. It matters here specifically because this file did not used to run at all:
# `runtests.jl` included only the two spline test files, so nothing depended on the RNG state at
# entry, and adding it to the suite made it sensitive to whatever the preceding testset had
# consumed.
#
# The 3σ tolerance is left as it was. Widening it to 5/sqrt(np) would make the tests robust
# rather than merely reproducible, but that is a change to what they assert.
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
