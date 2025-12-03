using Test
using VlasovMethods
using Statistics

npart_list = [10, 100, 1000, 10000, 100000]
domainx = (-5., 5.)


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
        for shift in -2.:1.0:+2.
        
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
        for shift in -2.:1.0:+2.
            dist = initialize!(ParticleDistribution(1, 1, np), DoubleMaxwellian(shift=shift))
            
            @test eltype(dist) == Float64
            @test size(dist) == (3, np)
            @test length(dist.particles.v) == np

            @test mean(dist.particles.v) ≈ 0.0 atol = 3.5/sqrt(np)
            @test sum(dist.particles.w) ≈ 1.0 atol = 1e-8
        end
    end

end

@testset "1D Particle Distribution Tests - Uniform" begin

    domain_list = [(-2., 2.), (-1., 1.), (-2., 0.), (0., +2.)]

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

    domain_list = [(-2., 2.), (-1., 1.), (-2., 0.), (0., +2.)]

    for np in npart_list
        for shift in -2.:1.0:+2.
            for domainv in domain_list

                dist = initialize!(ParticleDistribution(1, 1, np), ShiftedUniformDistribution(domainx, domainv, shift))
                
                @test eltype(dist) == Float64
                @test size(dist) == (3, np)
                @test length(dist.particles.v) == np

                @test mean(dist.particles.v) ≈ (domainv[1] + domainv[2])/2 + shift atol = 3.5/sqrt(np)
                @test sum(dist.particles.w) ≈ 1.0 atol = 1e-8
                
            end

        end
    end

end