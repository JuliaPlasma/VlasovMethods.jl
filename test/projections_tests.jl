using AdaptiveRejectionSampling
using PoissonSolvers
using SimpleSplines
using Test
using VlasovMethods
using VlasovMethods: projection!

@testset "Projections" begin
    npart = 1000000
    ncells = 32
    order = 5

    domain = (0.0, 1.0)
    pbasis = PeriodicBasisSpline(domain, order, ncells)
    potential = Potential(pbasis)

    μ = 0.0
    σ = 2.0
    f = x -> exp(-0.5 * (4π * x - μ - 2π)^2 / σ^2) * sqrt(π * σ^2) / sqrt(2)

    dist = ParticleDistribution(1, 1, npart)
    dist.particles.x .= run_sampler!(RejectionSampler(f, domain, max_segments = 5), npart)'
    dist.particles.w .= (ones(npart) ./ npart)'

    projection!(potential, dist)

    # Deposition yields the load vector ∫ρφᵢ. The density's spline coefficients are what the mass
    # matrix maps that to, which is the same step `l2_projection` takes after its own contraction.
    quadrature = SplineQuadrature(pbasis)
    ρ = Spline(pbasis, mass_factorization(quadrature) \ PoissonSolvers.rhs(potential))

    x = domain[begin]:0.1:domain[end]

    cutoff = 2
    @test f.(x)[(begin + cutoff):(end - cutoff)]≈ρ.(x)[(begin + cutoff):(end - cutoff)] atol=5e-2
end
