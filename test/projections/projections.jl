using PoissonSolvers
using Random
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

    # The sample is drawn from the global RNG, so without this the tolerance below is asserted
    # against a different sample on every run.
    Random.seed!(1234)

    μ = 0.0
    σ = 2.0
    f = x -> exp(-0.5 * (4π * x - μ - 2π)^2 / σ^2) * sqrt(π * σ^2) / sqrt(2)

    # Substituting 4πx - 2π = 2π(2x - 1) turns `f` into √(2π)·exp(-2π²(x - ½)²), which is the
    # normal density with mean 1/2 and standard deviation 1/(2π) — both the exponent and the
    # normalising constant match. So the sample is drawn from that normal directly. A draw outside
    # the domain is replaced; the edges sit at ±π standard deviations, which rejects about 0.17%,
    # and the sampler this replaces was truncated to the same domain.
    xmean = 0.5
    xstd = 1 / 2π
    samples = Vector{Float64}(undef, npart)
    for i in eachindex(samples)
        x = xmean + xstd * randn()
        while !(domain[begin] < x < domain[end])
            x = xmean + xstd * randn()
        end
        samples[i] = x
    end

    dist = ParticleDistribution(1, 1, npart)
    dist.particles.x .= samples'
    dist.particles.w .= (ones(npart) ./ npart)'

    projection!(potential, dist)

    # Deposition yields the load vector ∫ρφᵢ. The density's spline coefficients are what the mass
    # matrix maps that to, which is the same step `l2_projection` takes after its own contraction.
    quadrature = SplineQuadrature(pbasis)
    ρ = Spline(pbasis, mass_factorization(quadrature) \ PoissonSolvers.rhs(potential))

    x = domain[begin]:0.1:domain[end]

    # The error here is the sampling error of a million draws, measured at ~1.3e-2 across seeds
    # and Julia versions. The tolerance keeps a factor of two over that and no more.
    @test f.(x)≈ρ.(x) atol=2.5e-2

    @testset "deposition is exact" begin
        # The reconstruction above cannot see the wrapping: a deposition that drops the wrapped
        # contributions instead of folding them still meets that tolerance. A direct sum over
        # every basis function is the definition of the load vector, so it has no cancellation
        # to hide behind, and it covers the recombined basis, whose block is padded with zeros
        # at indices past `nbasis`.
        # Every second point of this range is a breakpoint, and the first and last are the
        # domain's own ends.
        nsmall = 2 * ncells + 1

        for b in (PeriodicBasisSpline(domain, order, ncells),
            DirichletBasisSpline(domain, order, ncells))
            small = ParticleDistribution(1, 1, nsmall)
            # The ends and the breakpoints are where the wrap and the padding decide the
            # answer, so they are sampled rather than avoided.
            small.particles.x .= collect(range(domain[begin], domain[end], length = nsmall))'
            small.particles.w .= (collect(1:nsmall) ./ nsmall)'

            p = Potential(b)
            projection!(p, small)

            reference = [sum(w * evaluate(b, i, x)
                         for (x, w) in zip(small.particles.x, small.particles.w))
                         for i in 1:nbasis(b)]

            @test PoissonSolvers.rhs(p) ≈ reference
        end
    end
end
