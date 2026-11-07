# The 2D reference cases of the Landau collision methods: the seeded bi-Maxwellian, double
# Maxwellian and BKW profile of `src/examples/`, the two closed forms that
# `scripts/reference_forms.jl` derives from the moment equations, and the isotropisation rate that
# the particle-pair scheme of `γ = 0` reproduces at the mesh of the plan's measurement.
#
# Every standard error is the convention of the spec — `σ/√N` for a mean, `σ²√(2/N)` for a variance
# and `σ₁σ₂/√N` for a covariance — and every bound on a sampled deviation is four of them; the
# fourth moment's bound is four of the standard error that the sample itself gives.
#
# The rate test draws its sample with `StableRNG`, as the probe of the plan's measurement does, and
# not through `AnisotropicBiMaxwellian`, whose `seed` streams `MersenneTwister`: the rate it compares
# with was measured on that sample, and `StableRNGs` is a test-only dependency, so the sample here is
# the one `bi_maxwellian(1600)` of the plan's checks draws. The type's own sample is checked by
# "the sample's moments" below.

using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, TensorProductBasis,
                     TensorProductQuadrature
using StableRNGs: StableRNG
using Test
using VlasovMethods

# The functions of the two derivations: `moment_rate` and the BKW profile, whose closed forms this
# file compares against its own.
include(joinpath(@__DIR__, "..", "..", "scripts", "reference_forms.jl"))

"The number of particles of a sampling test, at which a mean's standard error is `σ/316`."
const NPARTICLES = 100_000

"The seed of every sample this file draws."
const SEED = 20261107

# The first and second moments of the particles `v` (2 × N) with the weights `w`:
# `T_ij = Σ_α w_α (v_{i,α} - u_i)(v_{j,α} - u_j) / Σ_α w_α` and `u = Σ_α w_α v_α / Σ_α w_α`.
function moments(v::AbstractMatrix, w::AbstractVector)
    total = sum(w)
    mean = v * w ./ total
    centered = v .- mean
    return (mean = mean,
        T11 = sum(w .* centered[1, :] .^ 2) / total,
        T22 = sum(w .* centered[2, :] .^ 2) / total,
        T12 = sum(w .* centered[1, :] .* centered[2, :]) / total)
end

"The particles of a `ParticleDistribution` of `npart` particles and two velocity components."
function particles(npart)
    return ParticleDistribution(1, 2, npart)
end

# The particle-pair model at the Maxwell molecules `γ = 0` of the particles `v` (2 × N) with the
# weights `w`, on quadratic splines of `ncell × ncell` cells of `[-L, L]²`.
function reference_model(v::AbstractMatrix, w::AbstractVector; ncell, L)
    npart = size(v, 2)
    variables = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    list = ParticleList(zeros(eltype(v), 4, npart); variables)
    dist = ParticleDistribution(1, 2, list)
    dist.particles.v .= v
    dist.particles.w[1, :] .= w
    axis = BSplineBasis(UniformMesh(ncell, -L .. L), 2, Free())
    basis = TensorProductBasis((axis, axis))
    splines = SplineDistribution(1, 2, basis, TensorProductQuadrature(basis),
        zeros(eltype(v), size(basis)...))
    return PairwiseLandau(dist, CollisionEntropy(splines); γ = 0, ν = 1.0)
end

# The rate `-d ln|T₁₁ - T₂₂|/dt` of the anisotropy gap of `model` over `nsteps` steps of `Δt`, as
# the slope of `ln|T₁₁ - T₂₂|` against `t`, with the drift of the trace over the same run.
function fitted_rate(model, Δt, nsteps; abstol = 1e-7)
    w = collect(vec(model.landau.dist.particles.w))
    first_moment = moments(model.landau.dist.particles.v, w)
    trace = first_moment.T11 + first_moment.T22
    ts = Float64[0.0]
    gaps = Float64[first_moment.T11 - first_moment.T22]
    for k in 1:nsteps
        discrete_gradient_step!(model, Δt; abstol)
        m = moments(model.landau.dist.particles.v, w)
        push!(ts, k * Δt)
        push!(gaps, m.T11 - m.T22)
    end
    last_moment = moments(model.landau.dist.particles.v, w)
    ordinates = log.(abs.(gaps))
    mean_t = sum(ts) / length(ts)
    mean_ln = sum(ordinates) / length(ordinates)
    slope = sum((ts .- mean_t) .* (ordinates .- mean_ln)) / sum((ts .- mean_t) .^ 2)
    return (rate = -slope, trace_drift = abs(last_moment.T11 + last_moment.T22 - trace))
end

# The sample of the rate test: `bi_maxwellian(1600)` of the plan's checks, a centred sample of the
# two variances with equal weights, its velocities 2 × npart.
function rate_sample(npart; variance_one = 0.5, variance_two = 0.3, seed = SEED)
    rng = StableRNG(seed)
    v = zeros(2, npart)
    v[1, :] .= sqrt(variance_one) .* randn(rng, npart)
    v[2, :] .= sqrt(variance_two) .* randn(rng, npart)
    return v, fill(1 / npart, npart)
end

@testset "the sample's moments" begin
    variance_one, variance_two = 0.5, 0.3
    npart = NPARTICLES
    dist = particles(npart)
    # the weights the distribution holds are the ones the sample is measured with
    dist.particles.w[1, :] .= 1 / npart
    initialize!(dist, AnisotropicBiMaxwellian(variance_one, variance_two; seed = SEED))
    w = collect(vec(dist.particles.w))
    m = moments(dist.particles.v, w)

    @test sum(w) ≈ 1
    @test abs(m.mean[1]) ≤ 4 * sqrt(variance_one) / √npart
    @test abs(m.mean[2]) ≤ 4 * sqrt(variance_two) / √npart
    @test abs(m.T11 - variance_one) ≤ 4 * variance_one * √(2 / npart)
    @test abs(m.T22 - variance_two) ≤ 4 * variance_two * √(2 / npart)
    @test abs(m.T12) ≤ 4 * sqrt(variance_one * variance_two) / √npart

    # the four seeds of the measurement: no deviation of any exceeds four standard errors
    for seed in (SEED, 1, 2, 3)
        s = initialize!(particles(npart), AnisotropicBiMaxwellian(variance_one, variance_two; seed))
        ms = moments(s.particles.v, w)
        @test maximum(abs,
            (ms.mean[1] / (sqrt(variance_one) / √npart),
                ms.mean[2] / (sqrt(variance_two) / √npart),
                (ms.T11 - variance_one) / (variance_one * √(2 / npart)),
                (ms.T22 - variance_two) / (variance_two * √(2 / npart)),
                ms.T12 / (sqrt(variance_one * variance_two) / √npart))) ≤ 4
    end

    # the seed holds the sample, and another seed gives another one
    @test initialize!(particles(npart), AnisotropicBiMaxwellian(variance_one, variance_two;
        seed = SEED)).particles.v == dist.particles.v
    @test initialize!(particles(npart), AnisotropicBiMaxwellian(variance_one, variance_two;
        seed = SEED + 1)).particles.v != dist.particles.v

    # the weights the distribution holds are the ones the sample keeps
    held = particles(10)
    held.particles.w[1, :] .= 0.1 .* (1:10)
    weights = collect(vec(held.particles.w))
    initialize!(held, AnisotropicBiMaxwellian(variance_one, variance_two; seed = SEED))
    @test collect(vec(held.particles.w)) == weights

    # the control: a sample of the variance 0.35 in place of 0.3 breaks the bound of the second
    # variance, which a sampler of the wrong variance cannot pass
    wrong = initialize!(particles(npart),
        AnisotropicBiMaxwellian(variance_one, 0.35; seed = SEED))
    @test abs(moments(wrong.particles.v, w).T22 - variance_two) >
          4 * variance_two * √(2 / npart)

    # the reference case is two-dimensional
    @test_throws DimensionMismatch initialize!(ParticleDistribution(1, 1, 10),
        AnisotropicBiMaxwellian(variance_one, variance_two))
    @test_throws DimensionMismatch initialize!(ParticleDistribution(1, 3, 10),
        AnisotropicBiMaxwellian(variance_one, variance_two))
end

@testset "the double Maxwellian's moments" begin
    npart = NPARTICLES
    weights = (0.6, 0.4)
    variances = ((0.5, 0.3), (0.3, 0.5))
    means = ((0.4, 0.0), (-0.4, 0.0))
    # the mixture's exact mean and variances, about that mean
    exact_mean = [sum(weights[c] * means[c][d] for c in 1:2) for d in 1:2]
    exact_variance = [sum(weights[c] * (variances[c][d] + (means[c][d] - exact_mean[d])^2)
                      for c in 1:2) for d in 1:2]

    @test exact_mean[1] ≈ 0.08
    @test exact_mean[2] == 0
    @test exact_variance[1] ≈ 0.5736
    @test exact_variance[2] ≈ 0.38

    dist = initialize!(particles(npart),
        DoubleMaxwellian2D(; weights, variances, means, seed = SEED))
    w = collect(vec(dist.particles.w))
    m = moments(dist.particles.v, w)

    @test sum(w) ≈ 1
    @test abs(m.mean[1] - exact_mean[1]) ≤ 4 * sqrt(exact_variance[1]) / √npart
    @test abs(m.mean[2] - exact_mean[2]) ≤ 4 * sqrt(exact_variance[2]) / √npart
    @test abs(m.T11 - exact_variance[1]) ≤ 4 * exact_variance[1] * √(2 / npart)
    @test abs(m.T22 - exact_variance[2]) ≤ 4 * exact_variance[2] * √(2 / npart)

    # the four seeds of the measurement
    for seed in (SEED, 1, 2, 3)
        s = initialize!(particles(npart), DoubleMaxwellian2D(; weights, variances, means, seed))
        ms = moments(s.particles.v, w)
        @test maximum(abs,
            ((ms.mean[1] - exact_mean[1]) / (sqrt(exact_variance[1]) / √npart),
                (ms.mean[2] - exact_mean[2]) / (sqrt(exact_variance[2]) / √npart),
                (ms.T11 - exact_variance[1]) / (exact_variance[1] * √(2 / npart)),
                (ms.T22 - exact_variance[2]) / (exact_variance[2] * √(2 / npart)))) ≤ 4
    end

    # the seed holds the sample, and another seed gives another one
    @test initialize!(particles(npart),
        DoubleMaxwellian2D(; weights, variances, means,
            seed = SEED)).particles.v == dist.particles.v
    @test initialize!(particles(npart),
        DoubleMaxwellian2D(; weights, variances, means,
            seed = SEED + 1)).particles.v != dist.particles.v

    # the control: the two weights swapped, which moves the mean to -0.08 and breaks the bound of
    # the first component, which a sampler of the wrong weights cannot pass
    swapped = initialize!(particles(npart),
        DoubleMaxwellian2D(; weights = (0.4, 0.6), variances, means, seed = SEED))
    @test abs(moments(swapped.particles.v, w).mean[1] - exact_mean[1]) >
          4 * sqrt(exact_variance[1]) / √npart
end

@testset "the derivation and the moment equations" begin
    npart = NPARTICLES
    # the parameter of the profile, at the kernel U = 16 A of this package
    for t₀ in (0.1, 0.3, 0.6, 0.9)
        @test isapprox(bkw_parameter(t₀), 1 - exp(-2 * t₀) / 2; rtol = 1e-14)
    end

    # the profile's mass, energy, temperature and fourth moment, and the equation of its fourth
    # moment, against the moment equations at the four times
    for t₀ in (0.1, 0.3, 0.6, 0.9)
        K = bkw_parameter(t₀)
        squared = v -> v[1]^2 + v[2]^2
        @test isapprox(bkw_moment(K, v -> 1.0), 1.0; rtol = 1e-10)
        @test isapprox(bkw_moment(K, squared), 2.0; rtol = 1e-10)
        @test isapprox(bkw_moment(K, squared) / 2, 1.0; rtol = 1e-10)
        @test isapprox(bkw_moment(K, v -> squared(v)^2), 8 * K * (2 - K); rtol = 1e-10)
        @test isapprox(
            moment_rate(v -> (4 * squared(v) * v[1], 4 * squared(v) * v[2]),
                bkw_log_gradient(K), bkw_weight(K), K, K),
            32 * (1 - K)^2;
            rtol = 1e-10)
    end

    # the control: the Boltzmann rate of the paper's kernel in place of the sixteenfold one, which
    # the fourth moment's equation refuses
    t₀, K = 0.3, bkw_parameter(0.3)
    squared = v -> v[1]^2 + v[2]^2
    fourth_rate = moment_rate(v -> (4 * squared(v) * v[1], 4 * squared(v) * v[2]),
        bkw_log_gradient(K), bkw_weight(K), K, K)
    @test abs(16 * (1 - K) * exp(-t₀ / 8) - fourth_rate) > 1e-3 * abs(fourth_rate)

    # the type's sample of the profile at two times: its two means and variances, and its fourth
    # moment, each within four standard errors
    for t₀ in (0.2, 0.6)
        K = bkw_parameter(t₀)
        dist = initialize!(particles(npart), BKWDistribution(t₀; seed = SEED))
        w = collect(vec(dist.particles.w))
        m = moments(dist.particles.v, w)
        @test sum(w) ≈ 1
        @test abs(m.mean[1]) ≤ 4 / √npart
        @test abs(m.mean[2]) ≤ 4 / √npart
        @test abs(m.T11 - 1) ≤ 4 * √(2 / npart)
        @test abs(m.T22 - 1) ≤ 4 * √(2 / npart)

        # the fourth moment is the profile's, which the first two moments cannot tell from a
        # Gaussian's
        squared_radius = dist.particles.v[1, :] .^ 2 .+ dist.particles.v[2, :] .^ 2
        fourth = sum(w .* squared_radius .^ 2)
        fourth_error = sqrt(sum((w .* (squared_radius .^ 2 .- fourth)) .^ 2))
        @test abs(fourth - 8 * K * (2 - K)) ≤ 4 * fourth_error
    end

    # the control at t₀ = 0.2: a Gaussian of temperature one, whose first two moments are the
    # profile's and whose fourth moment is 8, breaks the bound of the fourth moment
    t₀, K = 0.2, bkw_parameter(0.2)
    dist = initialize!(particles(npart), BKWDistribution(t₀; seed = SEED))
    w = collect(vec(dist.particles.w))
    squared_radius = dist.particles.v[1, :] .^ 2 .+ dist.particles.v[2, :] .^ 2
    fourth = sum(w .* squared_radius .^ 2)
    fourth_error = sqrt(sum((w .* (squared_radius .^ 2 .- fourth)) .^ 2))
    gaussian = randn(StableRNG(SEED), 2, npart)
    gaussian_radius = gaussian[1, :] .^ 2 .+ gaussian[2, :] .^ 2
    @test abs(sum(w .* gaussian_radius .^ 2) - 8 * K * (2 - K)) > 4 * fourth_error
    @test abs(8 - 8 * K * (2 - K)) > 4 * fourth_error

    # the profile is no distribution at a negative time, where its quadratic factor changes sign
    @test_throws ArgumentError BKWDistribution(-0.1)

    # the script's own sample of the profile at t₀ = 0.2: each deviation of its first two moments,
    # and of its fourth, within four standard errors
    @test all(r -> abs(r.value) ≤ 4, sample_rows(t₀))
end

@testset "the rate's derivation" begin
    for (T₁, T₂) in ((1.2, 1.0), (1.5, 1.0), (2.0, 1.0), (3.0, 1.0))
        r = rate_checks(T₁, T₂)
        @test isapprox(r.rate₁₁, 4 * (T₂ - T₁); rtol = 1e-10)
        @test isapprox(r.rate₂₂, 4 * (T₁ - T₂); rtol = 1e-10)
        @test isapprox(r.gap, -8 * (T₁ - T₂); rtol = 1e-10)

        # the control: the constant 4 in place of 8 misses by more than 1e-3 relative
        @test abs(-4 * (T₁ - T₂) - r.gap) > 1e-3 * abs(r.gap)
    end
end

@testset "the rate at 16 × 16 cells, and the 8 × 8 control" begin
    v, w = rate_sample(1600)

    # the mesh of the plan's measurement, at which the scheme reproduces the closed rate 8
    fine = fitted_rate(reference_model(v, w; ncell = 16, L = 3.0), 1e-5, 30)
    @test abs(fine.rate - 8) ≤ 0.05 * 8
    @test fine.trace_drift ≤ 1e-12

    # the control: 8 × 8 cells of [-4, 4]², whose width is a Gaussian's, lose the rate, outside the
    # clamp of the fine mesh
    coarse = fitted_rate(reference_model(v, w; ncell = 8, L = 4.0), 1e-4, 30)
    @test coarse.rate < 5
    @test abs(coarse.rate - 8) > 0.05 * 8
end
