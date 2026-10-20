# The shared part of the drivers of the paper's two example sets, `scripts/bump_on_tail_common.jl`:
# the sweeps, the seeded initial condition, the field of one sample, the accelerations it applies
# and the stride of the analytic set's saves. No test runs a driver; each runs for minutes.

using LinearAlgebra
using ParticleMethods: ParticleList
using PoissonSolvers
using ReducedComplexityModeling: Parameter, ParameterSpace
using Test
using VlasovMethods
using VlasovMethods: coefficients, efield!, energy, update!

include(joinpath(@__DIR__, "..", "..", "scripts", "bump_on_tail_common.jl"))

const SCRIPTS = joinpath(@__DIR__, "..", "..", "scripts")

"The potential of the reference problem of the sweeps, on `16` cubic cells."
reference_potential() = Potential(PeriodicBasisSpline(domain(KAPPA_REF), 3, 16))

"The phrase of the exception a call raises, or a fixed word when it raises nothing."
function refusal(f)
    try
        f()
        return "no exception"
    catch thrown
        return sprint(showerror, thrown)
    end
end

"""
The snapshots of one full-order march of the particles `plt` on `field`, `nt` steps of `dt` at the
scale `χ`, saved at every step.
"""
function sweep_march(plt, field, nt, dt, χ)
    n, nh = size(plt.x, 2), length(coefficients(field))
    SS = Snapshots(Float64, 1, n, nh, nt + 1, ParameterSpace(Parameter(:χ, [χ])))
    integrate_vp(
        plt, (χ = χ,), field, SS, IntegratorParameters(dt, nt, nt + 1, nh, n, 1), 1)
    return SS
end

@testset "the initial condition of the set is seeded" begin
    a = initial_distribution(500, 1234)
    b = initial_distribution(500, 1234)
    c = initial_distribution(500, 4321)
    @test a.particles.x == b.particles.x
    @test a.particles.v == b.particles.v
    @test a.particles.x != c.particles.x

    # the positions lie on one wavelength at `KAPPA_REF`, and the weights add up to its length
    L = 2π / KAPPA_REF
    @test all(x -> 0 ≤ x ≤ L, a.particles.x)
    @test sum(a.particles.w) ≈ L rtol=1e-12

    # the bump: a share near `α = 0.1` of the particles is drawn around `v₀ = 4.5`
    fast = count(>(3), a.particles.v) / 500
    @test 0.05 < fast < 0.15
end

@testset "the sweeps" begin
    # the step count rounds the quotient, where the old script's `Int(div(25, 0.1))` is 249
    @test step_count(25.0, 0.1) == 250
    @test Int(div(25, 0.1)) == 249

    κ = training_kappas(0.1, 0.5, 10)
    @test length(κ) == 10
    @test first(κ) == 0.1
    @test last(κ) == 0.5
    @test allunique(κ)

    # ten seeded test values, sorted, inside the range, and the same for the same seed
    κt = testing_kappas(0.1, 0.5, 10, 1234)
    @test length(κt) == 10
    @test issorted(κt)
    @test all(k -> 0.1 ≤ k < 0.5, κt)
    @test κt == testing_kappas(0.1, 0.5, 10, 1234)
    @test κt != testing_kappas(0.1, 0.5, 10, 4321)

    # the scale of a sample is its wave number over the reference one
    @test march_parameters(KAPPA_REF) == (χ = 1.0,)
    @test march_parameters(0.15).χ == 0.15 / 0.3
end

@testset "the field of a sample is scaled once, at the march's χ" begin
    potential = reference_potential()
    p = (χ = 2.0,)

    # (i) one `ScaledField` over the Poisson field of the caller's potential, at `params.χ`
    f = vp_field(potential, p)
    @test f isa ScaledField
    @test f.field isa PoissonField
    @test f.field.potential === potential
    @test f.χ === p.χ

    # (ii) at `χ = 2` the field and the saved energy are a quarter of the plain field's, exactly,
    # because a division by a power of two is exact
    plt = initial_distribution(400, 1234).particles
    x, w = vec(plt.x), vec(plt.w)
    plain = PoissonField(reference_potential())
    update!(f, x, w, 0.0)
    update!(plain, x, w, 0.0)
    e, e₀ = similar(x), similar(x)
    efield!(f, e, x)
    efield!(plain, e₀, x)
    @test any(!iszero, e₀)
    @test e == e₀ ./ 4

    scaled = sweep_march(plt, vp_field(reference_potential(), p), 4, 0.1, 2.0)
    unscaled = sweep_march(plt, PoissonField(reference_potential()), 4, 0.1, 2.0)
    @test scaled.W[1, 1] > 0
    @test scaled.W[1, 1] == unscaled.W[1, 1] / 4

    # (iii) the field scale is part of the model: the march at `(dt, χ = 2)` with the scaled field
    # is not the plain field's march at `(2dt, χ = 1)`, which has the same step
    twice = sweep_march(plt, PoissonField(reference_potential()), 4, 0.2, 1.0)
    @test scaled.X != twice.X

    # (iv) the drivers build one `ScaledField`, the one of `vp_field`
    hits = [(name, i)
            for name in readdir(SCRIPTS) if endswith(name, ".jl")
            for (i, line) in enumerate(eachline(joinpath(SCRIPTS, name)))
            if occursin("ScaledField(", line)]
    @test length(hits) == 1
    @test first(only(hits)) == "bump_on_tail_common.jl"
end

@testset "the accelerations are the ones the march applied" begin
    plt = initial_distribution(300, 1234).particles
    w = vec(plt.w)
    χ = 2.0
    field = vp_field(reference_potential(), (χ = χ,))
    SS = sweep_march(plt, field, 3, 0.1, χ)
    @test all(iszero, SS.A)
    t = collect(0:3) .* 0.1
    @test accelerations!(SS, field, w, t, 1) === SS

    # each saved acceleration is the plain field at the saved positions, divided by `χ² = 4`
    plain = PoissonField(reference_potential())
    for ts in 1:4
        x = SS.X[1, :, ts, 1]
        update!(plain, x, w, t[ts])
        @test SS.A[1, :, ts, 1] == efield!(plain, similar(x), x) ./ 4
    end
    @test any(!iszero, SS.A[1, :, 4, 1])
end

@testset "the stride of the saved times is a whole number of reduced steps" begin
    @test save_stride(10, 1e-3, 1e-4) == 100
    @test save_stride(100, 1e-3, 1e-4) == 1000
    @test save_stride(3, 0.7, 0.1) == 21

    # the computed ratio is `24.999999999999993`, `1.28 eps` from `25`, relative: inside the bound
    # of `2 eps` and outside a bound half as wide
    @test 3 * 0.7 / 0.084 != 25
    @test abs(3 * 0.7 / 0.084 - 25) > eps(Float64) * 25
    @test save_stride(3, 0.7, 0.084) == 25

    # fractions at the paper's own stride of 100, and half a step and a quarter step at 10⁹
    phrase = "which is not a whole number"
    for args in ((10, 1e-3, 1.5e-4), (10, 1e-3, 3e-4), (1, 100.00000001, 1.0),
        (10, 1e-3, 1.0000000001e-4), (1, 1e9 + 0.25, 1.0), (1, 1e9 + 0.5, 1.0))
        @test_throws ArgumentError save_stride(args...)
        @test occursin(phrase, refusal(() -> save_stride(args...)))
    end

    # a save of no reduced step is refused by its own phrase
    for args in ((0, 1e-3, 1e-4), (1, 1e-5, 1e-4))
        m = refusal(() -> save_stride(args...))
        @test occursin("a save is at least one reduced step", m)
        @test !occursin(phrase, m)
    end
end
