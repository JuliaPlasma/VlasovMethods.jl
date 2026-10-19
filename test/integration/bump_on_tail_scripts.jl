# Check 5 of §5V.P38b: the file the three bump-on-tail scripts share, the field it hands a march,
# and the scale that march applies. The scripts are not a package, so the file is included rather
# than imported; it uses VlasovMethods and the standard library only, which is what lets this file
# include it.

using PoissonSolvers: PeriodicBasisSpline, Potential
using Random
using Test
using VlasovMethods
using VlasovMethods: PoissonField, coefficients
using ParticleMethods: ParticleList
using ReducedComplexityModeling: Parameter, ParameterSpace

include(joinpath(@__DIR__, "..", "..", "scripts", "bump_on_tail_common.jl"))

"""
    script_fixture(::Type{T}; n = 64, seed = 7, ncells = 16) -> (plt, field)

The bump-on-tail fixture of the check: `n` particles drawn under `seed` on the scripts' domain,
and the field that `vp_field` builds from a periodic spline potential there.
"""
function script_fixture(::Type{T}; n = 64, seed = 7, ncells = 16) where {T}
    Random.seed!(seed)
    x = T(domain_length()) .* rand(T, n)
    v = sqrt(T(2)) .* randn(T, n)
    for i in 1:n
        rand() < α_BOT && (v[i] = σ_BOT * v[i] + v₀_BOT)
    end
    w = fill(one(T) / n, n)
    plt = ParticleList(reshape(x, 1, n), reshape(v, 1, n), reshape(w, 1, n))
    basis = PeriodicBasisSpline((zero(T), T(domain_length())), DEGREE_BOT, ncells)
    return plt, vp_field(Potential(basis))
end

"One full-order march of `plt` on `field`, `nt` steps of `dt` at the scale `χ`, every step saved."
function script_march(::Type{T}, plt, field, nt; dt, χ) where {T}
    n = size(plt.x, 2)
    nh = length(coefficients(field))
    space = ParameterSpace(Parameter(:χ, T(χ), T(χ), 1))
    SS = Snapshots(T, 1, n, nh, nt + 1, space)
    IP = IntegratorParameters(dt, nt, nt + 1, nh, n, 1)
    integrate_vp(plt, march_parameters(χ), field, SS, IP, 1)
    return SS
end

@testset "the sweeps of the bump-on-tail scripts" begin
    # the training sweep is the paper's ten κ over the range, each over the reference wave number
    χ = training_chi()
    @test length(χ) == NSAMPLES_BOT == 10
    @test issorted(χ)
    @test first(χ) == κ_RANGE[1] / κ_REF
    @test last(χ) == κ_RANGE[2] / κ_REF
    @test χ == training_chi()

    # the test sweep is ten more, sorted and inside the range, and the seed reproduces it
    χt = test_chi()
    @test length(χt) == NSAMPLES_BOT
    @test issorted(χt)
    @test all(κ_RANGE[1] / κ_REF .≤ χt .≤ κ_RANGE[2] / κ_REF)
    @test χt == test_chi()
end

@testset "vp_field is the unscaled field" begin
    @test domain_length() == 2π / 0.3
    basis = PeriodicBasisSpline((0.0, domain_length()), DEGREE_BOT, NCELLS_BOT)
    potential = Potential(basis)
    field = vp_field(potential)
    @test field isa PoissonField
    @test field.potential === potential
    @test length(coefficients(field)) == nbasis(basis)
end

@testset "the initial condition of the set is seeded" begin
    # the same seed draws the same particles twice, and the second draw of a running generator
    # would not: the scripts' sample is reproducible, and the run carries its seed as an attribute
    first = initial_distribution(; n = 16)
    second = initial_distribution(; n = 16)
    @test first.particles.x == second.particles.x
    @test first.particles.v == second.particles.v
    @test first.particles.w == second.particles.w
    @test !isempty(commit_id())
end

@testset "the march the scripts feed applies params.χ once" begin
    T = Float64
    n, nt, dt, χ = 64, 20, 1e-3, 2
    plt, field = script_fixture(T; n = n)

    # `(dt, χ = 2)` is `(2dt, χ = 1)` state for state: the χ reaches the timestep and the field
    # once, and the field the scripts build carries no scale of its own
    scaled = script_march(T, plt, field, nt; dt = dt, χ = χ)
    plain = script_march(T, plt, field, nt; dt = 2dt, χ = 1)
    @test scaled.X == plain.X
    @test scaled.V == plain.V
    @test scaled.Φ == plain.Φ
    @test scaled.W == plain.W
    @test scaled.K == plain.K
    @test scaled.M == plain.M

    # a `ScaledField` at exactly `params.χ` keeps its own χ² divisor, applied once: its initial
    # save is the plain march's divided by χ², not divided twice
    at = script_march(T, plt, ScaledField(field, χ), 1; dt = dt, χ = χ)
    @test at.W[1, 1] == scaled.W[1, 1] / T(χ)^2
    @test at.W[1, 1] != scaled.W[1, 1]

    # a stored χ that differs is refused before the loop
    @test_throws ArgumentError script_march(
        T, plt, ScaledField(field, 5), nt; dt = dt, χ = χ)
end

@testset "no script builds a ScaledField" begin
    root = pkgdir(VlasovMethods)
    @test root !== nothing
    @test !success(`git -C $root grep -q ScaledField -- scripts`)
end
