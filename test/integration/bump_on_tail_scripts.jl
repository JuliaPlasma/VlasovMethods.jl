# Check 5 of §5V.P38b: the file the three bump-on-tail scripts share, the field it hands a march,
# and the scale that march applies. The scripts are not a package, so the file is included rather
# than imported; it uses VlasovMethods and the standard library only, which is what lets this file
# include it.

using PoissonSolvers: PeriodicBasisSpline, Potential
using Random
using Test
using VlasovMethods
using VlasovMethods: PoissonField, ScaledField, coefficients, efield!, update!
using ParticleMethods: ParticleList
using ReducedComplexityModeling: Parameter, ParameterSpace

include(joinpath(@__DIR__, "..", "..", "scripts", "bump_on_tail_common.jl"))

"""
    script_fixture(::Type{T}; n = 64, seed = 7, ncells = 16) -> (plt, potential)

The bump-on-tail fixture of the check: `n` particles drawn under `seed` on the scripts' domain,
and the periodic spline potential that `vp_field` takes there.
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
    return plt, Potential(basis)
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

@testset "vp_field is the sample's scaled field" begin
    @test domain_length() == 2π / 0.3
    basis = PeriodicBasisSpline((0.0, domain_length()), DEGREE_BOT, NCELLS_BOT)
    potential = Potential(basis)
    params = march_parameters(2.0)

    # the field of a sample is its scaled field, over the potential and at the χ of the same
    # parameter tuple the march reads
    field = vp_field(potential, params)
    @test field isa ScaledField
    @test field.field isa PoissonField
    @test field.field.potential === potential
    @test field.χ === params.χ
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
end

@testset "the commit attribute is the checkout's HEAD" begin
    # a checkout is read, and a directory that is not one falls back to `"unknown"`, which is what
    # the attribute of a file written outside a checkout holds
    root = mktempdir()
    run(`git init --quiet $root`)
    run(`git -C $root -c user.email=scripts@example.com -c user.name=scripts commit --quiet --allow-empty -m init`)
    sha = readchomp(`git -C $root rev-parse HEAD`)
    @test occursin(r"^[0-9a-f]{40}$", sha)
    @test commit_id(root) == sha
    @test commit_id(mktempdir()) == "unknown"
end

@testset "the scale of the scripts' field reaches the march once" begin
    T = Float64
    n, nt, dt, χ = 64, 20, 1e-3, 2
    plt, potential = script_fixture(T; n = n)
    params = march_parameters(T(χ))
    field = vp_field(potential, params)
    plain = PoissonField(potential)

    # the wrapper's field is the plain field's divided by χ², and its deposit is the plain deposit
    x, w = vec(plt.x), vec(plt.w)
    e_plain, e_scaled = zeros(T, n), zeros(T, n)
    update!(plain, x, w, 0.0)
    efield!(plain, e_plain, x)
    update!(field, x, w, 0.0)
    efield!(field, e_scaled, x)
    @test e_scaled == e_plain ./ T(χ)^2

    # the march the scripts run applies the scale once: its first save is the plain march's / χ²
    scaled = script_march(T, plt, field, nt; dt = dt, χ = χ)
    at1 = script_march(T, plt, plain, 1; dt = dt, χ = χ)
    @test scaled.W[1, 1] == at1.W[1, 1] / T(χ)^2

    # and it is not a change of the time unit alone: the `(dt, χ = 2)` march of the scaled field
    # differs from the `(2dt, χ = 1)` march of the plain one, which carries no field scale
    retimed = script_march(T, plt, plain, nt; dt = 2dt, χ = 1)
    @test scaled.X != retimed.X

    # a stored χ other than the parameters' is refused before the loop
    @test_throws ArgumentError script_march(
        T, plt, ScaledField(plain, 5), nt; dt = dt, χ = χ)
end

@testset "one script builds the scale" begin
    # read the directory rather than `git grep`: a mutated copy of the tree carries no `.git`, and
    # a git command there fails whatever the scripts hold
    root = pkgdir(VlasovMethods)
    @test root !== nothing
    scripts = joinpath(root, "scripts")
    @test isdir(scripts)
    files = sort(filter(endswith(".jl"), readdir(scripts)))
    for name in ("bump_on_tail_common.jl", "bump_on_tail_training.jl",
        "bump_on_tail_projections.jl", "bump_on_tail_testing.jl", "analytic_quadratic.jl")
        @test name in files
    end

    builds = [(name, count("ScaledField(", read(joinpath(scripts, name), String)))
              for name in files]
    @test sum(last, builds) == 1
    @test only(filter(!iszero ∘ last, builds)) == ("bump_on_tail_common.jl", 1)
end
