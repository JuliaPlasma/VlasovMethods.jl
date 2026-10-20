using ParticleMethods: ParticleList
using PoissonSolvers
using ReducedComplexityModeling: Parameter, ParameterSpace
using Test
using VlasovMethods
using VlasovMethods: PoissonField, ScaledField, coefficients, efield!, update!

include("../../scripts/bump_on_tail_common.jl")

# The potential of the bump-on-tail model, the particles of a small set and their weights. The
# positions are seeded, so that a march of them is the same on every run.
function fixture(::Type{T}; n = 32, ncells = 8) where {T}
    potential = Potential(PeriodicBasisSpline((zero(T), T(2π) / κ_ref), 3, ncells))
    x = T.(range(0.1, 6.0; length = n))
    v = T(0.5) .* x
    w = fill(T(1) / n, n)
    plt = ParticleList(reshape(x, 1, n), reshape(v, 1, n), reshape(w, 1, n))
    return potential, plt, x, w
end

# The parameter space of the march, which holds one sample.
space(::Type{T}) where {T} = ParameterSpace(Parameter(:χ, T(1), T(1), 1))

# One full-order march of the particles on the field: the scaled field of the sample when
# `scaled` is set, the plain field otherwise.
function march(::Type{T}, plt, potential, x, w; dt = T(1e-3), χ = 1, nt = 10,
        scaled = false) where {T}
    field = scaled ? vp_field(potential, (χ = χ,)) : PoissonField(potential)
    nh = length(coefficients(field))
    SS = Snapshots(T, 1, length(x), nh, nt + 1, space(T))
    ip = IntegratorParameters(dt, nt, nt + 1, nh, length(x), 1)
    integrate_vp(plt, (χ = χ,), field, SS, ip, 1)
    return SS
end

@testset "the field of a sample carries the scale of its march" begin
    potential, _, _, _ = fixture(Float64)
    params = (χ = 2.0,)

    field = vp_field(potential, params)
    @test field isa ScaledField
    @test field.field isa PoissonField
    @test field.field.potential === potential
    @test field.χ === params.χ
end

@testset "the field of a sample divides the field by χ²" begin
    T = Float64
    potential, plt, x, w = fixture(T)

    plain = PoissonField(potential)
    scaled = vp_field(potential, (χ = 2.0,))
    update!(plain, x, w, zero(T))
    update!(scaled, x, w, zero(T))

    e = similar(x)
    e′ = similar(x)
    efield!(plain, e′, x)
    efield!(scaled, e, x)
    @test e == e′ ./ 4

    # the energy of the first saved state, which is the initial state in both marches: the scaled
    # one is a quarter of the plain one
    SSs = march(T, plt, potential, x, w; χ = 2.0, scaled = true)
    SSp = march(T, plt, potential, x, w; χ = 1)
    @test SSs.W[1, 1] == SSp.W[1, 1] / 4
end

@testset "the field of a sample changes the march" begin
    T = Float64
    potential, plt, x, w = fixture(T)

    # the scaled field divides the acceleration by four, so the march at `(dt, χ = 2)` is not the
    # plain march at `(2dt, χ = 1)`, which steps as long with the whole acceleration
    SSs = march(T, plt, potential, x, w; dt = T(1e-3), χ = 2.0, scaled = true)
    SSp = march(T, plt, potential, x, w; dt = T(2e-3), χ = 1)
    @test SSs.X[1, :, 2, 1] != SSp.X[1, :, 2, 1]
end

@testset "the initial condition of the set is seeded" begin
    first = initial_distribution(32)
    second = initial_distribution(32)

    @test first.particles.x == second.particles.x
    @test first.particles.v == second.particles.v
    @test first.particles.w == second.particles.w
    @test all(!iszero, first.particles.x)
    @test all(>(0), first.particles.w)
end

@testset "the scale of a sample is its wave number over the reference" begin
    @test march_parameters(κ_ref).χ == 1
    @test march_parameters(0.5).χ == 0.5 / κ_ref
    @test march_parameters(0.1).χ == 0.1 / κ_ref
end

@testset "the wave numbers of a sweep" begin
    κ = κ_sweep(0.1, 0.5, 10)
    @test length(κ) == 10
    @test issorted(κ)
    @test first(κ) == 0.1
    @test last(κ) == 0.5

    # the test sweep is drawn from the seeded generator and sorted
    κt = κ_test(0.1, 0.5, 10)
    @test length(κt) == 10
    @test issorted(κt)
    @test all(k -> 0.1 ≤ k ≤ 0.5, κt)
    @test κt == κ_test(0.1, 0.5, 10)
end

@testset "the step count of a march" begin
    # the old script's floating-point `div` gave 249 here
    @test step_count(25.0, 0.1) == 250
    @test step_count(10.0, 1e-3) == 10000
    @test step_count(1.0, 0.25) == 4
end

@testset "the name of a run carries its parameters" begin
    @test bump_on_tail_runid(1000, 10, 25.0) != bump_on_tail_runid(2000, 10, 25.0)
    @test bump_on_tail_runid(1000, 10, 25.0) != bump_on_tail_runid(1000, 20, 25.0)
    @test bump_on_tail_runid(1000, 10, 25.0) != bump_on_tail_runid(1000, 10, 50.0)
end

@testset "the commit of the checkout" begin
    # a repository of its own, so that the commit is read and not the fallback of a directory that
    # is no checkout
    mktempdir() do dir
        run(`git -C $dir init --quiet`)
        run(`git -C $dir -c user.name=test -c user.email=test@example.com commit --quiet
            --allow-empty -m "the initial commit"`)
        @test length(repository_commit(dir)) == 40
    end

    # a directory that is no checkout has no commit
    mktempdir() do dir
        @test repository_commit(dir) == "unknown"
    end
end

@testset "the scripts build the scaled field in one place" begin
    # the same check as `git grep -n 'ScaledField(' -- scripts`: one construction, in the file the
    # three scripts share
    directory = joinpath(@__DIR__, "..", "..", "scripts")
    hits = [(name, count("ScaledField(", read(joinpath(directory, name), String)))
            for name in readdir(directory) if endswith(name, ".jl")]

    @test sum(last, hits) == 1
    @test first(filter(h -> last(h) > 0, hits))[1] == "bump_on_tail_common.jl"
end
