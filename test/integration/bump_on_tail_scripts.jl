using HDF5
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

# The exception a call raises, or `nothing` when it raises none: a test of a refusal names the
# phrase of the guard and not only the type of the exception.
function caught(f)
    try
        f()
    catch thrown
        return thrown
    end
    return nothing
end

# The parameter space of the march, which holds one sample.
space(::Type{T}) where {T} = ParameterSpace(Parameter(:χ, T(1), T(1), 1))

# A `git` that answers the commit of the checkout, so that a script runs where the tree is no
# checkout. The commit itself is read by the test of `repository_commit`, not by the scripts.
function stub_git()
    directory = mktempdir()
    path = joinpath(directory, "git")
    write(path, "#!/bin/sh\necho 0000000000000000000000000000000000000000\n")
    chmod(path, 0o755)
    return directory
end

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
    # the floating-point `div` of `25.0 / 0.1` gives 249, and `round` gives 250
    @test step_count(25.0, 0.1) == 250
    @test step_count(10.0, 1e-3) == 10000
    @test step_count(1.0, 0.25) == 4
end

@testset "the stride of the saved times is a whole number of reduced steps" begin
    @test save_stride(10, 1e-3, 1e-3) == 10
    @test save_stride(10, 1e-3, 1e-4) == 100

    # a ratio that is a little below a whole number is the whole number
    @test save_stride(3, 0.7, 0.1) == 21

    # a ratio that is a little below a whole number is the whole number, when the difference is
    # the round-off of the four roundings of the quotient
    @test save_stride(3, 0.7, 0.084) == 25

    # the bound of the guard is the round-off of the four roundings: at a stride of two the bound
    # is two units in the last place, a ratio that lies exactly at it is the whole number, and a
    # ratio that lies one unit beyond it is not
    @test save_stride(1, 2.0 - 2.0^-50, 1.0) == 2
    thrown = caught(() -> save_stride(1, 2.0 - 3 * 2.0^-51, 1.0))
    @test thrown isa ArgumentError
    @test occursin("which is not a whole number", sprint(showerror, thrown))

    # a ratio that is not a whole number is refused: the reduced states would not be the states
    # of the saved times
    for (save_step, dt, rdt) in ((10, 1e-3, 3e-3), (10, 1e-3, 9e-4),
        (1, 100.00000001, 1.0), (10, 1e-3, 1.0000000001e-4),
        (1, 1e9 + 0.25, 1.0), (1, 1e9 + 0.5, 1.0))
        thrown = caught(() -> save_stride(save_step, dt, rdt))
        @test thrown isa ArgumentError
        @test occursin("which is not a whole number", sprint(showerror, thrown))
    end
end

@testset "the analytic script refuses a stride that is not a whole number" begin
    # The script is run in its own process, because it is a script: it refuses before it marches,
    # and it writes its files into the directories of the run.
    script = joinpath(@__DIR__, "..", "..", "scripts", "analytic_quadratic.jl")
    mktempdir() do dir
        command = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())
            $script --final-time 0.03 --rdt 0.003 --runs-dir $dir/runs
            --results-dir $dir/results`
        err = IOBuffer()
        command = setenv(command, "PATH" => stub_git() * ":" * ENV["PATH"])
        process = run(pipeline(ignorestatus(command), stderr = err), wait = false)
        wait(process)
        @test process.exitcode != 0
        @test occursin("which is not a whole number", String(take!(err)))
    end
end

@testset "the training script fills the accelerations of the saved states" begin
    T = Float64
    potential, plt, x, w = fixture(T)
    params = (χ = 2.0,)
    field = vp_field(potential, params)
    nh = length(coefficients(field))

    # two samples, so that the fill of one of them is seen
    samples = ParameterSpace(Parameter(:χ, T(1), T(2), [T(1), T(2)]))
    SS = Snapshots(T, 1, length(x), nh, 3, samples)
    ip = IntegratorParameters(T(1e-3), 2, 3, nh, length(x), 2)
    integrate_vp(plt, params, field, SS, ip, 1)
    integrate_vp(plt, params, field, SS, ip, 2)

    # the accelerations of the saved states, which the march does not save, are the field of the
    # sample at the positions of the state, divided by χ²
    fill_accelerations!(SS, field, ip, vec(plt.w), 2)

    plain = PoissonField(potential)
    update!(plain, x, w, zero(T))
    e = similar(x)
    efield!(plain, e, x)
    @test !all(iszero, e)

    # the second sample alone carries them
    @test all(iszero, SS.A[1, :, 1, 1])
    @test SS.A[1, :, 1, 2] == e ./ 4

    # and every saved state carries them, the last one included
    xlast = SS.X[1, :, 3, 2]
    update!(plain, xlast, w, zero(T))
    elast = similar(x)
    efield!(plain, elast, xlast)
    @test SS.A[1, :, 3, 2] == elast ./ 4
end

@testset "the training run holds the accelerations of its saved states" begin
    # The script is run in its own process, because it is a script: it marches the sweep of the
    # wave number and writes the snapshots of the run.
    script = joinpath(@__DIR__, "..", "..", "scripts", "bump_on_tail_training.jl")
    mktempdir() do dir
        command = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())
            $script --npart 32 --final-time 0.4 --runs-dir $dir/runs
            --results-dir $dir/results`
        command = setenv(command, "PATH" => stub_git() * ":" * ENV["PATH"])
        process = run(pipeline(ignorestatus(command), stderr = devnull), wait = false)
        wait(process)
        @test process.exitcode == 0

        file = joinpath(dir, "runs", "BoT_Np32_k_010_050_ns10_T0.4.h5")
        @test isfile(file)
        snapshots = h5open(file, "r") do handle
            Snapshots(handle)
        end
        @test !all(iszero, snapshots.A)
    end
end

@testset "the name of a run carries its parameters" begin
    @test bump_on_tail_runid(1000, 10, 25.0) != bump_on_tail_runid(2000, 10, 25.0)
    @test bump_on_tail_runid(1000, 10, 25.0) != bump_on_tail_runid(1000, 20, 25.0)
    @test bump_on_tail_runid(1000, 10, 25.0) != bump_on_tail_runid(1000, 10, 50.0)

    # the defaults of the three scripts are the paper's parameters, and the run they name is the
    # paper's run
    @test DEFAULT_NPART == 10000
    @test DEFAULT_NSAMPLES == 10
    @test DEFAULT_FINAL_TIME == 25.0
    @test DEFAULT_RUNID == "BoT_Np10000_k_010_050_ns10_T25"
end

@testset "the commit of the checkout" begin
    # a repository of its own, so that the commit is read
    mktempdir() do dir
        run(`git -C $dir init --quiet`)
        run(`git -C $dir -c user.name=test -c user.email=test@example.com commit --quiet
            --allow-empty -m "the initial commit"`)
        @test length(repository_commit(dir)) == 40
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
