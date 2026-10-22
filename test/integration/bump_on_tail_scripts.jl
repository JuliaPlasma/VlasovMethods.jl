# The tests of the four bump-on-tail scripts and of the code they share: the field every march
# takes, the two sweeps, the growth-rate fit of item 28, the P39 writer, and the docs page. It
# includes the scripts' common file, which loads no Makie package, and it runs no script: the four
# defaults run from a clean checkout by hand (T8.8).
#
# The file reads `scripts/*.jl`, so a name or a path that a scenario scans for is asserted here
# rather than by a one-time command.

using HDF5
using LinearAlgebra
using Random
using ReducedBasisMethods: CotangentLiftEVD, PODEVD, ReducedBasis
using ReducedComplexityModeling: Parameter, ParameterSpace
using Test
using VlasovMethods
using VlasovMethods: efield!, update!

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SCRIPTS = joinpath(ROOT, "scripts")

# whether Makie is loaded before the include: another testset of the same process may have loaded
# it, and what the clause asks is that the include does not
makie_loaded() = any(k -> k.name == "CairoMakie", keys(Base.loaded_modules))
const MAKIE_BEFORE = makie_loaded()

include(joinpath(SCRIPTS, "bump_on_tail_common.jl"))

# the scripts of the paper's two example sets, in the order a run takes them
const PAPER_SCRIPTS = ("bump_on_tail_training.jl", "bump_on_tail_projections.jl",
    "bump_on_tail_testing.jl", "analytic_quadratic.jl")

"The source of every `.jl` file of the scripts directory, by name."
function script_sources()
    Dict(f => read(joinpath(SCRIPTS, f), String)
    for f in readdir(SCRIPTS) if endswith(f, ".jl"))
end

@testset "the common file" begin
    # the include above is the test of a load: a package that `test/Project.toml` does not hold, or
    # Makie, fails here rather than in a run of a script
    @test isdefined(@__MODULE__, :vp_field)
    @test isdefined(@__MODULE__, :training_sweep)
    @test makie_loaded() == MAKIE_BEFORE
    # and no line of the common file loads one, whatever the process has already loaded
    common = read(joinpath(SCRIPTS, "bump_on_tail_common.jl"), String)
    @test !occursin(r"(?m)^\s*(?:using|import)\s+[^\n]*Makie", common)

    for f in PAPER_SCRIPTS
        @test isfile(joinpath(SCRIPTS, f))
    end
end

@testset "vp_field" begin
    np = 2000
    dist = bump_on_tail_particles(np)
    potential = bot_potential(16, 3)
    x, w = vec(dist.particles.x), vec(dist.particles.w)
    params = (χ = 2.0,)

    field = vp_field(potential, params)
    @test field isa ScaledField
    @test field.field isa PoissonField
    @test field.χ == 2.0

    # the field of a solved deposit at `χ = 2` is one quarter of the plain field's, exactly
    update!(field, x, w, 0.0)
    plain = PoissonField(potential)
    scaled = efield!(field, zeros(np), x)
    unscaled = efield!(plain, zeros(np), x)
    @test any(!iszero, unscaled)
    @test scaled == unscaled ./ 4

    # one construction of a `ScaledField` in `scripts/`, and it is the common file's
    sources = script_sources()
    counts = Dict(f => count("ScaledField(", s) for (f, s) in sources)
    @test sum(values(counts)) == 1
    @test only([f for (f, n) in counts if n > 0]) == "bump_on_tail_common.jl"
end

@testset "paths" begin
    # every file a script writes goes through `runs_path` or `results_path`: no script names the
    # reference's own output directory, nor a run directory above the checkout
    for (f, source) in script_sources()
        @test !occursin("../runs", source)
        @test !occursin("/toks/", source)
    end
end

@testset "sweeps" begin
    # `round(Int, T/dt)` and not `Int(div(T, dt))`: the truncated count drops the final state of
    # a span whose quotient is not an integer
    @test step_count(25.0, 0.1) == 250
    @test step_count(25.0, 0.5) == 50

    a = test_sweep(1234)
    b = test_sweep(1)
    @test length(a.κ) == 10
    @test issorted(a.κ)
    @test all(κ -> 0.1 ≤ κ ≤ 0.5, a.κ)
    @test a.χ == a.κ ./ 0.3
    @test b.χ == b.κ ./ 0.3
    # the seed decides the sweep, and another seed is another sweep
    @test a.κ != b.κ
    @test test_sweep(1234).κ == a.κ

    # the training sweep is the ten values of the paper, explicitly
    @test length(TRAINING_KAPPA) == 10
    @test first(TRAINING_KAPPA) == 0.1 && last(TRAINING_KAPPA) == 0.5
    @test issorted(TRAINING_KAPPA)
end

@testset "growth rate" begin
    # the reference's `W` of cell 36: an exponential of `γ = 0.05` that is flat after its fourth
    # maximum, on a window of five samples
    t = collect(0:0.001:20.0)
    W = exp.(0.1 .* min.(t, 3π)) .* cos.(t) .^ 2
    @test length(VlasovMethods._find_maxima(W, 5)) == 7
    @test isapprox(fit_growth_rate(t, W), 0.05; rtol = 1e-2)
    # `growth_rate`'s own defaults are not the reference's fit: `window = 2` and `maxima = 2:5`
    # put a saturated maximum into the regression (measured `0.0350` on the base)
    @test growth_rate(t, W) < 0.045

    # three maxima at the window: every one of them enters the fit
    t3 = collect(0:0.001:2.5π)
    W3 = exp.(0.1 .* min.(t3, 3π)) .* cos.(t3) .^ 2
    @test length(VlasovMethods._find_maxima(W3, 5)) == 3
    @test isapprox(fit_growth_rate(t3, W3), 0.05; rtol = 1e-2)

    # `log W` is linear through the fourth maximum of the history above, so `1:3` and `1:4` fit the
    # same line there; a history that saturates before its fourth maximum separates them, and the
    # fit is the one over the four
    ts = collect(0:0.001:20.0)
    Ws = exp.(0.1 .* min.(ts, 3.0)) .* cos.(ts) .^ 2
    @test length(VlasovMethods._find_maxima(Ws, 5)) == 7
    @test fit_growth_rate(ts, Ws) == growth_rate(ts, Ws; window = 5, maxima = 1:4)
    @test fit_growth_rate(ts, Ws) != growth_rate(ts, Ws; window = 5, maxima = 1:3)

    # a history with fewer than two maxima has no fit
    @test_throws ArgumentError fit_growth_rate(t, ones(length(t)))
end

@testset "p39 files" begin
    dir = mktempdir()
    space = ParameterSpace(Parameter(:χ, [1.0, 2.0]))
    snap = randn(Xoshiro(7), 8, 6)
    particle = ReducedBasis(CotangentLiftEVD(), snap, (α = 0.1,), space; k = 2)
    field = ReducedBasis(PODEVD(), snap, (α = 0.1,), space; k = 2)
    k, nₛ, nparam = particle.k, 3, 2
    X = randn(Xoshiro(8), k, nₛ, nparam)
    A = randn(Xoshiro(9), k, nₛ, nparam)
    commit = repeat("0123456789", 4)

    first = (np = 1000, nh = 16, order = 3, dt = 0.1,
        T = 25.0, κlo = 0.1, κhi = 0.5, seed = 1234)
    second = merge(first, (np = 2000,))
    @test run_id(first) != run_id(second)

    # a seed that `BOT_SEED` does not hold, so that a writer which hardcodes the constant fails the
    # attribute check below
    seed = 7
    ppath, vpath = write_p39_files(dir, first, particle, field, X, A;
        commit = commit, command = "--np 1000", seed = seed)
    write_p39_files(dir, second, particle, field, X, A;
        commit = commit, command = "--np 2000", seed = seed)

    @test isfile(ppath) && isfile(vpath)
    @test length(readdir(dir)) == 4
    @test basename(ppath) == "$(run_id(first))_projections.h5"
    @test basename(vpath) == "$(run_id(first))_vectorfield.h5"

    # the three root attributes, on both files, and the commit is the caller's
    for path in (ppath, vpath)
        h5open(path) do io
            attrs = attributes(io)
            @test haskey(attrs, "commit") && read(attrs["commit"]) == commit
            @test haskey(attrs, "command") && read(attrs["command"]) isa String
            @test haskey(attrs, "seed") && read(attrs["seed"]) == seed
        end
    end

    # the reduced data is the size the reader of P39 expects, and the basis reads back
    h5open(vpath) do io
        @test size(read(io["X"])) == (k, nₛ, nparam)
        @test size(read(io["A"])) == size(read(io["X"]))
    end
    @test ReducedBasis(ppath, "particle").k == k
    @test ReducedBasis(ppath, "field").k == field.k

    # a reduced state of the wrong shape is refused
    @test_throws ArgumentError write_p39_files(
        dir, first, particle, field, X, zeros(k, nₛ, 1);
        commit = commit, command = "", seed = 1234)
    @test_throws ArgumentError write_p39_files(dir, first, particle, field,
        zeros(k + 1, nₛ, nparam), zeros(k + 1, nₛ, nparam); commit = commit, command = "",
        seed = 1234)
end

@testset "docs page" begin
    make = read(joinpath(ROOT, "docs", "make.jl"), String)
    @test occursin("\"paper.md\"", make)

    page = read(joinpath(ROOT, "docs", "src", "paper.md"), String)
    for script in PAPER_SCRIPTS
        @test occursin(script, page)

        # each script's own section gives its paper command, its run time in seconds and a file it
        # writes. The page holds the path with `/` on every platform, so the path is built in the
        # page's form rather than from the filesystem, whose separator is `\` on Windows
        path = "scripts/" * script
        @test occursin("### `$path`", page)
        section = last(split(page, "### `$path`"))
        @test occursin("julia --project=scripts " * path, section)
        @test occursin(r"\d+ s", section)
        @test occursin(r"\.(h5|csv|png)", section)
    end
    # the page describes the runs and embeds nothing: a docs build runs no script and draws no
    # figure
    @test !occursin("@example", page)
    @test !occursin("@eval", page)
    @test !occursin("![", page)
end

@testset "figures" begin
    src = read(joinpath(SCRIPTS, "bump_on_tail_testing.jl"), String)

    # item 24 draws the four components of `trajectory_errors` that D-P38-1 names: the root mean
    # squares `x2` and `v2` and the `L∞` maxima `xmax` and `vmax`, each against time and against κ
    # at `T`, `0.75T` and `0.5T`
    chunks = split(src, "figure = plot_series")
    @test length(chunks) > 1
    plotting = join(chunks[2:end])
    for component in ("x2", "v2", "xmax", "vmax")
        @test occursin(component, plotting)
    end

    # items 26 and 27 plot the paper's `W`, the energy of the plain field: a march stores the
    # energy of the scaled field, `W/χ²`, and the script multiplies it back
    @test occursin("HFM.W[:, 1] .* χ^2", src)
    @test occursin("RSS.W[:, 1] .* χ^2", src)
end
