# The four drivers in `scripts/`, each run in a new Julia process of this environment at small
# options: training, projections and testing with the same options, in the order of the docs page,
# and the analytic model. Each writes into new directories, and the tests read what it wrote. The
# drivers' functions are tested in `bump_on_tail_scripts.jl`.

using HDF5
using Test
using VlasovMethods

const REPOSITORY = joinpath(@__DIR__, "..", "..")

"""
    driver(script, args...) -> (exit code, output)

Run the driver `script` of `scripts/` with the arguments `args` in a new Julia process of the
active environment, and return its exit code and its output.
"""
function driver(script, args...)
    output = IOBuffer()
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())
               $(joinpath(REPOSITORY, "scripts", script)) $args`
    process = run(pipeline(ignorestatus(command); stdout = output, stderr = output))
    return process.exitcode, String(take!(output))
end

"Whether the driver ran to exit 0; its output is printed where it did not."
function succeeded((code, output))
    code == 0 || println(output)
    return code == 0
end

"The root attributes `commit`, `command` and `seed` of the HDF5 file `file`, those that it has."
function provenance(file)
    h5open(file) do io
        Dict(k => read(attributes(io)[k])
        for k in ("commit", "command", "seed")
        if haskey(attributes(io), k))
    end
end

"The PNG files of `dir` that start with `prefix`, each of `filesize > 0`."
function pngs(dir, prefix)
    filter(
        f -> startswith(f, prefix) && endswith(f, ".png") &&
             filesize(joinpath(dir, f)) > 0, readdir(dir))
end

const COMMIT = readchomp(`git -C $REPOSITORY rev-parse HEAD`)
const STATUS = read(`git -C $REPOSITORY status --porcelain --ignored`, String)

mktempdir() do dir
    runs, results = joinpath(dir, "runs"), joinpath(dir, "results")
    options = ("--particles", "1000", "--final-time", "10.0", "--seed", "7",
        "--runs-dir", runs, "--results-dir", results)

    @testset "training" begin
        @test succeeded(driver("bump_on_tail_training.jl", options...))
        file = only(f for f in readdir(runs; join = true) if endswith(f, ".h5"))
        @test provenance(file) == Dict("commit" => COMMIT,
            "command" => "bump_on_tail_training.jl " * join(options, " "), "seed" => 7)

        # every sample of the sweep has its own accelerations
        SS = h5open(io -> Snapshots(io), file, "r")
        @test size(SS.A, 4) == 10
        @test all(p -> any(!iszero, SS.A[:, :, :, p]), axes(SS.A, 4))
        @test all(p -> SS.A[:, :, :, p] != SS.A[:, :, :, p + 1], 1:9)
    end

    @testset "projections" begin
        @test succeeded(driver("bump_on_tail_projections.jl", options...))
        projections = only(f
        for f in readdir(runs; join = true)
        if endswith(f, "_projections.h5"))
        vectorfield = only(f
        for f in readdir(runs; join = true)
        if endswith(f, "_vectorfield.h5"))
        for f in (projections, vectorfield)
            @test provenance(f)["commit"] == COMMIT
            @test startswith(provenance(f)["command"], "bump_on_tail_projections.jl --particles")
            @test provenance(f)["seed"] == 7
        end
        X, A = h5open(io -> (read(io["X"]), read(io["A"])), vectorfield)
        @test size(X) == size(A)
        @test size(X, 3) == 10
        @test all(p -> any(!iszero, A[:, :, p]), axes(A, 3))
        @test length(pngs(results, "item11_")) == 2
    end

    @testset "testing" begin
        @test succeeded(driver("bump_on_tail_testing.jl", options..., "--samples", "3"))
        rates = readlines(joinpath(results, "growth_rates.csv"))
        @test rates[1] == "kappa,gamma_hfm,gamma_rom"
        @test length(rates) == 4
        @test length(pngs(results, "item24_")) == 8
        @test length(pngs(results, "item25_")) == 1
        @test "item26_field_energy.png" in pngs(results, "item26_")
        @test length(pngs(results, "item27_")) == 3
        @test length(pngs(results, "item28_")) == 1
        @test length(pngs(results, "item29_")) == 3
    end

    @testset "analytic" begin
        args = ("--final-time", "1.0", "--saves", "10", "--particles", "100",
            "--runs-dir", runs, "--results-dir", results)
        @test succeeded(driver("analytic_quadratic.jl", args...))
        for name in ("analytic_full", "analytic_ode_k5", "analytic_pode_k20")
            file = joinpath(runs, name * ".h5")
            @test provenance(file) == Dict("commit" => COMMIT,
                "command" => "analytic_quadratic.jl " * join(args, " "))
            h5open(file) do io
                @test size(read(io["q"])) == (100, 11)
                @test all(isfinite, read(io["H"]))
            end
        end
        @test length(pngs(results, "item11_analytic_")) == 2
        @test length(pngs(results, "item24_analytic_")) == 2
        @test length(pngs(results, "item27_analytic_")) == 1

        # a run with no save is refused before it computes anything
        code, output = driver("analytic_quadratic.jl", "--saves", "0",
            "--runs-dir", runs, "--results-dir", results)
        @test code != 0
        @test occursin("a run saves at least once", output)
    end
end

@testset "the checkout is unchanged" begin
    @test read(`git -C $REPOSITORY status --porcelain --ignored`, String) == STATUS
end
