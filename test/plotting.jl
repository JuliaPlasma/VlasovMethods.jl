using Test
using TOML
using VlasovMethods

# `plot_distribution` and `plot_particles` are declared in `src/plotting.jl` and every method of
# them is in `ext/VlasovMethodsMakieExt.jl`. This file asserts the half of that split a caller
# without CairoMakie has: neither name has a method, and each call raises the plain `MethodError`
# of a function with no methods.
#
# That half needs a fresh environment, because the environment this file runs in holds CairoMakie,
# which is a test dependency. The check therefore runs in a subprocess whose own manifest has no
# Makie UUID, and this file claims nothing about the extension in the process it runs in — a
# runner may have loaded CairoMakie in an earlier file.
#
# The rendering tests are in `test/integration/makie_extension.jl`, which loads CairoMakie.

@testset "a fresh environment with no Makie loads the package without the extension" begin
    program = """
        import Pkg
        Pkg.activate(; temp = true, io = devnull)
        Pkg.develop(path = $(repr(pkgdir(VlasovMethods))); io = devnull)

        # a weak dependency of the package is not in the manifest of an environment that names
        # neither it nor a package that depends on it, so this environment has no Makie at all
        dependencies = Pkg.dependencies()
        @assert !any(occursin("Makie", info.name) for info in values(dependencies)) "a Makie package is in the environment"
        @assert !haskey(dependencies, Base.UUID("ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a")) "Makie is in the environment"
        @assert !haskey(dependencies, Base.UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0")) "CairoMakie is in the environment"

        using VlasovMethods

        @assert Base.get_extension(VlasovMethods, :VlasovMethodsMakieExt) === nothing
        @assert isempty(methods(VlasovMethods.plot_distribution))
        @assert isempty(methods(VlasovMethods.plot_particles))
        @assert isempty(methods(VlasovMethods.plot_errors))
        @assert isempty(methods(VlasovMethods.plot_energy))
        @assert isempty(methods(VlasovMethods.plot_growth_rate))

        grid = VlasovMethods.GridDistribution(4, 3, (0.0, 1.0), (-1.0, 1.0))
        particles = VlasovMethods.ParticleDistribution(1, 1, 1)
        for call in (
            () -> VlasovMethods.plot_distribution(grid),
            () -> VlasovMethods.plot_distribution(particles; xdomain = (0.0, 1.0), vdomain = (-1.0, 1.0)),
            () -> VlasovMethods.plot_particles(particles; xdomain = (0.0, 1.0), vdomain = (-1.0, 1.0)),
            () -> VlasovMethods.plot_errors([0.0, 1.0], ([0.5, 0.5],)),
            () -> VlasovMethods.plot_energy([0.0, 1.0], ([0.5, 0.5],)),
            () -> VlasovMethods.plot_growth_rate([0.0, 1.0], ([0.5, 0.5],)))
            exception = try
                call()
                nothing
            catch thrown
                thrown
            end
            @assert exception isa MethodError "an unloaded plotting call threw \$(typeof(exception))"
        end

        print("no Makie in the environment: the extension is not loaded and neither name has a method")
        """

    # `JULIA_LOAD_PATH` without `@v1.x`, so that a CairoMakie of the shared environment cannot
    # satisfy a stray `using` in the program; its entries are separated by `;` on Windows
    load_path = join(("@", "@stdlib"), Sys.iswindows() ? ';' : ':')
    command = addenv(`$(Base.julia_cmd()) --startup-file=no -e $program`,
        "JULIA_LOAD_PATH" => load_path)

    # `success` reports only the exit status, so the subprocess's output is shown on a failure
    output = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = output, stderr = output))
    success(process) || print(String(take!(output)))
    @test success(process)
end

@testset "the docstrings name the extension and CairoMakie" begin
    @test occursin("package extension", string(@doc VlasovMethods.plot_distribution))
    @test occursin("load `CairoMakie`", string(@doc VlasovMethods.plot_distribution))
    @test occursin("package extension", string(@doc VlasovMethods.plot_particles))
    @test occursin("load `CairoMakie`", string(@doc VlasovMethods.plot_particles))
end

@testset "the three series plots name the extension and CairoMakie" begin
    # `@doc` reads the name literally, so each of the three is spelled out
    for doc in (string(@doc VlasovMethods.plot_errors),
        string(@doc VlasovMethods.plot_energy),
        string(@doc VlasovMethods.plot_growth_rate))
        @test occursin("package extension", doc)
        @test occursin("load `CairoMakie`", doc)
    end
end

@testset "CairoMakie is a weak dependency" begin
    project = TOML.parsefile(joinpath(pkgdir(VlasovMethods), "Project.toml"))

    @test get(project["weakdeps"], "CairoMakie", nothing) ==
          "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
    @test get(project["extensions"], "VlasovMethodsMakieExt", nothing) == "CairoMakie"
    @test haskey(project["compat"], "CairoMakie")

    # No Makie of any kind is a dependency, so a fresh environment that holds VlasovMethods and
    # nothing else has no Makie UUID in its manifest, and the extension stays unloaded there.
    @test isempty([name for name in keys(project["deps"]) if occursin("Makie", name)])

    # Legacy plotting stacks and unrelated modelling packages are not hard dependencies.
    for name in ("Plots", "StatsPlots", "LaTeXStrings", "Optimisers", "DifferentialEquations")
        @test !haskey(project["deps"], name)
    end
end
