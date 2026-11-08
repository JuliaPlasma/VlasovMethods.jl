using Test
using TOML
using VlasovMethods

# Every plotting function is declared in `src/plotting.jl` and every method of it is in
# `ext/VlasovMethodsMakieExt.jl`. This file asserts the half of that split a caller without
# CairoMakie has: no plotting name has a method, and each call raises the plain `MethodError` of a
# function with no methods.
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
        const PLOTTING = (VlasovMethods.plot_distribution, VlasovMethods.plot_particles,
            VlasovMethods.plot_series, VlasovMethods.plot_groups, VlasovMethods.plot_damping,
            VlasovMethods.plot_phase_space, VlasovMethods.plot_collision_distribution,
            VlasovMethods.plot_drift, VlasovMethods.plot_entropy, VlasovMethods.plot_error,
            VlasovMethods.plot_convergence, VlasovMethods.plot_cost,
            VlasovMethods.plot_velocity_projections)
        for f in PLOTTING
            @assert isempty(methods(f)) "\$f has a method"
        end

        grid = VlasovMethods.GridDistribution(4, 3, (0.0, 1.0), (-1.0, 1.0))
        particles = VlasovMethods.ParticleDistribution(1, 1, 1)
        for call in (
            () -> VlasovMethods.plot_distribution(grid),
            () -> VlasovMethods.plot_distribution(particles; xdomain = (0.0, 1.0), vdomain = (-1.0, 1.0)),
            () -> VlasovMethods.plot_particles(particles; xdomain = (0.0, 1.0), vdomain = (-1.0, 1.0)),
            () -> VlasovMethods.plot_series([0.0, 1.0], ([0.5, 0.5],)),
            () -> VlasovMethods.plot_groups([0.0, 1.0], ["a" => [[0.5, 0.5]]]),
            () -> VlasovMethods.plot_damping(
                [(t = [0.0, 1.0], y = [0.5, 0.4], peak_t = [0.0], peak_y = [0.5], label = "a")],
                (t = [0.0, 1.0], y = [0.5, 0.3], label = "b")),
            () -> VlasovMethods.plot_phase_space(zeros(2, 1, 1), [1]; xdomain = (0.0, 1.0),
                titles = ["a"]),
            () -> VlasovMethods.plot_collision_distribution(particles, particles;
                xdomain = (0.0, 1.0), vdomain = (-1.0, 1.0)),
            () -> VlasovMethods.plot_drift([0.0, 1.0], ([1.0, 1.0], [1.0, 1.0], [1.0, 1.0])),
            () -> VlasovMethods.plot_entropy([0.0, 1.0], [0.5, 0.4]),
            () -> VlasovMethods.plot_error([0.0, 1.0], [0.1, 0.2], [0.1, 0.2]),
            () -> VlasovMethods.plot_convergence([0.2, 0.1], [0.04, 0.01]),
            () -> VlasovMethods.plot_cost([1.0, 2.0], [0.04, 0.01]),
            () -> VlasovMethods.plot_velocity_projections([0.0, 1.0], [0.5, 0.4],
                [(v = [0.0, 1.0], f = identity, df = identity, title = "a")];
                vdomain = (-1.0, 1.0)))
            exception = try
                call()
                nothing
            catch thrown
                thrown
            end
            @assert exception isa MethodError "an unloaded plotting call threw \$(typeof(exception))"
        end

        print("no Makie in the environment: the extension is not loaded and no plotting name has a method")
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
    @test occursin("package extension", string(@doc VlasovMethods.plot_series))
    @test occursin("load `CairoMakie`", string(@doc VlasovMethods.plot_series))
    for doc in (@doc(VlasovMethods.plot_groups), @doc(VlasovMethods.plot_damping),
        @doc(VlasovMethods.plot_phase_space), @doc(VlasovMethods.plot_collision_distribution),
        @doc(VlasovMethods.plot_drift), @doc(VlasovMethods.plot_entropy),
        @doc(VlasovMethods.plot_error), @doc(VlasovMethods.plot_convergence),
        @doc(VlasovMethods.plot_cost), @doc(VlasovMethods.plot_velocity_projections))
        @test occursin("package extension", string(doc))
        @test occursin("load `CairoMakie`", string(doc))
        @test occursin("-> Figure", string(doc))
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
