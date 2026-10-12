using Test
using TOML
using VlasovMethods

# `plot_distribution` and `plot_particles` are declared in `src/plotting.jl` and every method of
# them is in `ext/VlasovMethodsMakieExt.jl`. This file asserts the half of that split that needs
# no Makie: each name is a function with no method until CairoMakie is loaded, and `Project.toml`
# keeps CairoMakie a weak dependency, so an environment that installs VlasovMethods alone has no
# Makie UUID in its manifest.
#
# The rendering tests are in `test/integration/makie_extension.jl`, which loads CairoMakie. No
# other test file may load it: the first testset here holds only while it is not loaded.

@testset "the plotting names have no method without Makie" begin
    @test Base.get_extension(VlasovMethods, :VlasovMethodsMakieExt) === nothing
    @test isempty(methods(plot_distribution))
    @test isempty(methods(plot_particles))

    grid = GridDistribution(4, 3, (0.0, 1.0), (-1.0, 1.0))
    @test_throws MethodError plot_distribution(grid)

    particles = ParticleDistribution(1, 1, 2)
    @test_throws MethodError plot_distribution(particles; xdomain = (0.0, 1.0),
        vdomain = (-1.0, 1.0))
    @test_throws MethodError plot_particles(particles; xdomain = (0.0, 1.0),
        vdomain = (-1.0, 1.0))
end

@testset "the docstrings name the extension and CairoMakie" begin
    @test occursin("package extension", string(@doc VlasovMethods.plot_distribution))
    @test occursin("load `CairoMakie`", string(@doc VlasovMethods.plot_distribution))
    @test occursin("package extension", string(@doc VlasovMethods.plot_particles))
    @test occursin("load `CairoMakie`", string(@doc VlasovMethods.plot_particles))
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

    # The plotting stacks this work replaces, and the two the plan forbids in every package.
    for name in ("Plots", "StatsPlots", "LaTeXStrings", "Optimisers", "DifferentialEquations")
        @test !haskey(project["deps"], name)
    end
end
