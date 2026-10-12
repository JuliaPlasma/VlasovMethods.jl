using CairoMakie
using ParticleMethods: ParticleList
using Test
using VlasovMethods

# The rendering tests of the plotting extension. Loading CairoMakie activates
# `VlasovMethodsMakieExt`, where every method of `plot_distribution` and `plot_particles` is
# defined; `test/plotting.jl` asserts the other half of that split, that neither name has a
# method until this file's `using` has run.
#
# The axis holds the caller's plots in `ax.scene.plots` and its own furniture -- background,
# spines, ticks -- in `ax.blockscene.plots`, so the tests read the former.

const XDOMAIN = (0.0, 1.0)
const VDOMAIN = (-4.0, 4.0)

# one 1D1V particle of weight `w` at position `x` and velocity `v`
particle(x, v; w = 1.0) = ParticleDistribution(1, 1, ParticleList([x], [v], [w]))

# the one axis of the figure, which holds the one plot the caller asked for
function only_axis(fig)
    @test fig isa Figure
    @test length(fig.content) == 1
    ax = fig.content[1]
    @test ax isa Axis
    @test length(ax.scene.plots) == 1
    return ax
end

# the matrix the heatmap of that axis is drawn from
function plotted_matrix(fig)
    plot = only_axis(fig).scene.plots[1]
    @test plot isa Heatmap
    return plot[3][]
end

@testset "CairoMakie activates the extension" begin
    @test Base.get_extension(VlasovMethods, :VlasovMethodsMakieExt) isa Module
end

@testset "a grid distribution is one untransposed heatmap" begin
    nx, nv = 8, 5
    values = reshape(1.0:(nx * nv), nx, nv)

    fig = plot_distribution(GridDistribution(values, XDOMAIN, VDOMAIN))

    # the columns of the matrix are the `v` nodes, so its transpose fails both of these
    @test size(plotted_matrix(fig)) == (nx, nv)
    @test plotted_matrix(fig) == values

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0
end

@testset "one particle in the middle of the first cell" begin
    nx, nv = 8, 8
    hx = (last(XDOMAIN) - first(XDOMAIN)) / nx
    hv = (last(VDOMAIN) - first(VDOMAIN)) / nv

    dist = particle(first(XDOMAIN) + 0.5 * hx, first(VDOMAIN) + 0.5 * hv)
    F = plotted_matrix(plot_distribution(dist;
        xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv))

    # the normalisation `1/(hx hv)` is what makes the histogram a density
    @test F[1, 1] == 1 / (hx * hv)
    @test count(!iszero, F) == 1
end

@testset "a particle one period away lands in the same cell" begin
    period = last(XDOMAIN) - first(XDOMAIN)
    nx, nv = 8, 8
    hx = period / nx
    hv = (last(VDOMAIN) - first(VDOMAIN)) / nv
    v = first(VDOMAIN) + 0.5 * hv

    wrapped = plotted_matrix(plot_distribution(
        particle(first(XDOMAIN) + 0.5 * hx + period, v);
        xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv))
    centre = plotted_matrix(plot_distribution(
        particle(first(XDOMAIN) + 0.5 * hx, v);
        xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv))

    @test wrapped == centre
    @test wrapped[1, 1] == 1 / (hx * hv)
end

@testset "a position below the domain lands in the last cell" begin
    nx, nv = 8, 8
    hx = (last(XDOMAIN) - first(XDOMAIN)) / nx
    hv = (last(VDOMAIN) - first(VDOMAIN)) / nv

    # `mod(-1e-20, 1.0) == 1.0` in Float64, so the cell index of this particle is one past the
    # last cell without the clamp
    dist = particle(first(XDOMAIN) - 1e-20, first(VDOMAIN) + 0.5 * hv)
    F = plotted_matrix(plot_distribution(dist;
        xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv))

    @test F[nx, 1] == 1 / (hx * hv)
    @test count(!iszero, F) == 1
end

@testset "a velocity outside the domain is dropped, not clamped" begin
    nx, nv = 8, 8
    hx = (last(XDOMAIN) - first(XDOMAIN)) / nx
    hv = (last(VDOMAIN) - first(VDOMAIN)) / nv
    w_inside = 0.25

    # `v = vb` and `v = va - hv` both lie outside the half-open cells; clamping either of them
    # into an end cell, or admitting `v = vb`, changes the sum below
    dist = ParticleDistribution(1, 1,
        ParticleList([0.25, 0.5, 0.75],
            [last(VDOMAIN), first(VDOMAIN) - hv, 0.0],
            [1.0, 1.0, w_inside]))
    F = plotted_matrix(plot_distribution(dist;
        xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv))

    @test sum(F) * hx * hv == w_inside
    @test count(!iszero, F) == 1
end

@testset "plot_particles wraps x into the domain" begin
    period = last(XDOMAIN) - first(XDOMAIN)
    dist = ParticleDistribution(1, 1,
        ParticleList([0.25, 0.25 + period], [-1.0, 1.0], [0.5, 0.5]))

    fig = plot_particles(dist; xdomain = XDOMAIN, vdomain = VDOMAIN)
    ax = only_axis(fig)
    plot = ax.scene.plots[1]
    @test plot isa Scatter

    points = plot[1][]
    @test length(points) == 2
    @test all(p -> first(XDOMAIN) ≤ p[1] < last(XDOMAIN), points)
    @test first.(points) == [0.25, 0.25]
    @test last.(points) == [-1.0, 1.0]
    @test ax.limits[] == (first(XDOMAIN), last(XDOMAIN), first(VDOMAIN), last(VDOMAIN))
end

@testset "a distribution with no particles is refused" begin
    empty = ParticleDistribution(1, 1, 0)
    @test iszero(length(empty.particles))
    @test_throws ArgumentError plot_distribution(empty; xdomain = XDOMAIN, vdomain = VDOMAIN)
    @test_throws ArgumentError plot_particles(empty; xdomain = XDOMAIN, vdomain = VDOMAIN)
end
