using CairoMakie
using ParticleMethods: ParticleList
using ReducedComplexityModeling: ParameterSpace, Parameter
using Test
using VlasovMethods

# The rendering tests of the plotting extension. Loading CairoMakie activates
# `VlasovMethodsMakieExt`, where every method of the plotting functions is defined;
# `test/plotting.jl` asserts the other half of that split, that no plotting name has a method
# until this file's `using` has run.
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

# the one plot the axis of `fig` holds
function only_plot(fig)
    return only_axis(fig).scene.plots[1]
end

# the matrix the heatmap of that axis is drawn from
function plotted_matrix(fig)
    plot = only_plot(fig)
    @test plot isa Heatmap
    return plot[3][]
end

# the phrase of the exception a call raises, so that two guards of the same kind are told apart
function message_of(f)
    try
        f()
        return "no exception"
    catch thrown
        return sprint(showerror, thrown)
    end
end

@testset "CairoMakie activates the extension" begin
    @test Base.get_extension(VlasovMethods, :VlasovMethodsMakieExt) isa Module
end

@testset "a grid distribution is one untransposed heatmap" begin
    nx, nv = 8, 5
    values = reshape(1.0:(nx * nv), nx, nv)
    dist = GridDistribution(values, XDOMAIN, VDOMAIN)

    fig = plot_distribution(dist)
    plot = only_plot(fig)

    # the columns of the matrix are the `v` nodes, so its transpose fails both of these
    @test size(plotted_matrix(fig)) == (nx, nv)
    @test plotted_matrix(fig) == values

    # the heatmap is drawn on the grid's own nodes, not on the array's indices
    @test plot.args[][1] == dist.x
    @test plot.args[][2] == dist.v

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0
end

@testset "the histogram is drawn on the cell edges" begin
    nx, nv = 8, 5
    dist = particle(first(XDOMAIN) + 0.5 * (last(XDOMAIN) - first(XDOMAIN)) / nx, 0.0)

    plot = only_plot(plot_distribution(dist;
        xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv))

    # a matrix of `nx × nv` cell values is drawn on the `nx + 1` and `nv + 1` cell edges, so a
    # heatmap of the same matrix on `nx` nodes is not this one
    @test size(plot[3][]) == (nx, nv)
    @test plot[1][] == collect(range(first(XDOMAIN), last(XDOMAIN); length = nx + 1))
    @test plot[2][] == collect(range(first(VDOMAIN), last(VDOMAIN); length = nv + 1))
end

@testset "the histogram's cell counts default to 100 by 100" begin
    plot = only_plot(plot_distribution(particle(0.5, 0.0);
        xdomain = XDOMAIN, vdomain = VDOMAIN))

    @test size(plot[3][]) == (100, 100)
    @test length(plot[1][]) == 101
    @test length(plot[2][]) == 101
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

@testset "a velocity one rounding step below the upper end lands in the last cell" begin
    nx, nv = 8, 8
    for T in (Float32, Float64)
        xa, xb, va, vb = T(0), T(1), T(-4), T(4)
        hx = (xb - xa) / nx
        hv = (vb - va) / nv
        v = prevfloat(vb)
        @test va ≤ v < vb

        # the half-open cells put this particle in the last cell, but `(v - va)/hv` rounds up to
        # `nv` in both precisions, so the index needs the same clamp as the position's
        dist = ParticleDistribution(1, 1, ParticleList(T[0.5], T[v], T[1]))
        F = plotted_matrix(plot_distribution(dist;
            xdomain = (xa, xb), vdomain = (va, vb), nx = nx, nv = nv))

        @test F[5, nv] == 1 / (hx * hv)
        @test count(!iszero, F) == 1
    end
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

@testset "a non-finite position is dropped from the histogram" begin
    # the wrap maps `NaN` and `Inf` to `NaN`, which no cell index can be taken from
    w_inside = 0.25
    dist = ParticleDistribution(1, 1,
        ParticleList([NaN, Inf, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, w_inside]))
    nx, nv = 8, 8
    hx = (last(XDOMAIN) - first(XDOMAIN)) / nx
    hv = (last(VDOMAIN) - first(VDOMAIN)) / nv
    F = plotted_matrix(plot_distribution(dist;
        xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv))

    @test sum(F) * hx * hv == w_inside
    @test count(!iszero, F) == 1
end

@testset "a non-finite weight is dropped from the histogram" begin
    # the particles share one cell with a finite one, whose weight a `NaN` or `Inf` would hide
    w_inside = 0.25
    dist = ParticleDistribution(1, 1,
        ParticleList([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [NaN, Inf, w_inside]))
    nx, nv = 8, 8
    hx = (last(XDOMAIN) - first(XDOMAIN)) / nx
    hv = (last(VDOMAIN) - first(VDOMAIN)) / nv
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

@testset "a position below the domain is plotted inside the half-open domain" begin
    # `mod(-1e-20, 1.0) == 1.0` in Float64, which is the upper end the half-open domain excludes
    points = only_plot(plot_particles(particle(first(XDOMAIN) - 1e-20, 0.0);
        xdomain = XDOMAIN, vdomain = VDOMAIN))[1][]

    @test length(points) == 1
    @test all(p -> first(XDOMAIN) ≤ p[1] < last(XDOMAIN), points)
end

@testset "a non-finite position is plotted as NaN" begin
    # the wrap cannot place these: `mod` maps them to `NaN`, the `xb` test is false for it, and
    # the axis skips it; neither is moved onto an end
    for bad in (NaN, Inf)
        points = only_plot(plot_particles(particle(bad, 0.0);
            xdomain = XDOMAIN, vdomain = VDOMAIN))[1][]

        @test length(points) == 1
        @test isnan(points[1][1])
    end
end

@testset "a distribution with no particles is refused" begin
    empty = ParticleDistribution(1, 1, 0)
    @test iszero(length(empty.particles))
    @test_throws ArgumentError plot_distribution(empty; xdomain = XDOMAIN, vdomain = VDOMAIN)
    @test_throws ArgumentError plot_particles(empty; xdomain = XDOMAIN, vdomain = VDOMAIN)
end

@testset "plot_distribution refuses a domain or a cell count that has no cells" begin
    dist = particle(0.5, 0.0)
    # a reversed or empty domain would give a negative or an infinite density
    for (xdomain, vdomain) in (((1.0, 0.0), VDOMAIN), (XDOMAIN, (4.0, -4.0)),
        (XDOMAIN, (4.0, 4.0)))
        @test occursin("each domain must run from its lower to its upper end",
            message_of(() -> plot_distribution(dist; xdomain, vdomain)))
    end
    # an infinite end gives a `NaN` cell index, or edges that `range` refuses without naming the
    # domain
    for (xdomain, vdomain) in ((XDOMAIN, (-Inf, 1.0)), ((0.0, Inf), VDOMAIN),
        ((NaN, 1.0), VDOMAIN), (XDOMAIN, (-1.0, NaN)))
        @test occursin("each domain must have finite ends",
            message_of(() -> plot_distribution(dist; xdomain, vdomain)))
    end
    for (nx, nv) in ((0, 4), (4, 0))
        @test occursin("the cell counts must be positive",
            message_of(() -> plot_distribution(dist; xdomain = XDOMAIN, vdomain = VDOMAIN,
                nx, nv)))
    end
end

@testset "the axes are labelled x and v and carry no title" begin
    grid = GridDistribution(reshape(1.0:15.0, 5, 3), XDOMAIN, VDOMAIN)
    middle = particle(0.5, 0.0)

    for fig in (plot_distribution(grid),
        plot_distribution(middle; xdomain = XDOMAIN, vdomain = VDOMAIN),
        plot_particles(middle; xdomain = XDOMAIN, vdomain = VDOMAIN))
        ax = only_axis(fig)

        @test ax.xlabel[] == "x"
        @test ax.ylabel[] == "v"
        @test ax.title[] == ""
    end
end

# `plot_series` holds the `(x, y)` pairs a `lines!` was handed, so the axis carries one plot per
# series and the read-back of plot `i` is the `Point`s of that series against its times.
function series_axis(fig, nplots)
    @test fig isa Figure
    @test length(fig.content) == 1
    ax = fig.content[1]
    @test ax isa Axis
    @test length(ax.scene.plots) == nplots
    return ax
end
function series_points(fig, i)
    plot = series_axis(fig, 2).scene.plots[i]
    @test plot isa Lines
    return plot[1][]
end

@testset "plot_series draws each series against the times, $T" for T in (Float64, Float32)
    t = collect(T, 0:0.25:2)
    s1, s2 = sin.(t), cos.(t)

    # the axis labels default to `t` and nothing, and the keywords set them
    fig = plot_series(t, (s1, s2))
    ax = series_axis(fig, 2)
    @test ax.xlabel[] == "t"
    @test ax.ylabel[] == ""
    @test ax.title[] == ""
    p1, p2 = series_points(fig, 1), series_points(fig, 2)
    @test first.(p1) == t && last.(p1) == s1
    @test first.(p2) == t && last.(p2) == s2

    labelled = plot_series(t, (s1, s2); xlabel = "parameter", ylabel = "γ")
    @test series_axis(labelled, 2).xlabel[] == "parameter"
    @test series_axis(labelled, 2).ylabel[] == "γ"

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0

    # a single series is one plot, and a vector of series is as good as a tuple
    @test length(series_axis(plot_series(t, (s1,)), 1).scene.plots) == 1
    @test length(series_axis(plot_series(t, [s1, s2]), 2).scene.plots) == 2

    # the labels name the lines in a legend, in the order of the series
    named = plot_series(t, (s1, s2); labels = ["sin", "cos"])
    @test length(named.content) == 2
    ax, legend = named.content
    @test ax isa Axis
    @test legend isa Legend
    @test [p.label[] for p in ax.scene.plots] == ["sin", "cos"]
    @test [e.label[] for e in only(legend.entrygroups[])[2]] == ["sin", "cos"]
    save(file, named)
    @test filesize(file) > 0

    # no labels, no legend
    @test length(fig.content) == 1

    # a label count that differs from the series count would name the wrong lines
    @test_throws DimensionMismatch plot_series(t, (s1, s2); labels = ["sin"])
    @test_throws DimensionMismatch plot_series(t, (s1,); labels = ["sin", "cos"])

    # an empty collection holds no line, and a time vector of the wrong length would pair the
    # values with the wrong abscissae
    @test_throws ArgumentError plot_series(t, ())
    @test_throws DimensionMismatch plot_series(t, (s1, s2[1:(end - 1)]))
    @test_throws DimensionMismatch plot_series(t, (vcat(s1, s1),))

    # an empty line holds no value and is refused before its length is compared with `t`, so an
    # empty line under an empty `t` is refused too, while a non-empty line under an empty `t` is
    # the length mismatch. The phrases tell the two `ArgumentError`s apart, which the type alone
    # does not.
    @test occursin("no series to plot", message_of(() -> plot_series(t, ())))
    @test occursin("a line holds no value", message_of(() -> plot_series(T[], (T[],))))
    @test occursin("a line holds no value", message_of(() -> plot_series(t, (T[],))))
    @test_throws DimensionMismatch plot_series(T[], ([one(T)],))
end

@testset "plot_particles draws a saved state of the snapshots, $T" for T in (Float64, Float32)
    n, nₛ, nh = 5, 3, 2
    space = ParameterSpace(Parameter(:β, T(6), T(7), 2))
    SS = Snapshots(T, 1, n, nh, nₛ, space)

    # a state that lies inside the domains, so the plotted data are the input itself
    SS.X[1, :, 2, 1] .= T[0.1, 0.2, 0.3, 0.4, 0.5]
    SS.V[1, :, 2, 1] .= T[-1, -2, -3, -4, -5]

    fig = plot_particles(SS, 2, 1; xdomain = XDOMAIN, vdomain = VDOMAIN)
    ax = only_axis(fig)
    plot = ax.scene.plots[1]
    @test plot isa Scatter

    points = plot[1][]
    @test length(points) == n
    @test first.(points) == SS.X[1, :, 2, 1]
    @test last.(points) == SS.V[1, :, 2, 1]
    @test ax.limits[] == (first(XDOMAIN), last(XDOMAIN), first(VDOMAIN), last(VDOMAIN))

    # the plotted state is that of the parameter `p` the caller asked for, not of the first one
    SS.X[1, :, 2, 2] .= T[0.5, 0.4, 0.3, 0.2, 0.1]
    SS.V[1, :, 2, 2] .= T[1, 2, 3, 4, 5]
    other = only_plot(plot_particles(SS, 2, 2; xdomain = XDOMAIN, vdomain = VDOMAIN))
    @test first.(other[1][]) == SS.X[1, :, 2, 2]
    @test last.(other[1][]) == SS.V[1, :, 2, 2]
    @test first.(other[1][]) != first.(points)

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0

    # `ts` and `p` index the saved times and the parameters of the snapshots
    @test_throws ArgumentError plot_particles(
        SS, 0, 1; xdomain = XDOMAIN, vdomain = VDOMAIN)
    @test_throws ArgumentError plot_particles(
        SS, nₛ + 1, 1; xdomain = XDOMAIN, vdomain = VDOMAIN)
    @test_throws ArgumentError plot_particles(
        SS, 1, 0; xdomain = XDOMAIN, vdomain = VDOMAIN)
    @test_throws ArgumentError plot_particles(
        SS, 1, 3; xdomain = XDOMAIN, vdomain = VDOMAIN)
    @test occursin("1:3", message_of(() -> plot_particles(
        SS, 4, 1; xdomain = XDOMAIN, vdomain = VDOMAIN)))
    @test occursin("1:2", message_of(() -> plot_particles(
        SS, 1, 3; xdomain = XDOMAIN, vdomain = VDOMAIN)))

    # `x` wraps into `xdomain`, as the particle method's does: one period away is the same point
    period = last(XDOMAIN) - first(XDOMAIN)
    SS.X[1, :, 3, 1] .= T[0.25 + period, 0.5, 0.75, 0.25, 0.5]
    wrapped = only_plot(plot_particles(SS, 3, 1; xdomain = XDOMAIN, vdomain = VDOMAIN))[1][]
    @test first.(wrapped) == T[0.25, 0.5, 0.75, 0.25, 0.5]
end

# The two panels of a collision figure, each an axis that holds a heatmap and a contour.
function collision_panels(fig)
    @test fig isa Figure
    @test length(fig.content) == 2
    return map(fig.content) do ax
        @test ax isa Axis
        @test length(ax.scene.plots) == 2
        heat, contour = ax.scene.plots
        @test heat isa Heatmap
        @test contour isa Contour
        (ax, heat, contour)
    end
end

@testset "plot_collision_distribution draws both histograms and the Maxwellian of f" begin
    nx, nv = 8, 16
    hx = (last(XDOMAIN) - first(XDOMAIN)) / nx
    hv = (last(VDOMAIN) - first(VDOMAIN)) / nv

    # weights that sum to 3.5, at `v = ±1` with the mean zero: the Maxwellian of `f` has
    # `n = 3.5` on the unit domain, `u = 0` and `T = 1`
    f = ParticleDistribution(1, 1,
        ParticleList([0.1, 0.4, 0.6, 0.9], [-1.0, 1.0, -1.0, 1.0], [0.75, 0.75, 1.0, 1.0]))
    g = ParticleDistribution(1, 1,
        ParticleList([0.2, 0.5, 0.7], [0.0, 0.5, -0.5], [1.0, 1.0, 1.0]))
    keywords = (; xdomain = XDOMAIN, vdomain = VDOMAIN, nx = nx, nv = nv)

    fig = plot_collision_distribution(f, g; keywords...)
    (initial, hf, cf), (final, hg, cg) = collision_panels(fig)
    @test (initial.title[], final.title[]) == ("initial", "final")

    # each histogram is the one `plot_distribution` draws of the same input, and integrates to
    # the particles' total weight
    @test hf[3][] == plotted_matrix(plot_distribution(f; keywords...))
    @test hg[3][] == plotted_matrix(plot_distribution(g; keywords...))
    @test sum(hf[3][]) * hx * hv ≈ 3.5
    @test sum(hg[3][]) * hx * hv ≈ 3.0

    # the Maxwellian of `f` on the cell centres, on both panels
    xc = first(XDOMAIN) .+ hx .* ((1:nx) .- 0.5)
    vc = first(VDOMAIN) .+ hv .* ((1:nv) .- 0.5)
    M = [3.5 / sqrt(2π) * exp(-v^2 / 2) for x in xc, v in vc]
    @test cf[3][] ≈ M
    @test cg[3][] ≈ M

    # weights that differ and a mean that is not zero, so that the weighted moments differ from
    # the unweighted ones and from a temperature about zero: `n = W = 4` on the unit domain,
    # `u = (3·0 + 1·2)/4 = 1/2` and `T = (3·(1/2)² + 1·(3/2)²)/4 = 3/4`, where the unweighted
    # mean, the unweighted temperature and the temperature about zero are each `1`
    skewed = ParticleDistribution(1, 1, ParticleList([0.25, 0.75], [0.0, 2.0], [3.0, 1.0]))
    (_, _, cs), _ = collision_panels(plot_collision_distribution(skewed, g; keywords...))
    @test cs[3][] ≈
          [4 / sqrt(2π * 0.75) * exp(-(v - 0.5)^2 / (2 * 0.75)) for x in xc, v in vc]

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0

    # no particle leaves nothing to draw, and the error names the empty distribution; one
    # particle has no temperature and no Maxwellian
    empty = ParticleDistribution(1, 1, 0)
    @test_throws ArgumentError plot_collision_distribution(empty, g; keywords...)
    @test occursin("`f` holds no particle",
        message_of(() -> plot_collision_distribution(empty, g; keywords...)))
    @test_throws ArgumentError plot_collision_distribution(f, empty; keywords...)
    @test occursin("`g` holds no particle",
        message_of(() -> plot_collision_distribution(f, empty; keywords...)))
    @test_throws ArgumentError plot_collision_distribution(particle(0.5, 0.0), g; keywords...)
    @test occursin("no positive temperature",
        message_of(() -> plot_collision_distribution(particle(0.5, 0.0), g; keywords...)))

    # a particle with a non-finite velocity takes no part in the moments, so a distribution of
    # such particles alone has no mass
    unplaced = ParticleDistribution(1, 1, ParticleList([0.5, 0.25], [NaN, Inf], [1.0, 1.0]))
    @test_throws ArgumentError plot_collision_distribution(unplaced, g; keywords...)
    @test occursin("no positive mass",
        message_of(() -> plot_collision_distribution(unplaced, g; keywords...)))

    # a particle whose weight is not finite is dropped from the histogram and from the moments,
    # so both panels are those of the distributions without it
    for w in (NaN, Inf)
        f_w = ParticleDistribution(1, 1,
            ParticleList([0.1, 0.4, 0.6, 0.9, 0.3], [-1.0, 1.0, -1.0, 1.0, 0.5],
                [0.75, 0.75, 1.0, 1.0, w]))
        g_w = ParticleDistribution(1, 1,
            ParticleList([0.2, 0.5, 0.7, 0.3], [0.0, 0.5, -0.5, 0.5], [1.0, 1.0, 1.0, w]))
        (_, hf_w, cf_w), (_, hg_w, cg_w) = collision_panels(
            plot_collision_distribution(f_w, g_w; keywords...))
        @test hf_w[3][] == hf[3][]
        @test hg_w[3][] == hg[3][]
        @test cf_w[3][] ≈ M
        @test cg_w[3][] ≈ M
    end

    # a reversed or empty domain would give a negative or an infinite density, and no cell
    # holds a particle where a cell count is zero
    for (xdomain, vdomain) in (((1.0, 0.0), VDOMAIN), (XDOMAIN, (4.0, -4.0)),
        (XDOMAIN, (4.0, 4.0)))
        @test_throws ArgumentError plot_collision_distribution(f, g; xdomain, vdomain)
        @test occursin("each domain must run from its lower to its upper end",
            message_of(() -> plot_collision_distribution(f, g; xdomain, vdomain)))
    end
    for (nx, nv) in ((0, 4), (4, 0))
        @test_throws ArgumentError plot_collision_distribution(f, g;
            xdomain = XDOMAIN, vdomain = VDOMAIN, nx, nv)
        @test occursin("the cell counts must be positive",
            message_of(() -> plot_collision_distribution(f, g;
                xdomain = XDOMAIN, vdomain = VDOMAIN, nx, nv)))
    end
    # an infinite end gives a `NaN` cell index, or edges that `range` refuses without naming the
    # domain; a `NaN` end is no end, and its orientation is not defined
    for (xdomain, vdomain) in ((XDOMAIN, (-Inf, Inf)), ((0.0, Inf), VDOMAIN),
        ((NaN, 1.0), VDOMAIN))
        @test_throws ArgumentError plot_collision_distribution(f, g; xdomain, vdomain)
        @test occursin("each domain must have finite ends",
            message_of(() -> plot_collision_distribution(f, g; xdomain, vdomain)))
    end
end

@testset "plot_drift draws the relative drift on a log axis, $T" for T in (Float64, Float32)
    t = collect(T, 0:0.5:2)
    mass = T[2, 2, 2, 2, 2]
    momentum = T[1, 1.5, 0.5, 1, 3]
    energy = T[4, 5, 6, 2, 4]

    fig = plot_drift(t, (mass, momentum, energy))
    @test fig isa Figure
    @test length(fig.content) == 2
    ax, legend = fig.content
    @test ax isa Axis
    @test legend isa Legend
    @test ax.yscale[] === log10
    @test length(ax.scene.plots) == 3
    @test [p.label[] for p in ax.scene.plots] == ["mass", "momentum", "energy"]

    # `|Q - Q(t₁)| / |Q(t₁)|`: a zero is dropped from a line with a positive drift, and the line
    # of the conserved mass, zero at every time, is a decade below the smallest drift, 0.25
    floor = T(0.25) / 10
    expected = ((t, fill(floor, 5)), (t[[2, 3, 5]], T[0.5, 0.5, 2]),
        (t[[2, 3, 4]], T[0.25, 0.5, 0.5]))
    for (plot, (times, drift)) in zip(ax.scene.plots, expected)
        @test plot isa Lines
        points = plot[1][]
        @test first.(points) == times
        @test last.(points) == drift
    end

    # the axis spans the drawn values, and not the hundreds of decades down to `floatmin`
    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0
    @test ax.finallimits[].origin[2] > floor / 10

    # a non-finite value or time drops its pair, and the rest of the line is drawn
    dropped = plot_drift(t, (T[2, 3, NaN, 4, 3],); labels = ["mass"])
    @test first.(only(dropped.content[1].scene.plots)[1][]) == t[[2, 4, 5]]
    t_inf = copy(t)
    t_inf[2] = Inf
    dropped = plot_drift(t_inf, (mass,); labels = ["mass"])
    @test first.(only(dropped.content[1].scene.plots)[1][]) == t[[1, 3, 4, 5]]

    # a series of another length, and a label count that differs, name the shorter one
    @test_throws DimensionMismatch plot_drift(t, (mass, momentum[1:(end - 1)], energy))
    @test occursin("`momentum` is shorter than `t`",
        message_of(() -> plot_drift(t, (mass, momentum[1:(end - 1)], energy))))
    @test_throws DimensionMismatch plot_drift(t[1:(end - 1)], (mass, momentum, energy))
    @test occursin("`t` is shorter than `mass`",
        message_of(() -> plot_drift(t[1:(end - 1)], (mass, momentum, energy))))
    two_labels = ["mass", "momentum"]
    @test_throws DimensionMismatch plot_drift(t, (mass, momentum, energy); labels = two_labels)
    @test occursin("`labels` is shorter than `quantities`",
        message_of(() -> plot_drift(t, (mass, momentum, energy); labels = two_labels)))

    # nothing to draw, an initial value with no relative drift, and the single drift zero
    @test_throws ArgumentError plot_drift(t, ())
    @test occursin("`quantities` holds no quantity", message_of(() -> plot_drift(t, ())))
    @test_throws ArgumentError plot_drift(T[], (T[], T[], T[]))
    @test occursin("`t` holds no value", message_of(() -> plot_drift(T[], (T[], T[], T[]))))
    @test_throws ArgumentError plot_drift(t, (T[0, 1, 1, 1, 1],); labels = ["momentum"])
    @test occursin("relative drift is not defined",
        message_of(() -> plot_drift(t, (T[0, 1, 1, 1, 1],); labels = ["momentum"])))
    @test_throws ArgumentError plot_drift(t, (T[NaN, 1, 1, 1, 1],); labels = ["mass"])
    @test occursin("relative drift is not defined",
        message_of(() -> plot_drift(t, (T[NaN, 1, 1, 1, 1],); labels = ["mass"])))
    @test_throws ArgumentError plot_drift(T[0], (T[2],); labels = ["mass"])
    @test occursin("single value zero",
        message_of(() -> plot_drift(T[0], (T[2],); labels = ["mass"])))
end

@testset "plot_entropy draws S_h on a linear axis, $T" for T in (Float64, Float32)
    t = collect(T, 0:0.5:2)
    S = T[1, 0.8, 0.7, 0.65, 0.62]

    fig = plot_entropy(t, S)
    ax = only_axis(fig)
    @test ax.yscale[] === identity
    plot = ax.scene.plots[1]
    @test plot isa Lines
    @test first.(plot[1][]) == t
    @test last.(plot[1][]) == S

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0

    # a non-finite value drops its pair, and a series of one value is a marker
    dropped = only_plot(plot_entropy(t, T[1, Inf, 0.7, NaN, 0.62]))
    @test first.(dropped[1][]) == t[[1, 3, 5]]
    single = only_plot(plot_entropy(T[0], T[1]))
    @test single isa Scatter
    @test single[1][] == [Point2(zero(T), one(T))]

    @test_throws DimensionMismatch plot_entropy(t, S[1:(end - 1)])
    @test occursin("`S` is shorter than `t`", message_of(() -> plot_entropy(t, S[1:(end - 1)])))
    @test_throws ArgumentError plot_entropy(T[], T[])
    @test occursin("`t` holds no value", message_of(() -> plot_entropy(T[], T[])))
    @test_throws ArgumentError plot_entropy(t, fill(T(NaN), 5))
    @test occursin("no pair of `t` and `S` is finite",
        message_of(() -> plot_entropy(t, fill(T(NaN), 5))))
end

@testset "plot_error draws the error and the closed form, $T" for T in (Float64, Float32)
    t = collect(T, 0:0.5:2)
    closed = exp.(-t)
    err = T[0, 1e-3, 2e-3, 1e-3, 5e-4]

    fig = plot_error(t, err, closed)
    @test fig isa Figure
    @test length(fig.content) == 2
    ax = fig.content[1]
    @test ax isa Axis
    @test ax.yscale[] === log10
    @test length(ax.scene.plots) == 2
    @test all(p -> p isa Lines, ax.scene.plots)
    @test [p.label[] for p in ax.scene.plots] == ["error", "closed form"]
    error_line, closed_line = ax.scene.plots
    # the zero error at `t₁` is dropped from the error's line
    @test first.(error_line[1][]) == t[2:end]
    @test last.(error_line[1][]) == err[2:end]
    @test last.(closed_line[1][]) == closed

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0

    # a series of one value is a marker, unless that value is zero
    single = plot_error(T[0], T[1e-3], T[1])
    @test all(p -> p isa Scatter, single.content[1].scene.plots)
    @test_throws ArgumentError plot_error(T[0], T[0], T[1])
    @test occursin("`err` holds the single value zero",
        message_of(() -> plot_error(T[0], T[0], T[1])))

    # a negative value has no logarithm
    @test_throws ArgumentError plot_error(t, -err, closed)
    @test occursin("`err` holds a negative value", message_of(() -> plot_error(t, -err, closed)))

    # a non-finite value drops its pair
    dropped = plot_error(t, T[NaN, 1e-3, 2e-3, 1e-3, 5e-4], closed)
    @test first.(dropped.content[1].scene.plots[1][1][]) == t[2:end]

    @test_throws DimensionMismatch plot_error(t, err, closed[1:(end - 1)])
    @test occursin("`closed` is shorter than `t`",
        message_of(() -> plot_error(t, err, closed[1:(end - 1)])))
    @test_throws ArgumentError plot_error(T[], T[], T[])
    @test occursin("`t` holds no value", message_of(() -> plot_error(T[], T[], T[])))
end

@testset "plot_convergence fits C h^p by least squares, $T" for T in (Float64, Float32)
    h = T[0.4, 0.2, 0.1, 0.05]
    err = T[0.17, 0.035, 0.011, 0.0026]

    fig = plot_convergence(h, err)
    @test fig isa Figure
    @test length(fig.content) == 1
    ax = fig.content[1]
    @test ax isa Axis
    @test ax.xscale[] === log10
    @test ax.yscale[] === log10
    @test length(ax.scene.plots) == 2
    markers, line = ax.scene.plots
    @test markers isa Scatter
    @test line isa Lines
    @test first.(markers[1][]) == h
    @test last.(markers[1][]) == err

    # the least-squares line of `ln err` against `ln h`, here from the normal equations of the
    # design matrix rather than from the centred sums of the code
    A = [ones(length(h)) log.(Float64.(h))]
    lnC, p = A \ log.(Float64.(err))
    # the line through the first and the last pair has another slope on these data, so a label
    # that held it fails the comparison below
    p_ends = log(err[end] / err[1]) / log(h[end] / h[1])
    @test abs(p_ends - p) > 0.02
    label = match(r"p = (\S+)$", ax.ylabel[])
    @test !isnothing(label)
    @test abs(parse(Float64, label.captures[1]) - p) ≤ 0.006
    @test first.(line[1][]) == sort(h)
    @test last.(line[1][]) ≈ exp(lnC) .* sort(Float64.(h)) .^ p rtol = 1e-4

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0

    # a zero error is neither drawn nor fitted
    zeroed = plot_convergence(h, T[0.17, 0.035, 0.011, 0])
    zax = zeroed.content[1]
    @test first.(zax.scene.plots[1][1][]) == h[1:3]
    @test last.(zax.scene.plots[1][1][]) == T[0.17, 0.035, 0.011]
    p3 = last(A[1:3, :] \ log.(Float64.(err[1:3])))
    @test abs(parse(Float64, match(r"p = (\S+)$", zax.ylabel[]).captures[1]) - p3) ≤ 0.006

    # fewer than two finite pairs, or one width only, leave no slope
    @test_throws ArgumentError plot_convergence(T[0.1], T[0.01])
    @test occursin("fewer than two pairs", message_of(() -> plot_convergence(T[0.1], T[0.01])))
    @test occursin("fewer than two pairs",
        message_of(() -> plot_convergence(T[0.1, NaN], T[0.01, 0.02])))
    @test_throws ArgumentError plot_convergence(T[0.1, 0.1], T[0.01, 0.02])
    @test occursin("every `h` is the same",
        message_of(() -> plot_convergence(T[0.1, 0.1], T[0.01, 0.02])))
    # the rounded mean of three `log(0.06)` in `Float64` misses each of them by an ulp, so a test
    # of the deviations from the mean fits a slope of about -9e15
    @test occursin("every `h` is the same",
        message_of(() -> plot_convergence(T[0.06, 0.06, 0.06], T[0.01, 0.02, 0.03])))
    # a negative error is named as such, also where it leaves fewer than two pairs to fit
    @test_throws ArgumentError plot_convergence(T[0.2, 0.1], T[-1, 0.01])
    @test occursin("`err` holds a negative value",
        message_of(() -> plot_convergence(T[0.2, 0.1], T[-1, 0.01])))

    @test_throws DimensionMismatch plot_convergence(h, err[1:(end - 1)])
    @test occursin("`err` is shorter than `h`",
        message_of(() -> plot_convergence(h, err[1:(end - 1)])))
    @test_throws ArgumentError plot_convergence(T[], T[])
    @test occursin("`h` holds no value", message_of(() -> plot_convergence(T[], T[])))
end

@testset "plot_cost draws one marker per run on logarithmic axes, $T" for T in (Float64, Float32)
    cost = T[0.5, 2, 8, 30]
    err = T[1e-2, 3e-3, 1e-3, 2e-4]

    fig = plot_cost(cost, err)
    ax = only_axis(fig)
    @test ax.xscale[] === log10
    @test ax.yscale[] === log10
    plot = ax.scene.plots[1]
    @test plot isa Scatter
    @test first.(plot[1][]) == cost
    @test last.(plot[1][]) == err

    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0

    # a run with a non-finite cost is dropped, and a single run is one marker
    @test first.(only_plot(plot_cost(T[0.5, Inf, 8, 30], err))[1][]) == cost[[1, 3, 4]]
    @test length(only_plot(plot_cost(T[2], T[1e-3]))[1][]) == 1

    @test_throws DimensionMismatch plot_cost(cost, err[1:(end - 1)])
    @test occursin("`err` is shorter than `cost`",
        message_of(() -> plot_cost(cost, err[1:(end - 1)])))
    @test_throws ArgumentError plot_cost(T[], T[])
    @test occursin("`cost` holds no value", message_of(() -> plot_cost(T[], T[])))
    @test_throws ArgumentError plot_cost(-cost, err)
    @test occursin("`cost` holds a negative value", message_of(() -> plot_cost(-cost, err)))

    # a run with a zero error is dropped where another error is positive, and no run is left
    # where the zero cost and the zero error fall on different runs
    @test first.(only_plot(plot_cost(cost, T[1e-2, 0, 1e-3, 2e-4]))[1][]) == cost[[1, 3, 4]]
    @test_throws ArgumentError plot_cost(T[1, 0], T[0, 1])
    @test occursin("no run has a positive `cost` and a positive `err`",
        message_of(() -> plot_cost(T[1, 0], T[0, 1])))
end

# The limits that are set on the axis of `fig`, `(x, y)` with `nothing` where Makie chooses them,
# after `fig` is saved.
function saved_log_limits(fig)
    ax = fig.content[1]
    file = tempname() * ".png"
    save(file, fig)
    @test filesize(file) > 0
    return ax.limits[]
end

# The limits of an axis whose every value is `c`: a range of positive width around `c`, within
# the positive finite `Float64`s, so that neither end is zero or infinite.
function test_flat_limits(limits, c)
    @test limits isa Tuple{Real, Real}
    lo, hi = limits
    @test floatmin(Float64) ≤ lo ≤ c ≤ hi ≤ floatmax(Float64)
    @test lo < hi
end

@testset "a log axis whose every value is one value is saved, $T" for T in (Float64, Float32)
    t = collect(T, 0:0.5:1)
    tiny = Float64(floatmin(T))

    # exactly conserved quantities: every drift is zero, drawn at `floatmin`
    xlims, ylims = saved_log_limits(plot_drift(t, (T[2, 2, 2],); labels = ["mass"]))
    @test isnothing(xlims)
    test_flat_limits(ylims, tiny)
    _, ylims = saved_log_limits(plot_drift(t, (T[2, 2, 2], T[1, 1, 1], T[4, 4, 4])))
    test_flat_limits(ylims, tiny)

    # an error and a closed form that are zero, or one value of each that is the same
    _, ylims = saved_log_limits(plot_error(t, zeros(T, 3), zeros(T, 3)))
    test_flat_limits(ylims, tiny)
    _, ylims = saved_log_limits(plot_error(T[0], T[1e-30], T[1e-30]))
    test_flat_limits(ylims, Float64(T(1e-30)))

    # the costs and the errors of runs, each axis on its own
    xlims, ylims = saved_log_limits(plot_cost(T[1, 2, 4], zeros(T, 3)))
    @test isnothing(xlims)
    test_flat_limits(ylims, tiny)
    xlims, ylims = saved_log_limits(plot_cost(zeros(T, 3), T[1e-2, 1e-3, 1e-4]))
    test_flat_limits(xlims, tiny)
    @test isnothing(ylims)
    _, ylims = saved_log_limits(plot_cost(T[1, 2, 4], fill(floatmax(T), 3)))
    test_flat_limits(ylims, Float64(floatmax(T)))

    # an error that is the same at each width, which the fitted line of slope zero meets; at
    # `1e-170` Makie's own limits reach zero
    xlims, ylims = saved_log_limits(plot_convergence(T[0.5, 2], T[1, 1]))
    @test isnothing(xlims)
    test_flat_limits(ylims, 1.0)
    if T == Float64
        _, ylims = saved_log_limits(plot_convergence([0.5, 2.0], [1e-170, 1e-170]))
        test_flat_limits(ylims, 1e-170)
    end

    # where the values differ, Makie's own limits stay: a conserved mass beside a drifting
    # energy, an error beside a closed form of another value, and errors that fall with `h`
    drifting = (T[2, 2, 2], T[4, 5, 4])
    _, ylims = saved_log_limits(plot_drift(t, drifting; labels = ["m", "e"]))
    @test isnothing(ylims)
    _, ylims = saved_log_limits(plot_error(t, zeros(T, 3), T[1, 0.5, 0.25]))
    @test isnothing(ylims)
    _, ylims = saved_log_limits(plot_convergence(T[0.4, 0.2, 0.1], T[0.1, 0.02, 0.005]))
    @test isnothing(ylims)
end
