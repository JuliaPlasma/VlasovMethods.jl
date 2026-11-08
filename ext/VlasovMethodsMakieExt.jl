module VlasovMethodsMakieExt

using CairoMakie
using VlasovMethods

# `x` is periodic on `xdomain`, which is any object with a `first` and a `last`.
@inline function _wrap(x, xdomain)
    xa = first(xdomain)
    return xa + mod(x - xa, last(xdomain) - xa)
end

function VlasovMethods.plot_distribution(dist::GridDistribution)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "v")
    heatmap!(ax, dist.x, dist.v, dist.values)
    return fig
end

# The weighted particle histogram of `plot_distribution`: the cell edges in `x` and `v`, and the
# density of each of the `nx × nv` cells.
function _histogram(dist::ParticleDistribution{DT, 1, 1}, xdomain, vdomain,
        nx::Integer, nv::Integer) where {DT}
    # a domain end that is not finite leaves a cell index `NaN` or makes `range` refuse the
    # edges, and a reversed domain gives a negative density
    nx ≥ 1 && nv ≥ 1 || throw(ArgumentError("the cell counts must be positive"))
    all(isfinite, (first(xdomain), last(xdomain), first(vdomain), last(vdomain))) ||
        throw(ArgumentError("each domain must have finite ends"))
    first(xdomain) < last(xdomain) && first(vdomain) < last(vdomain) ||
        throw(ArgumentError("each domain must run from its lower to its upper end"))
    # `isempty` iterates a `ParticleList` and raises `BoundsError` on an empty list, so the
    # particle methods check its length directly. ParticleMethods K2 records the defect; the
    # fix is on its `main` and is not released.
    iszero(length(dist.particles)) && throw(ArgumentError("no particles to plot"))

    xa, xb = first(xdomain), last(xdomain)
    va, vb = first(vdomain), last(vdomain)
    hx = (xb - xa) / nx
    hv = (vb - va) / nv

    F = zeros(DT, nx, nv)
    for k in eachindex(dist.particles)
        x = dist.particles.x[1, k]
        v = dist.particles.v[1, k]
        w = dist.particles.w[1, k]
        # the wrap cannot place a non-finite position, no cell holds a `v` outside the domain,
        # and a non-finite weight would be the value of its whole cell
        isfinite(x) && va ≤ v < vb && isfinite(w) || continue
        # `mod(-1e-20, 2π) == 2π` in `Float64`, so the last cell needs the clamp.
        i = min(floor(Int, (_wrap(x, xdomain) - xa) / hx) + 1, nx)
        # `(v - va)/hv` can round up to `nv` for `v = prevfloat(vb)`, so the last cell needs the
        # clamp as well.
        j = min(floor(Int, (v - va) / hv) + 1, nv)
        F[i, j] += w / (hx * hv)
    end

    return range(xa, xb; length = nx + 1), range(va, vb; length = nv + 1), F
end

function VlasovMethods.plot_distribution(dist::ParticleDistribution{DT, 1, 1};
        xdomain, vdomain, nx::Integer = 100, nv::Integer = 100) where {DT}
    xedges, vedges, F = _histogram(dist, xdomain, vdomain, nx, nv)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "v")
    heatmap!(ax, xedges, vedges, F)
    return fig
end

function VlasovMethods.plot_particles(dist::ParticleDistribution{DT, 1, 1};
        xdomain, vdomain) where {DT}
    iszero(length(dist.particles)) && throw(ArgumentError("no particles to plot"))

    xa, xb = first(xdomain), last(xdomain)

    # `mod(-1e-20, 1.0) == 1.0` in `Float64` leaves a position at `xb`, which the half-open
    # domain excludes; the two ends of a periodic domain are one point, so it is plotted at `xa`.
    # A position that the wrap cannot place, a non-finite one, leaves `mod` as `NaN`, and neither
    # the `xb` test nor the axis places it: it is plotted as `NaN`.
    x = map(vec(dist.particles.x)) do x₀
        x̃ = _wrap(x₀, xdomain)
        x̃ == xb ? xa : x̃
    end

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = "x", ylabel = "v",
        limits = (xa, xb, first(vdomain), last(vdomain)))
    scatter!(ax, x, vec(dist.particles.v))
    return fig
end

function VlasovMethods.plot_particles(SS::Snapshots, ts::Integer, p::Integer;
        xdomain, vdomain)
    1 ≤ ts ≤ size(SS.X, 3) || throw(ArgumentError(
        "the saved times are 1:$(size(SS.X, 3)); got ts = $ts"))
    1 ≤ p ≤ size(SS.X, 4) || throw(ArgumentError(
        "the parameters are 1:$(size(SS.X, 4)); got p = $p"))

    xa, xb = first(xdomain), last(xdomain)

    # the wrap of the particle method above: the two ends of a periodic domain are one point, so a
    # position that wraps exactly to `xb` is plotted at `xa`
    x = map(SS.X[1, :, ts, p]) do x₀
        x̃ = _wrap(x₀, xdomain)
        x̃ == xb ? xa : x̃
    end

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = "x", ylabel = "v",
        limits = (xa, xb, first(vdomain), last(vdomain)))
    scatter!(ax, x, SS.V[1, :, ts, p])
    return fig
end

# The series of `plot_series` and of each group of `plot_groups`, and the labels of `plot_series`,
# checked before anything is drawn: an empty collection has no line to draw, an empty line holds
# no value, a line of another length would pair its values with the wrong abscissae, and a label
# count that differs from the series count would name the wrong lines.
function _check_series(t, series, labels)
    isempty(series) && throw(ArgumentError("no series to plot"))
    for y in series
        isempty(y) && throw(ArgumentError("a line holds no value"))
        length(y) == length(t) || throw(DimensionMismatch(
            "a series holds $(length(y)) values; the abscissae hold $(length(t))"))
    end
    isnothing(labels) || length(labels) == length(series) ||
        throw(DimensionMismatch(
            "the labels hold $(length(labels)) names; the series hold $(length(series)) lines"))
    return nothing
end

function VlasovMethods.plot_series(t::AbstractVector, series;
        labels = nothing, xlabel = "t", ylabel = "", yscale = identity)
    _check_series(t, series, labels)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, yscale = yscale)
    for (i, y) in enumerate(series)
        if isnothing(labels)
            lines!(ax, t, y)
        else
            lines!(ax, t, y; label = labels[i])
        end
    end
    isnothing(labels) || axislegend(ax)
    return fig
end

# The density, the mean velocity and the temperature of the Maxwellian that carries the mass, the
# momentum and the energy of `dist`, over the particles whose position, velocity and weight are
# finite.
function _maxwellian_moments(dist::ParticleDistribution{DT, 1, 1}, xdomain) where {DT}
    finite = isfinite.(dist.particles.x[1, :]) .& isfinite.(dist.particles.v[1, :]) .&
             isfinite.(dist.particles.w[1, :])
    w = dist.particles.w[1, finite]
    v = dist.particles.v[1, finite]
    W = sum(w; init = zero(DT))
    W > 0 || throw(ArgumentError("the initial distribution has no positive mass"))
    u = sum(w .* v) / W
    T = sum(w .* (v .- u) .^ 2) / W
    T > 0 || throw(ArgumentError("the initial distribution has no positive temperature"))
    return W / (last(xdomain) - first(xdomain)), u, T
end

function VlasovMethods.plot_collision_distribution(f::ParticleDistribution{DT, 1, 1},
        g::ParticleDistribution{DT, 1, 1};
        xdomain, vdomain, nx::Integer = 100, nv::Integer = 100) where {DT}
    # `_histogram` refuses an empty distribution without naming it, and `isempty` raises
    # `BoundsError` on an empty `ParticleList` (ParticleMethods K2)
    iszero(length(f.particles)) && throw(ArgumentError("`f` holds no particle"))
    iszero(length(g.particles)) && throw(ArgumentError("`g` holds no particle"))
    histograms = map(dist -> _histogram(dist, xdomain, vdomain, nx, nv), (f, g))
    n, u, T = _maxwellian_moments(f, xdomain)

    # the Maxwellian on the cell centres, where the histograms hold their values; the two
    # histograms share their cells
    xedges, vedges = first(histograms)
    xc = (xedges[1:(end - 1)] .+ xedges[2:end]) ./ 2
    vc = (vedges[1:(end - 1)] .+ vedges[2:end]) ./ 2
    M = [n / sqrt(2 * DT(π) * T) * exp(-(v - u)^2 / (2 * T)) for _ in xc, v in vc]

    fig = Figure()
    for (column, (_, _, F), title) in zip(1:2, histograms, ("initial", "final"))
        ax = Axis(fig[1, column]; xlabel = "x", ylabel = "v", title = title)
        heatmap!(ax, xedges, vedges, F)
        contour!(ax, xc, vc, M; color = :white)
    end
    return fig
end

# The error of two series of different lengths, which names the shorter one: `a` and `b` pair a
# name with a length, and `what` names the entries that are counted.
function _mismatch(a::Pair, b::Pair, what)
    (short, nshort), (long, nlong) = last(a) < last(b) ? (a, b) : (b, a)
    return DimensionMismatch(
        "`$short` is shorter than `$long`: it holds $nshort $what, and `$long` holds $nlong")
end

# The arguments of a collision figure's series, checked before anything is drawn: each pairs its
# name with its values, an empty one holds nothing to draw, and one of another length than the
# first would pair its values with the wrong abscissae.
function _check_lengths(named::Pair...)
    for (name, y) in named
        isempty(y) && throw(ArgumentError("`$name` holds no value"))
    end
    first_name, first_values = first(named)
    for (name, y) in named
        length(y) == length(first_values) ||
            throw(_mismatch(name => length(y), first_name => length(first_values), "values"))
    end
    return nothing
end

# The pairs of `x` and `y` that have a place on the axes: a pair with a non-finite value is
# dropped. `xname` and `yname` name the two in the error of a series that keeps no pair.
function _finite_pairs(x, y, xname, yname)
    keep = isfinite.(x) .& isfinite.(y)
    any(keep) || throw(ArgumentError("no pair of `$xname` and `$yname` is finite"))
    return x[keep], y[keep]
end

# The values of `y` on a logarithmic axis. `log10` has no value at a negative number, nor at the
# zero of a series of one value, so both are refused.
function _log_values(y, name)
    any(<(0), y) && throw(ArgumentError(
        "`$name` holds a negative value, which a logarithmic axis cannot show"))
    if length(y) == 1 && iszero(only(y))
        throw(ArgumentError(
            "`$name` holds the single value zero, which a logarithmic axis cannot show"))
    end
    return float.(y)
end

# The series `ys` that share one logarithmic axis, each as the mask of its pairs that are drawn
# and the values they are drawn at. A zero has no logarithm: in a series with a positive value it
# is dropped, as a non-finite value is, such as a drift at `t = 0`. A series that is zero at every
# pair, such as the drift of an exactly conserved quantity, is drawn one decade below the smallest
# positive value on the axis, or at `floatmin` of its element type where the axis holds none.
function _log_axis(ys)
    positive = (v for y in ys for v in y if v > 0)
    floor = isempty(positive) ? nothing : minimum(positive) / 10
    return map(ys) do y
        any(>(0), y) && return (y .> 0, y)
        T = eltype(y)
        level = isnothing(floor) ? floatmin(T) : max(T(floor), floatmin(T))
        return (trues(length(y)), fill(level, length(y)))
    end
end

# Makie widens a logarithmic range of zero width, every value at one `c`, to the decades between
# `c²` and one. The square underflows to zero for a `c` below about `1e-162`, such as the
# `floatmin` of an axis whose every series is zero, and overflows for a `c` above about `1e154`;
# Makie then refuses the limits when the figure is saved. `setlims!` gives such an axis one decade on each
# side of `c` instead, within the positive finite `Float64`s. `plotted` holds the values of each
# series on the axis.
function _flat_log_limits!(setlims!, ax, plotted)
    c = Float64(first(first(plotted)))
    all(y -> all(==(c), y), plotted) || return nothing
    setlims!(ax, max(c / 10, floatmin(Float64)), min(c * 10, floatmax(Float64)))
    return nothing
end

# One series on `ax`: a line, or a marker where the series holds one value, which a line cannot
# show.
function _series!(ax, x, y; kwargs...)
    length(x) == 1 ? scatter!(ax, x, y; kwargs...) : lines!(ax, x, y; kwargs...)
    return nothing
end

function VlasovMethods.plot_drift(t::AbstractVector, quantities;
        labels = ("mass", "momentum", "energy"))
    isempty(quantities) && throw(ArgumentError("`quantities` holds no quantity"))
    length(labels) == length(quantities) ||
        throw(_mismatch(
            "labels" => length(labels), "quantities" => length(quantities), "entries"))
    _check_lengths("t" => t, (string(l) => Q for (l, Q) in zip(labels, quantities))...)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "relative drift", yscale = log10)
    series = map(zip(labels, quantities)) do (label, Q)
        Q₀ = first(Q)
        isfinite(Q₀) && !iszero(Q₀) || throw(ArgumentError(
            "the initial value of `$label` is $Q₀, so its relative drift is not defined"))
        x, drift = _finite_pairs(t, abs.(Q .- Q₀) ./ abs(Q₀), "t", label)
        (string(label), x, _log_values(drift, label))
    end
    plotted = map(series, _log_axis(map(last, series))) do (label, x, _), (keep, y)
        _series!(ax, x[keep], y[keep]; label = label)
        y[keep]
    end
    _flat_log_limits!(ylims!, ax, plotted)
    axislegend(ax)
    return fig
end

function VlasovMethods.plot_entropy(t::AbstractVector, S::AbstractVector)
    _check_lengths("t" => t, "S" => S)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "S_h")
    _series!(ax, _finite_pairs(t, S, "t", "S")...)
    return fig
end

function VlasovMethods.plot_error(t::AbstractVector, err::AbstractVector,
        closed::AbstractVector)
    _check_lengths("t" => t, "err" => err, "closed" => closed)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "error", yscale = log10)
    series = (("err", err, "error"), ("closed", closed, "closed form"))
    checked = map(series) do (name, values, label)
        x, finite = _finite_pairs(t, values, "t", name)
        (label, x, _log_values(finite, name))
    end
    plotted = map(checked, _log_axis(map(last, checked))) do (label, x, _), (keep, y)
        _series!(ax, x[keep], y[keep]; label = label)
        y[keep]
    end
    _flat_log_limits!(ylims!, ax, plotted)
    axislegend(ax)
    return fig
end

function VlasovMethods.plot_convergence(h::AbstractVector, err::AbstractVector)
    _check_lengths("h" => h, "err" => err)
    x, y = _finite_pairs(h, err, "h", "err")
    # the values on the logarithmic axes, before the fit, so that a negative value is named as
    # such and not as a pair that the fit cannot use
    xs, ys = _log_values(x, "h"), _log_values(y, "err")

    # the least-squares line through the pairs whose logarithms are finite; a pair with a zero
    # has no logarithm, and is neither drawn nor fitted
    fitted = findall((xs .> 0) .& (ys .> 0))
    length(fitted) ≥ 2 || throw(ArgumentError(
        "fewer than two pairs of positive finite `h` and `err`, so no slope can be fitted"))
    lnh, lnerr = log.(x[fitted]), log.(y[fitted])
    # a test of the deviations from the mean cannot see equal widths: the rounded mean of three
    # or more equal values can miss each of them by an ulp
    allequal(lnh) &&
        throw(ArgumentError("every `h` is the same, so no slope can be fitted"))
    dh = lnh .- sum(lnh) / length(lnh)
    Shh = sum(abs2, dh)
    p = sum(dh .* lnerr) / Shh
    C = exp(sum(lnerr) / length(lnerr) - p * sum(lnh) / length(lnh))

    fig = Figure()
    ylabel = "error ∝ h^p, p = $(round(p; sigdigits = 3))"
    ax = Axis(fig[1, 1]; xlabel = "h", ylabel = ylabel, xscale = log10, yscale = log10)
    scatter!(ax, xs[fitted], ys[fitted])
    hs = sort(x[fitted])
    lines!(ax, hs, C .* hs .^ p)
    # the widths differ, so only the `y` axis can be flat: an error that is the same at each
    # width, which the fitted line of slope zero meets within rounding
    _flat_log_limits!(ylims!, ax, (ys[fitted],))
    return fig
end

function VlasovMethods.plot_cost(cost::AbstractVector, err::AbstractVector)
    _check_lengths("cost" => cost, "err" => err)
    x, y = _finite_pairs(cost, err, "cost", "err")

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "wall-clock time", ylabel = "error",
        xscale = log10, yscale = log10)
    (xkeep, xs), = _log_axis((_log_values(x, "cost"),))
    (ykeep, ys), = _log_axis((_log_values(y, "err"),))
    keep = xkeep .& ykeep
    any(keep) || throw(ArgumentError(
        "no run has a positive `cost` and a positive `err`, so none has a place on both axes"))
    scatter!(ax, xs[keep], ys[keep])
    _flat_log_limits!(xlims!, ax, (xs[keep],))
    _flat_log_limits!(ylims!, ax, (ys[keep],))
    return fig
end

function VlasovMethods.plot_groups(t::AbstractVector, groups; xlabel = "t", ylabel = "",
        yscale = identity, colors = (:blue, :orange, :green))
    isempty(groups) && throw(ArgumentError("no group to plot"))
    length(colors) ≥ length(groups) || throw(DimensionMismatch(
        "the colours hold $(length(colors)) colours; the groups hold $(length(groups)) groups"))
    for (_, series) in groups
        _check_series(t, series, nothing)
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, yscale = yscale)
    # every line of a group is drawn half transparent, so that the lines of a sweep read as a band
    shades = [(colour, 0.5) for (colour, _) in zip(colors, groups)]
    for ((_, series), shade) in zip(groups, shades)
        for y in series
            lines!(ax, t, y; color = shade)
        end
    end
    axislegend(ax, [LineElement(color = shade) for shade in shades],
        [name for (name, _) in groups])
    return fig
end

function VlasovMethods.plot_damping(runs, reference; title = "",
        colors = (:steelblue, :firebrick))
    isempty(runs) && throw(ArgumentError("no run to plot"))
    length(colors) ≥ length(runs) || throw(DimensionMismatch(
        "the colours hold $(length(colors)) colours; the runs hold $(length(runs)) runs"))
    for run in runs
        _check_lengths("t" => run.t, "y" => run.y)
        length(run.peak_y) == length(run.peak_t) || throw(_mismatch(
            "peak_y" => length(run.peak_y), "peak_t" => length(run.peak_t), "values"))
    end
    _check_lengths("t" => reference.t, "y" => reference.y)

    fig = Figure(size = (1200, 500))
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "log|E₁|", title = title)
    for (run, colour) in zip(runs, colors)
        lines!(ax, run.t, run.y; color = colour, label = run.label)
        scatter!(ax, run.peak_t, run.peak_y; color = colour, markersize = 8)
    end
    lines!(ax, reference.t, reference.y; color = :black, linestyle = :dash,
        label = reference.label)
    axislegend(ax; position = :lb)
    return fig
end

function VlasovMethods.plot_phase_space(z::AbstractArray{<:Real, 3}, slices; xdomain, titles)
    size(z, 1) == 2 || throw(ArgumentError(
        "`z` holds $(size(z, 1)) rows; a 1D1V trajectory holds two, `x` and `v`"))
    size(z, 2) ≥ 1 || throw(ArgumentError("`z` holds no particle"))
    size(z, 3) ≥ 1 || throw(ArgumentError("`z` holds no time slice"))
    isempty(slices) && throw(ArgumentError("`slices` holds no slice"))
    for slice in slices
        1 ≤ slice ≤ size(z, 3) || throw(ArgumentError(
            "the time slices of `z` are 1:$(size(z, 3)); got $slice"))
    end
    length(titles) == length(slices) || throw(DimensionMismatch(
        "the titles hold $(length(titles)) titles; the slices hold $(length(slices)) slices"))

    # each particle keeps the colour of the sign of its first velocity, so that each half of the
    # cloud can be followed through the panels
    colours = ifelse.(z[2, :, 1] .>= 0, :steelblue, :firebrick)
    vmax = ceil(maximum(abs, z[2, :, 1]))
    limits = ((first(xdomain), last(xdomain)), (-vmax, vmax))

    fig = Figure(size = (1500, 400))
    for (panel, (slice, title)) in enumerate(zip(slices, titles))
        ylabel = panel == 1 ? "v" : ""
        ax = Axis(fig[1, panel]; xlabel = "x", ylabel, title, limits)
        x = [_wrap(x₀, xdomain) for x₀ in z[1, :, slice]]
        scatter!(ax, x, z[2, :, slice]; color = colours, markersize = 2)
    end
    return fig
end

end
