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

function VlasovMethods.plot_distribution(dist::ParticleDistribution{DT, 1, 1};
        xdomain, vdomain, nx::Integer = 100, nv::Integer = 100) where {DT}
    # `isempty` iterates a `ParticleList` and raises `BoundsError` on an empty list, so both
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
        # the wrap cannot place a non-finite position, and no cell holds a `v` outside the domain
        isfinite(x) && va ≤ v < vb || continue
        # `mod(-1e-20, 2π) == 2π` in `Float64`, so the last cell needs the clamp.
        i = min(floor(Int, (_wrap(x, xdomain) - xa) / hx) + 1, nx)
        # `(v - va)/hv` can round up to `nv` for `v = prevfloat(vb)`, so the last cell needs the
        # clamp as well.
        j = min(floor(Int, (v - va) / hv) + 1, nv)
        F[i, j] += dist.particles.w[1, k] / (hx * hv)
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "v")
    heatmap!(ax, range(xa, xb; length = nx + 1), range(va, vb; length = nv + 1), F)
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

# The arguments of `plot_series`, checked before anything is drawn: an empty collection has no
# line to draw, an empty line holds no value, a line of another length would pair its values with
# the wrong abscissae, and a label count that differs from the series count would name the wrong
# lines.
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
        labels = nothing, xlabel = "t", ylabel = "")
    _check_series(t, series, labels)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel)
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

end
