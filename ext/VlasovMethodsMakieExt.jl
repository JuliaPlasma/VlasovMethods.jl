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

# A series plot's pair of arguments, checked before anything is drawn: an empty collection has no
# line to draw, and a series of another length would pair its values with the wrong abscissae.
function _check_series(t, series)
    isempty(series) && throw(ArgumentError("no series to plot"))
    for y in series
        length(y) == length(t) || throw(DimensionMismatch(
            "a series holds $(length(y)) values; the abscissae hold $(length(t))"))
    end
    return nothing
end

# The one-axis figure of a series plot: one line per series, against the common abscissae.
function _series_figure(t, series; xlabel, ylabel)
    _check_series(t, series)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel)
    for y in series
        lines!(ax, t, y)
    end
    return fig
end

function VlasovMethods.plot_errors(t::AbstractVector, series)
    return _series_figure(t, series; xlabel = "t", ylabel = "error")
end

function VlasovMethods.plot_energy(t::AbstractVector, series)
    return _series_figure(t, series; xlabel = "t", ylabel = "energy")
end

function VlasovMethods.plot_growth_rate(t::AbstractVector, series)
    return _series_figure(t, series; xlabel = "parameter", ylabel = "γ")
end

end
