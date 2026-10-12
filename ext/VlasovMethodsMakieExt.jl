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
    # `isempty` of a `ParticleList` iterates it, and iterating an empty list raises a
    # `BoundsError` in ParticleMethods 0.1.2, so the guards of both particle methods ask the list
    # for its length instead.
    iszero(length(dist.particles)) && throw(ArgumentError("no particles to plot"))

    xa, xb = first(xdomain), last(xdomain)
    va, vb = first(vdomain), last(vdomain)
    hx = (xb - xa) / nx
    hv = (vb - va) / nv

    F = zeros(DT, nx, nv)
    for k in eachindex(dist.particles)
        v = dist.particles.v[1, k]
        va ≤ v < vb || continue
        # `mod(-1e-20, 2π) == 2π` in `Float64`, so the last cell needs the clamp.
        i = min(floor(Int, (_wrap(dist.particles.x[1, k], xdomain) - xa) / hx) + 1, nx)
        # `(v - va)/hv` rounds up to `nv` at the upper end, as `prevfloat(vb)` does, so the last
        # cell needs the clamp as well.
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

end
