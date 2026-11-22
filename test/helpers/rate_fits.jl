# The rate fits of VlasovReductionBenchmark (`src/full_order.jl` at 3e4fa5e), ported unchanged:
# the least-squares slope of `y = log|E₁|` through its local maxima for a damped mode, and
# through every point of the window for a growing mode.

"The slope of the least-squares line through the points `(t, y)`."
function slope(t, y)
    tm, ym = sum(t) / length(t), sum(y) / length(y)

    return sum((t .- tm) .* (y .- ym)) / sum(abs2, t .- tm)
end

"""
    landau_rate_fit(t, y, window)

The named tuple `(γ, ω, npeaks)`: `γ` is the slope of the least-squares line through the local
maxima `y[n - 1] < y[n] ≥ y[n + 1]` of `y = log|E₁|` with `t[n]` in `window`, `ω` is `π` over
the mean spacing of those maxima, two per period of ``E_1``, and `npeaks` their number. It throws
an error for fewer than four maxima.
"""
function landau_rate_fit(t, y, window)
    peaks = [n
             for n in 2:(length(t) - 1)
             if window[1] ≤ t[n] ≤ window[2] && y[n - 1] < y[n] ≥ y[n + 1]]
    length(peaks) ≥ 4 ||
        error("the Landau fit found $(length(peaks)) maxima of log|E₁| in " *
              "[$(window[1]), $(window[2])]; it needs at least four")
    tp = t[peaks]

    return (γ = slope(tp, y[peaks]), ω = π / (sum(diff(tp)) / (length(tp) - 1)),
        npeaks = length(peaks))
end

"The slope of the least-squares line through `y = log|E₁|` at every time `t` in `window`."
function growth_rate_fit(t, y, window)
    w = [window[1] ≤ tn ≤ window[2] for tn in t]

    return slope(t[w], y[w])
end
