# The diagnostics of a sweep: the norms of an HFM state against a ROM state, the relative
# field-energy error, the energy drift, the save indices of a sweep, and the growth rate of a
# field-energy history. Each takes the arrays themselves rather than a `Snapshots`, so a model
# whose data is no snapshots uses the same functions.

"The root sum of squares of `D` along its first axis, with the trailing axes kept. Unlike `norm`,
it does not rescale, so a difference above about `sqrt(floatmax(T))` overflows to `Inf`."
_norm_reduce(D::AbstractArray) = dropdims(sqrt.(sum(abs2, D; dims = 1)); dims = 1)

"The largest absolute value of `D` along its first axis, with the trailing axes kept."
_maxabs_reduce(D::AbstractArray) = dropdims(maximum(abs, D; dims = 1); dims = 1)

"""
    trajectory_errors(Xh, Vh, Xr, Vr)

The errors of a full-order state `(Xh, Vh)` against a reduced state `(Xr, Vr)`: the `NamedTuple`
`(x2, xmax, v2, vmax)` of the position and the velocity errors. `x2` is the root sum of squares of
the position difference over the particle axis, divided by the particle count, and `xmax` is the
largest absolute position difference, the ``L^\\infty`` error; `v2` and `vmax` are the same for
the velocities.

Each argument holds one saved state per `(ts, p)`, as a `(nₚ, nₛ, nparam)` array or as a
`(nₚ, nₛ)` matrix, and each result has the trailing shape, `(nₛ, nparam)` or `(nₛ,)`, of the
element type of its argument. `NaN` and `Inf` in the input propagate into the result.

The two arrays of one pair must have the same size; a difference raises a `DimensionMismatch`
naming both.
"""
function trajectory_errors(Xh::AbstractArray, Vh::AbstractArray,
        Xr::AbstractArray, Vr::AbstractArray)
    size(Xh) == size(Xr) || throw(DimensionMismatch(
        "the positions hold $(size(Xh)) and $(size(Xr))"))
    size(Vh) == size(Vr) || throw(DimensionMismatch(
        "the velocities hold $(size(Vh)) and $(size(Vr))"))

    nₚ = size(Xh, 1)
    Dx = Xh .- Xr
    Dv = Vh .- Vr
    return (x2 = _norm_reduce(Dx) ./ nₚ, xmax = _maxabs_reduce(Dx),
        v2 = _norm_reduce(Dv) ./ nₚ, vmax = _maxabs_reduce(Dv))
end

"""
    potential_error(Φh, Φr)

The error of a reduced potential against a full-order potential: the root sum of squares of
`Φh - Φr` over the coefficient axis, per saved time and parameter. It is the absolute
``L^2``-type error of the potential, with no divisor.

`Φh` and `Φr` hold the coefficients of one potential per `(ts, p)`, as `(nh, nₛ, nparam)` arrays,
and the result is `(nₛ, nparam)` of the element type of the input. `NaN` and `Inf` propagate.

The arrays must have the same size; a difference raises a `DimensionMismatch` naming both.
"""
function potential_error(Φh::AbstractArray, Φr::AbstractArray)
    size(Φh) == size(Φr) || throw(DimensionMismatch(
        "the potential coefficients hold $(size(Φh)) and $(size(Φr))"))

    return _norm_reduce(Φh .- Φr)
end

"""
    field_energy_error(Wh, Wr)

The relative error of a reduced field energy against a full-order one, `abs.(Wh - Wr) ./ Wh`, per
saved time and parameter.

A field energy is `½ϕᵀKϕ ≥ 0` by construction, so the division is by a non-negative number and no
guard is needed: a `Wh` of `0` gives `Inf`, or `NaN` where `Wr` is also `0`, and both propagate.

`Wh` and `Wr` must have the same size; a difference raises a `DimensionMismatch` naming both.
"""
function field_energy_error(Wh::AbstractArray, Wr::AbstractArray)
    size(Wh) == size(Wr) || throw(DimensionMismatch(
        "the field energies hold $(size(Wh)) and $(size(Wr))"))

    return abs.(Wh .- Wr) ./ Wh
end

"""
    energy_drift(H)

The relative drift of an energy history `H`: `(H .- H[1]) ./ H[1]`, the energy at each saved time
against the initial one.

The energies of a conservative model are constant, so the drift is its discretisation error and
`H[1]` is its reference. `NaN` and `Inf` in `H` propagate into the result.

# Throws

- `ArgumentError` when `H` is empty, and when `H[1]` is zero, because a relative drift needs a
  non-zero reference.
"""
function energy_drift(H::AbstractVector)
    isempty(H) && throw(ArgumentError("an energy history holds at least one energy"))
    H[1] == 0 && throw(ArgumentError("the initial energy is zero"))
    return (H .- H[1]) ./ H[1]
end

"""
    saved_indices(nₛ)

The three saved states of a sweep of `nₛ` saved times that a plot of one trend uses: the final one
`nₛ`, the one nearest three quarters of the way, `div(3(nₛ - 1), 4) + 1`, and the one nearest
halfway, `div(nₛ - 1, 2) + 1`. The remainders take the floor, so the index never reaches past
`nₛ`. For `nₛ ≤ 3` the last two indices are the same state: `saved_indices(3) == (3, 2, 2)`.

# Throws

- `ArgumentError` when `nₛ < 2`: a sweep saves the initial state and at least one more, and a
  single saved time has no trend to plot.
"""
function saved_indices(nₛ::Integer)
    nₛ ≥ 2 || throw(ArgumentError(
        "a sweep saves the initial state and at least one more, got nₛ = $nₛ"))
    return (nₛ, div(3 * (nₛ - 1), 4) + 1, div(nₛ - 1, 2) + 1)
end

# The indices of the local maxima of `w`: a point no lower than the `window` points on each side.
# A point equal to the sample on its left is no maximum of its own, so a flat plateau counts once,
# at its first index, and the fit does not weight the plateau twice. An equal height elsewhere in
# the window does not drop a peak, so two equal peaks that are apart are two maxima.
function _find_maxima(w::AbstractVector, window::Int)
    inds = Int[]
    for i in (1 + window):(length(w) - window)
        left = @view w[(i - window):(i - 1)]
        right = @view w[(i + 1):(i + window)]
        if w[i] ≥ maximum(left) && w[i] ≥ maximum(right) && w[i - 1] != w[i]
            push!(inds, i)
        end
    end
    return inds
end

# The slope of the least-squares lines through the points `(t[ind], log(w[ind]))`, in closed form:
# `Σ(t - t̄)(y - ȳ) / Σ(t - t̄)²`. The arithmetic is in `eltype(w)`, so a `Float32` history gives a
# `Float32` slope.
function _log_slope(t::AbstractVector, w::AbstractVector{T}, inds) where {T}
    τ = T.(t[inds])
    d = τ .- sum(τ) / length(τ)
    y = log.(w[inds])
    return dot(d, y .- sum(y) / length(y)) / dot(d, d)
end

"""
    growth_rate(t, W; window = 2, maxima = 2:5)

The growth rate ``γ`` of the field energy history `W` at the times `t`: half the least-squares
slope of ``\\log W`` over the selected local maxima, because ``W ∝ |E|² ∝ e^{2γt}``. The times
are physical time, so this is the physical growth rate. The notebook this package's reduced-order
work follows plotted that slope itself under the name ``γ``, so the rates here are half the figures
it drew.

A maximum is a point no lower than the `window` points on each side, and a flat maximum counts
once, at its first index. `maxima` selects which of them the fit uses, by their order along `t`,
which drops the first few maxima; the fit runs in `eltype(W)`, so it takes no dependency on
`Statistics`.

A matrix `W` gives one `γ` per column, as a `Vector` of the element type of `W`.

# Throws

- `ArgumentError` when `window < 1`; when `maxima` names an index twice, naming it, because one
  maximum may not enter the fit twice; when it selects fewer than two maxima, or reaches past the
  maxima found, each naming the number found; when it selects an index below `1`; and when `W` is
  not positive at a selected maximum, where ``\\log W`` is not defined.
- `DimensionMismatch` when `t` and `W` differ in length, naming both.
"""
function growth_rate(t::AbstractVector, W::AbstractVector;
        window::Int = 2, maxima = 2:5)
    length(t) == length(W) || throw(DimensionMismatch(
        "the times hold $(length(t)) values; the field energies hold $(length(W))"))
    window ≥ 1 || throw(ArgumentError("a maximum needs a window of 1 or more, got $window"))
    inds = _find_maxima(W, window)

    if !allunique(maxima)
        repeated = first(m for m in maxima if count(==(m), maxima) > 1)
        throw(ArgumentError("the maxima range repeats the index $repeated"))
    end
    length(maxima) ≥ 2 || throw(ArgumentError(
        "a regression needs two maxima, got $(length(maxima)) of the " *
        "$(length(inds)) maxima found"))
    maximum(maxima) ≤ length(inds) || throw(ArgumentError(
        "the maxima range selects index $(maximum(maxima)), beyond the " *
        "$(length(inds)) maxima found"))
    minimum(maxima) ≥ 1 || throw(ArgumentError(
        "the maxima range selects index $(minimum(maxima)); the maxima are counted from 1"))

    sel = inds[maxima]
    any(≤(0), @view W[sel]) && throw(ArgumentError(
        "the field energy is not positive at the selected maxima, got $(view(W, sel))"))

    return _log_slope(t, W, sel) / 2
end

function growth_rate(t::AbstractVector, W::AbstractMatrix; kwargs...)
    return [growth_rate(t, w; kwargs...) for w in eachcol(W)]
end

"""
    best_and_worst(γh, γr)

The indices of the equation parameter whose reduced growth rate is nearest to, and furthest from,
its full-order value: the `argmin` and the `argmax` of `abs.(γr .- γh)`. Each is the first index
of a tie.
"""
function best_and_worst(γh::AbstractVector, γr::AbstractVector)
    d = abs.(γr .- γh)
    return argmin(d), argmax(d)
end
