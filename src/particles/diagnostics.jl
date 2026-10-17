# The diagnostics of a sweep: the norms of an HFM state against a ROM state, the relative
# field-energy error, the energy drift, the save indices of a sweep, and the growth rate of a
# field-energy history. Each takes the arrays themselves rather than a `Snapshots`, so a model
# whose data is no snapshots uses the same functions.

"The root sum of squares of each column of `D`, whose first axis is the reduced one, with the
trailing axes kept."
function _norm_reduce(D::AbstractArray)
    n = size(D, 1)
    A = reshape(D, n, :)
    return reshape([norm(@view A[:, j]) for j in axes(A, 2)], size(D)[2:end])
end

"The largest absolute value of each column of `D` along its first axis, with the trailing axes
kept."
function _maxabs_reduce(D::AbstractArray)
    n = size(D, 1)
    A = reshape(D, n, :)
    return reshape([maximum(abs, @view A[:, j]) for j in axes(A, 2)], size(D)[2:end])
end

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

The relative error of a reduced field energy against a full-order one, `abs.(W - Wr) ./ W`, per
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
`nₛ`.

# Throws

- `ArgumentError` when `nₛ < 2`: a sweep saves the initial state and at least one more, and a
  single saved time has no trend to plot.
"""
function saved_indices(nₛ::Integer)
    nₛ ≥ 2 || throw(ArgumentError(
        "a sweep saves the initial state and at least one more, got nₛ = $nₛ"))
    return (nₛ, div(3 * (nₛ - 1), 4) + 1, div(nₛ - 1, 2) + 1)
end

# The indices of the local maxima of `w`: a point that is above the `window` points on its left
# and no lower than the `window` points on its right. Only a point with a full window on both
# sides can be one, and a flat maximum counts once, at its first index, so that the fit does not
# weight it twice.
function _find_maxima(w::AbstractVector, window::Int)
    inds = Int[]
    for i in (1 + window):(length(w) - window)
        left = @view w[(i - window):(i - 1)]
        right = @view w[(i + 1):(i + window)]
        if w[i] > maximum(left) && w[i] ≥ maximum(right)
            push!(inds, i)
        end
    end
    return inds
end

# The slope of the least-squares lines through the points `(t[ind], log(w[ind]))`, in closed form:
# `Σ(t - t̄)(y - ȳ) / Σ(t - t̄)²`. The arithmetic is in `eltype(w)`, so a `Float32` history gives a
# `Float32` slope.
function _log_slope(t::AbstractVector, w::AbstractVector{T}, inds) where {T}
    τ = T[t[i] for i in inds]
    y = T[log(w[i]) for i in inds]
    τ̄ = sum(τ) / length(τ)
    ȳ = sum(y) / length(y)

    num = zero(T)
    den = zero(T)
    for j in eachindex(τ)
        d = τ[j] - τ̄
        num += d * (y[j] - ȳ)
        den += d^2
    end
    return num / den
end

"""
    growth_rate(t, W; window = 2, maxima = 2:5)

The growth rate ``γ`` of the field energy history `W` at the times `t`: half the least-squares
slope of ``\\log W`` over the selected local maxima, because ``W ∝ |E|² ∝ e^{2γt}``. The times
are physical time, so this is the physical growth rate.

A maximum is a point no lower than the `window` points on each side, and a flat maximum counts
once, at its first index. `maxima` selects which of them the fit uses, by their order along `t`,
which drops the first few maxima; the fit runs in `eltype(W)`, so it takes no dependency on
`Statistics`.

A matrix `W` gives one `γ` per column, as a `Vector` of the element type of `W`.

# Throws

- `ArgumentError` when `maxima` selects fewer than two maxima, when it reaches past the maxima
  found, naming that count, and when `W` is not positive at a selected maximum, where ``\\log W``
  is not defined.
"""
function growth_rate(t::AbstractVector, W::AbstractVector;
        window::Int = 2, maxima = 2:5)
    inds = _find_maxima(W, window)

    length(maxima) ≥ 2 || throw(ArgumentError(
        "a regression needs two maxima, got $(length(maxima))"))
    maximum(maxima) ≤ length(inds) || throw(ArgumentError(
        "the maxima range selects index $(maximum(maxima)), beyond the " *
        "$(length(inds)) maxima found"))

    sel = inds[maxima]
    any(≤(0), @view W[sel]) && throw(ArgumentError(
        "the field energy is not positive at the selected maxima, got $(view(W, sel))"))

    return _log_slope(t, W, sel) / 2
end

function growth_rate(t::AbstractVector, W::AbstractMatrix;
        window::Int = 2, maxima = 2:5)
    return [growth_rate(t, @view W[:, j]; window = window, maxima = maxima)
            for j in axes(W, 2)]
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
