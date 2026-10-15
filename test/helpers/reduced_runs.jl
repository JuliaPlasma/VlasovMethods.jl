# The helpers of the `T = 1000` runs: the timespan of a chunk and the step-count contract that the
# emitted problems of P37 assert, and the comparator that finds the first save at which a drift
# exceeds its bound. It is included by `test/particles/reduced_equations.jl` — the check-4 file,
# whose problem spans one chunk — and by `test/verification/reduced_equations.jl`, which is check 3's
# run. It is listed in no `runtests.jl` and exports nothing.

"""
    chunk_span(::Type{T}, nsave, dt)

The timespan of one chunk of `nsave` steps of `dt`, taken **from zero** and rounded **down** by one
`eps`: the largest double below `nsave·dt`.

GeometricEquations counts `div(tend - tbegin, tstep, RoundUp)` steps (`equation_problem.jl:125`),
and `RoundUp` is what makes a chunk of this length need care. The double nearest `nsave·dt` sits
**above** `nsave·dt` — measured: the true `1000·1e-4` is `0.100000000000000004792` and the double
`0.1` is `0.100000000000000005551` — so the naive span counts `nsave + 1` steps in `Float64`, while
in `Float32` the nearest double sits below and it counts `nsave`. A chunk built from the run's own
clock, the float span `(t_{j-1}, t_j)`, is worse still: the difference of two consecutive times is
not exactly `nsave·dt`, and the count then comes out `nsave` at some `j` and `nsave + 1` at others.
One span, one count, for every `j` and both precisions, is what the run's save grid needs: a chunk
that overshoots by one step puts every later save at the wrong time.

Taken from zero, the span is the same value every chunk, and the analytic model's equation is
autonomous — neither the PODE nor the ODE form reads `t` — so restarting the clock loses nothing.

See [`chunk_steps`](@ref) for the count this span has.
"""
chunk_span(::Type{T}, nsave, dt::T) where {T} = prevfloat(T(nsave) * dt)

"""
    chunk_steps(::Type{T}, nsave, dt)

The number of steps GeometricEquations counts for [`chunk_span`](@ref): `div(span, dt, RoundUp)`,
which is exactly `nsave` in `Float64` and in `Float32` for `nsave = 1000` and `dt = 1e-4` — the
intended count of a chunk, and the count the emitted problems assert.
"""
chunk_steps(::Type{T}, nsave, dt::T) where {T} = div(chunk_span(T, nsave, dt), dt, RoundUp)

"""
    first_exceedance(r, bound)

The 1-based position of the first entry of `r` whose relative drift is not `≤ bound`, and `nothing`
if every entry is. The test is written `value ≤ bound || return i` and never `value > bound`, so
that a `NaN` counts as an exceedance instead of being dropped in silence: a state that has
overflowed, or a ratio of a `NaN` Hamiltonian, must end the run rather than pass it.
"""
function first_exceedance(r, bound)
    for (i, value) in enumerate(r)
        value ≤ bound || return i
    end
    return nothing
end
