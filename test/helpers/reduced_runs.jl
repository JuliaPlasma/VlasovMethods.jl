# The helpers of the `T = 1000` runs: the timespan of a chunk, the number of steps the integrator
# takes over it, and the comparator that finds the first save at which a drift exceeds its bound. It
# is included by `test/particles/reduced_equations.jl` — the check-4 and check-5 file, whose problems
# span one chunk — and by `test/verification/reduced_equations.jl`, which is check 3's run. It is
# listed in no `runtests.jl` and exports nothing.

"""
    chunk_span(::Type{T}, nsave, dt)

The timespan of one chunk of `nsave` steps of `dt`, taken **from zero**: `T(nsave)·dt`, with no
`prevfloat` and no precision branch.

The count of a problem and the count of its solution are not the same number, and the solution's is
the steps taken. GeometricEquations counts `div(tend - tbegin, tstep, RoundUp)` for a problem
(`problems/equation_problem.jl:138-144`), while the integrator steps over the solution's
`TimeSeries`, `tbegin:Δt:tend`, whose step count is `length - 1` (GeometricSolutions
`geometric_solution.jl:77`, `timeseries.jl:5-13,48-49`). Measured with a real `integrate`: the span
`(zero(T), T(1000)·T(1e-4))` takes **1000** steps in `Float64` and in `Float32`, ending at the span's
end, against an `ntime(prob)` of `1001` in `Float64`, and `nsave` steps for each of
`nsave ∈ {10, 100, 999, 1000, 1001, 5000}` in both; `prevfloat` of that span reports
`ntime(prob) = 1000` and takes **999** in both. So every problem this part integrates asserts
`ntime(sol)`, and `ntime(prob)` is asserted nowhere.

Taken from zero, the span is the same value every chunk, and the analytic model's equation is
autonomous — neither the PODE nor the ODE form reads `t` — so restarting the clock loses nothing. A
chunk built from the run's own clock, the float span `(t_{j-1}, t_j)`, is worse than either: the
difference of two consecutive times is not exactly `nsave·dt`, so its solution takes `nsave` steps at
some `j` and `nsave - 1` at others, and a chunk that falls one step short puts every later save at
the wrong time.
"""
chunk_span(::Type{T}, nsave, dt::T) where {T} = T(nsave) * dt

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
