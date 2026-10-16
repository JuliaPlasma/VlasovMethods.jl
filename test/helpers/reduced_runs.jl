# The helpers of the reduced-form runs: the timespan of a chunk, the solver tolerance of the
# Lobatto IIIA–IIIB(2) runs, and the comparator that finds the first save at which a drift exceeds
# its bound. It is included by `test/particles/reduced_equations.jl` and
# `test/verification/reduced_equations.jl`, and is listed in no `runtests.jl`.

"""
    chunk_span(::Type{T}, nsave, dt)

The timespan of one chunk of `nsave` steps of `dt`, taken from zero: `T(nsave)·dt`.

The integrator takes `nsave` steps over this span, while the problem's own `ntime` can count one
more, so a test asserts `ntime(sol)` and never `ntime(prob)`. The analytic model is autonomous, so
every chunk restarts the clock at zero. A span built from the run's own clock, `(t_{j-1}, t_j)`, is
not exactly `nsave·dt` long, and its solution takes `nsave - 1` steps at some `j`.
"""
chunk_span(::Type{T}, nsave, dt::T) where {T} = T(nsave) * dt

"""
    lobatto_abstol(T, k)

The solver tolerance of a Lobatto IIIA–IIIB(2) run in a basis of size `k`: `1e-14` in `Float64` and
`max(8, 2k)·eps(T)` otherwise, so that no run depends on the default of the integrator package.
"""
lobatto_abstol(::Type{Float64}, k) = 1e-14
lobatto_abstol(::Type{T}, k) where {T} = max(8, 2k) * eps(T)

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
