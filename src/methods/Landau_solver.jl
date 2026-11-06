# The guess of the first three steps: a midpoint extrapolation of the vector field from
# `(t - Δt, vn)` to `t`, written into `v`. It allocates. The `EquationProblem` constructor of
# `GeometricEquations` sets the first type parameter of the problem with `eval`, so the
# assertion `::GeometricEquations.ODEProblem` states it, and the `extrapolate!` call dispatches
# statically.
function _landau_startup_guess!(v, vn, t, Δt, params, landau)
    problemGNI = GeometricEquations.ODEProblem(
        (v̇, t, v, params) -> collisional_vectorfield!(v̇, v, params, landau),
        # a copy, because the `StateVariable` of `GeometricBase` takes no view
        (t, t+Δt), Δt, copy(vn); parameters = params)::GeometricEquations.ODEProblem
    extrapolate!(t - Δt, vn, t, v, problemGNI, MidpointExtrapolation(5))
end

"""
    Picard_iterate_Landau_nls!(landau, Δt, ti, t, v_prev, v_prev_2, rhs_prev; abstol,
        max_iters = 100)

Take one implicit-midpoint step of the Landau model from the velocities `v_prev` by a Picard
iteration to the tolerance `abstol`, write the accepted velocities into the particles of
`landau.dist`, and return `(iterations = …, residual = …)`.

From the guess `v₀`, iteration `k = 0, 1, …` evaluates the vector field `f_k` at the midpoint
`(v_k + v_prev)/2` and the residual ``r_k = v_k - (v_{prev} + Δt \\, f_k)``. The first iterate
with ``\\|r_k\\|_∞ ≤`` `abstol` is accepted: the particles take `v_prev + Δt f_k`, so that the
accepted velocities and the stored vector field `f_k` belong together. Otherwise
`v_{k+1} = v_prev + Δt f_k`. `iterations` is `k + 1`, the number of vector-field evaluations,
and `residual` is ``\\|r_k\\|_∞``. The guess of the first three steps evaluates the vector field
too; those evaluations are not counted.

`abstol` must be finite and positive and `max_iters` at least 1, or an `ArgumentError` is thrown
before any vector field is evaluated. When no iterate meets `abstol` within `max_iters`
evaluations, an `ErrorException` is thrown that names the number of evaluations and the last
residual, and neither the particles nor `rhs_prev` are written. A residual that is not finite
never meets `abstol`, and an `abstol` below the rounding error of the residual, about the
machine epsilon of the velocities' element type times their size, may never be met. An
exception of the vector field passes through: the projection throws a `DomainError` once an
iterate moves particles out of the velocity domain, which can happen before `max_iters`
evaluations.

`ti` counts the steps from 1, and `t` is the time of the new velocities. From the fourth step on,
the guess is the Hermite extrapolation from `v_prev_2` and `v_prev` with the vector fields in
`rhs_prev[:, :, 2]` and `rhs_prev[:, :, 1]`; the first three steps (`ti < 4`) take it from a
midpoint extrapolation of the vector field. On return, `rhs_prev[:, :, 1]` holds the accepted
vector field `f_k` and `rhs_prev[:, :, 2]` the one it held before.

The residual ``\\|r_k\\|_∞`` of each iterate is sent as the `@debug` message
"Landau Picard residual" with the keyword `iteration = k`; nothing is written to `stdout`.

# Allocations

From the fourth step on, a call allocates at most `iterations` times the bytes of one
vector-field evaluation, plus 1 KiB. A vector-field evaluation allocates no more than its two
least-squares solves `K1 \\ LJ` and `K2 \\ LJ` where the number of particles is not the number of
basis functions (Np ≠ M); for a square `K` (Np == M) it takes the pivoted QR, where `\\` takes
LU, and that bound does not hold. The guess of the first three steps allocates at most the bytes
of its 37 vector-field evaluations plus 200 kB, for the arrays of GeometricIntegrators'
extrapolation and the problem it is given. `test/integration/landau_solver.jl` asserts both
bounds.
"""
function Picard_iterate_Landau_nls!(
        landau, Δt, ti, t, v_prev, v_prev_2, rhs_prev; abstol, max_iters::Integer = 100)
    isfinite(abstol) && abstol > 0 ||
        throw(ArgumentError("abstol must be finite and positive, got $abstol"))
    max_iters ≥ 1 || throw(ArgumentError("max_iters must be at least 1, got $max_iters"))

    # ti is the time index at which v_new is being computed, i.e. for t = ti * Δt
    # v_prev is v at the previous timestep
    # v_prev_2 is v at two timesteps prior to t
    # rhs_prev[:,:,1] is the rhs at t - Δt, and rhs_prev[:,:,2] is the rhs at t - 2Δt

    dist = landau.dist
    ent = landau.entropy

    # the iterate, in the cache of the particles' element type; both extrapolations below
    # overwrite it
    cache = landau.cache[eltype(dist.particles.v)]
    v_guess = cache.guess

    params = (dist = dist, ent = ent)

    # use Hermite extrapolation to get an initial guess
    if ti ≥ 4
        # the extrapolation takes its three times in one type
        tp, Δtp = promote(t, Δt)
        extrapolate!(
            tp - 2Δtp, v_prev_2, view(rhs_prev, :, :, 2), tp - Δtp, v_prev,
            view(rhs_prev, :, :, 1), tp, v_guess, HermiteExtrapolation())
    else
        _landau_startup_guess!(v_guess, v_prev, t, Δt, params, landau)
    end

    v_midpoint = cache.v
    v̇_midpoint = cache.v̇

    residual = zero(eltype(v_midpoint))
    for k in 0:(max_iters - 1)
        v_midpoint .= (v_guess .+ v_prev) ./ 2
        collisional_vectorfield!(v̇_midpoint, v_midpoint, params, landau)
        # the residual, in `v_midpoint`, which the next iteration overwrites
        v_midpoint .= v_guess .- (v_prev .+ Δt .* v̇_midpoint)
        residual = maximum(abs, v_midpoint)
        @debug "Landau Picard residual" iteration=k residual

        if residual ≤ abstol
            # accept the update of the converged vector field, so that the state and the
            # stored derivative belong together
            dist.particles.v .= v_prev .+ Δt .* v̇_midpoint

            # update rhs storage, through views: `rhs_prev[:, :, 1]` on the right would be a
            # copy
            view(rhs_prev, :, :, 2) .= view(rhs_prev, :, :, 1)
            view(rhs_prev, :, :, 1) .= v̇_midpoint

            return (iterations = k + 1, residual = residual)
        end

        v_guess .= v_prev .+ Δt .* v̇_midpoint
    end

    error("the Landau Picard iteration did not converge in $max_iters vector-field ",
        "evaluations: the last residual is $residual, above abstol = $abstol")
end
