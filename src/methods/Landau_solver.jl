# The guess of the first three steps: a midpoint extrapolation of the vector field from
# `(t - Δt, vn)` to `t`, written into `v`. It allocates. The `EquationProblem` constructor of
# `GeometricEquations` sets the first type parameter of the problem with `eval`, so the
# assertion `::GeometricEquations.ODEProblem` states it, and the `extrapolate!` call dispatches
# statically.
function _landau_startup_guess!(v, vn, t, Δt, params, landau)
    problemGNI = GeometricEquations.ODEProblem(
        (v̇, t, v, params) -> collisional_vectorfield!(v̇, v, params, landau),
        (t, t+Δt), Δt, vn; parameters = params)::GeometricEquations.ODEProblem
    extrapolate!(t - Δt, vn, t, v, problemGNI, MidpointExtrapolation(5))
end

"""
    Picard_iterate_Landau_nls!(landau, tol, ftol, β, Δt, ti, t, v_prev, v_prev_2, rhs_prev, m,
        chunksize)

Take one implicit-midpoint step of the Landau model by five Picard iterations from the
velocities `v_prev`, write the last iterate into the particles of `landau.dist`, and return the
iterates and their midpoint vector fields as the histories `(v = …, v̇ = …)`, six of each.

`ti` counts the steps from 1, and `t` is the time of the guess. From the fourth step on, the
guess is the Hermite extrapolation from `v_prev_2` and `v_prev` with the vector fields in
`rhs_prev[:, :, 2]` and `rhs_prev[:, :, 1]`; the first three steps (`ti < 4`) take it from a
midpoint extrapolation of the vector field. On return, `rhs_prev[:, :, 1]` holds the last
midpoint vector field and `rhs_prev[:, :, 2]` the one it held before. `tol`, `ftol`, `β`, `m`
and `chunksize` are not read.

The residual ``|v - (v_{prev} + Δt \\, \\dot{v})|`` of the guess and of each iterate is sent as
a `@debug` message; nothing is written to `stdout`.

# Allocations

From the fourth step on, a call allocates at most the bytes of its six vector-field
evaluations, plus the `Base.summarysize` of the histories it returns, plus 1 KiB. A
vector-field evaluation allocates no more than its two least-squares solves `K1 \\ LJ` and
`K2 \\ LJ` where the number of particles is not the number of basis functions (Np ≠ M); for a
square `K` (Np == M) it takes the pivoted QR, where `\\` takes LU, and that bound does not
hold. The guess of the first three steps allocates at most the bytes of its 37
vector-field evaluations plus 200 kB, for the arrays of GeometricIntegrators' extrapolation
and the problem it is given. `test/integration/landau_solver.jl` asserts both bounds.
"""
function Picard_iterate_Landau_nls!(
        landau, tol, ftol, β, Δt, ti, t, v_prev, v_prev_2, rhs_prev, m, chunksize)
    # β, the damping parameter of a damped Picard iteration (β = 1 is plain Picard), is not read
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
        extrapolate!(
            t - 2Δt, v_prev_2, view(rhs_prev, :, :, 2), t - Δt, v_prev,
            view(rhs_prev, :, :, 1), t, v_guess, HermiteExtrapolation())
    else
        _landau_startup_guess!(v_guess, v_prev, t, Δt, params, landau)
    end

    v_midpoint = cache.v
    v̇_midpoint = cache.v̇

    v_midpoint .= (v_guess .+ v_prev) ./ 2
    collisional_vectorfield!(v̇_midpoint, v_midpoint, params, landau)
    @debug "Landau Picard residual" iteration=0 residual=norm(
        v_guess .- (v_prev .+ Δt .* v̇_midpoint))

    v_sol = Vector{typeof(v_guess)}(undef, 6)
    v̇_sol = Vector{typeof(v̇_midpoint)}(undef, 6)
    v_sol[1] = copy(v_guess)
    v̇_sol[1] = copy(v̇_midpoint)

    for i in 1:5
        v_guess .= v_prev .+ Δt .* v̇_midpoint

        v_midpoint .= (v_guess .+ v_prev) ./ 2

        collisional_vectorfield!(v̇_midpoint, v_midpoint, params, landau)

        @debug "Landau Picard residual" iteration=i residual=norm(
            v_guess .- (v_prev .+ Δt .* v̇_midpoint))

        v_sol[i + 1] = copy(v_guess)
        v̇_sol[i + 1] = copy(v̇_midpoint)
    end

    sol = (v = v_sol, v̇ = v̇_sol)

    # update solution array
    dist.particles.v .= v_guess

    # update rhs storage, through views: `rhs_prev[:, :, 1]` on the right would be a copy
    view(rhs_prev, :, :, 2) .= view(rhs_prev, :, :, 1)
    view(rhs_prev, :, :, 1) .= v̇_midpoint

    # return solution at t
    return sol
end
