function Picard_iterate_Landau_nls!(
        landau, tol, ftol, β, Δt, ti, t, v_prev, v_prev_2, rhs_prev, m, chunksize)
    # β is the damping parameter for damped Picard iterations, with β = 1 yielding regular Picard iterations
    # ti is the time index at which v_new is being computed, i.e. for t = ti * Δt
    # v_prev is v at the previous timestep 
    # v_prev_2 is v at two timesteps prior to t
    # rhs_prev[:,:,1] is the rhs at t - Δt, and rhs_prev[:,:,2] is the rhs at t - 2Δt

    dist = landau.dist
    ent = landau.entropy

    # creating this to store the guess for the moment, for diagnostic purposes
    v_guess = copy(dist.particles.v)

    params = (dist = dist, ent = ent)

    # use Hermite extrapolation to get an initial guess
    if ti ≥ 4
        extrapolate!(
            t - 2Δt, v_prev_2, view(rhs_prev, :, :, 2), t - Δt, v_prev,
            view(rhs_prev, :, :, 1), t, v_guess, HermiteExtrapolation())
    else
        problemGNI = GeometricEquations.ODEProblem(
            (v̇, t, v, params) -> collisional_vectorfield!(v̇, v, params, landau),
            (t, t+Δt), Δt, v_prev; parameters = params)
        extrapolate!(
            t - Δt, v_prev, t, v_guess, problemGNI, MidpointExtrapolation(5))
    end

    v_midpoint = landau.cache[eltype(v_guess)].v
    v̇_midpoint = landau.cache[eltype(v_guess)].v̇

    v_midpoint .= (v_guess .+ v_prev) ./ 2
    collisional_vectorfield!(v̇_midpoint, v_midpoint, params, landau)
    println("   |f(v_0)| = ", norm(v_guess .- (v_prev .+ Δt .* v̇_midpoint)))

    v_sol = [copy(v_guess)]
    v̇_sol = [copy(v̇_midpoint)]

    for i in 1:5
        v_guess .= v_prev .+ Δt .* v̇_midpoint

        v_midpoint .= (v_guess .+ v_prev) ./ 2

        collisional_vectorfield!(v̇_midpoint, v_midpoint, params, landau)

        println("   |f(v_$i)| = ", norm(v_guess .- (v_prev .+ Δt .* v̇_midpoint)))

        push!(v_sol, copy(v_guess))
        push!(v̇_sol, copy(v̇_midpoint))
    end

    sol = (v = v_sol, v̇ = v̇_sol)

    println()

    # update solution array
    dist.particles.v .= v_guess

    # update rhs storage
    rhs_prev[:, :, 2] .= rhs_prev[:, :, 1]
    rhs_prev[:, :, 1] .= v̇_midpoint

    # return solution at t
    return sol
end
