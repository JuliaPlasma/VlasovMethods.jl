# function IM_rule!(g, x, y, f, Δt)
#     g .= x .- y .- Δt * f(0.5 * (y .+ x))
# end


function IM_rule!(g, x, y, f, Δt, params)
    g .= x .- y .- Δt * f(0.5 * (y .+ x), params)
end

function IM_update!(y_new, y_new_guess, yₙ, f, Δt, params)
    # println("define F!")
    F!(g, x) = IM_rule!(g, x, yₙ, f, Δt, params)
    obj = zero(yₙ)
    opt = Options(x_reltol = 1e-8, f_reltol = 1e-8, x_abstol = 1e-12, f_abstol = 1e-12)
    # println("construct NewtonSolver")
    nl = NewtonSolver(y_new_guess, obj, config = opt)
    # println("Newton solve")
    SimpleSolvers.solve!(y_new, F!, nl)

    return y_new
end

function explicit_update!(v_new, v, Δt, rhs)
    @. v_new = v + Δt * rhs
end

function Picard_iterate_over_particles(f::Function, dist, sdist, tol, Δt)
    # j is the Picard iteration index
    # create vectors to store the current iteration and the two previous iterations
    v_new = zero(dist.particles.v)
    v_prev = copy(dist.particles.v) #TODO: this should be updated to computing an initial guess using Hermite extrapolation
    
    err = 1.
    j = 0
    println("entering Picard loop")
    while err > tol
        println("j=",j)

        S = projection(v_prev, dist, sdist)
        # L = VlasovMethods.compute_L_Maxwellian(sdist)
        L = VlasovMethods.compute_L(sdist)

        params = (dist = dist, sdist = sdist, B = sdist.basis, L = L, v_array = v_prev)

        @time Threads.@threads for i in axes(v_prev,2)
            # @show i
            v_new[:,i] .= IM_update!(v_new[:,i], v_prev[:,i], dist.particles.v[:,i], f, Δt, params)
        end

        @show err = norm(v_new .- v_prev)^2
        v_prev .= v_new
        j += 1
    end

    dist.particles.v .= v_new
    return v_new
end

function Picard_iterate_Landau!(dist, sdist, tol, Δt, ti, t, v_prev_2, rhs_prev, sdist2, max_iters)
    # ti is the time index at which v_new is being computed, i.e. for t = ti * Δt
    # dist.particles.v is v at the previous timestep 
    # v_prev_2 is v at two timesteps prior to t
    # rhs_prev[:,:,1] is the rhs at t - Δt, and rhs_prev[:,:,2] is the rhs at t - 2Δt
    
    # create vectors to store the current and previous iterations
    v_new = zero(dist.particles.v)
    v_prev = copy(dist.particles.v) #TODO: this should be updated to computing an initial guess using Hermite extrapolation
    
    # creating this to store the guess for the moment, for diagnostic purposes
    v_guess = copy(dist.particles.v) 
    
    if ti ≥ 4
        Extrapolators.extrapolate!(t - 2Δt, v_prev_2, rhs_prev[:,:,2], t - Δt, dist.particles.v, rhs_prev[:,:,1], t, v_guess, Extrapolators.HermiteExtrapolation())
    end

    rhs_prev[:,:,2] .= rhs_prev[:,:,1]

    v_prev .= v_guess

    # j is the Picard iteration index
    j = 0
    # println("entering Picard loop")
    err = 1.
    while err > tol
        println("j=",j)
        v_midpoint = (v_prev .+ dist.particles.v) / 2
        # Landau_rhs_2!(rhs_prev[1,:,:])
        S = projection(v_midpoint, dist, sdist)

        # TODO: this should also be done using the rhs function, not implemented here. 
        # println("computing J")
        # @time J = VlasovMethods.compute_L(sdist)
        J = VlasovMethods.compute_L(sdist)

        # println("computing L_ij")
        # @time Lij = VlasovMethods.compute_L_ij(sdist)
        Lij = VlasovMethods.compute_L_ij(sdist)

        # println("computing K")
        # @time K1_plus, K2_plus = VlasovMethods.compute_K_plus(v_midpoint, dist, sdist)
        K1_plus, K2_plus = VlasovMethods.compute_K_plus(v_midpoint, dist, sdist)

        rhs_prev[1,:,1] .= K1_plus * Lij * J  
        rhs_prev[2,:,1] .= K2_plus * Lij * J  

        # println("integrating")
        # @time VlasovMethods.explicit_update!(view(v_new,1,:), dist.particles.v[1,:], Δt, rhs_prev[1,:,1])
        # @time VlasovMethods.explicit_update!(view(v_new,2,:), dist.particles.v[2,:], Δt, rhs_prev[1,:,2])

        VlasovMethods.explicit_update!(view(v_new,1,:), dist.particles.v[1,:], Δt, rhs_prev[1,:,1])
        VlasovMethods.explicit_update!(view(v_new,2,:), dist.particles.v[2,:], Δt, rhs_prev[2,:,1])

        err = norm(v_new .- v_prev)
        err_max = norm(v_new .- v_prev, Inf)

        println("L2 norm of (v_j+1 - v_j) = ", err)
        println("Linf norm of (v_j+1 - v_j) = ", err_max)
        println("")

        v_prev .= v_new

        if j == max_iters
            @warn " Maximum iterations reached, aborting Picard iterations and storing result"
            break
        end

        j += 1
    end

    # # updated rhs_prev to store the new rhs
    # rhs_prev[:,:,2] .= rhs_prev[:,:,1]
    # rhs_prev[1,:,1] .= rhs_1
    # rhs_prev[2,:,1] .= rhs_2

    # update solution array
    dist.particles.v .= v_new

    @show v_guess[:,1:5]
    @show v_new[:,1:5]

    S_guess = projection(v_guess, dist, sdist)
    coefs_guess = copy(sdist.coefficients)
    S_final = projection(v_new, dist, sdist2)
    coefs_final = copy(sdist2.coefficients)

    diff(v) = (S_final(v) .- S_guess(v)).^2

    L2_err_f, _ = hcubature(diff, [-3., -3.], [3., 3.])

    println("L2 error in f_final - f_initial: ", L2_err_f)

    max_err_v = norm(v_new .- v_guess, Inf)
    max_err_coefs = norm(coefs_final .- coefs_guess, Inf)

    println("L_inf error in v: ", max_err_v)
    println("L_inf error in coefficients: ", max_err_coefs)
    println("")

    # return solution at t
    return v_new
end
