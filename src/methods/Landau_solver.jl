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
        L = VlasovMethods.compute_J(sdist)

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

function f!(F, v_nplus1 , v_n, rhs, params, Δt)
    v_midpoint = (v_nplus1 .+ v_n) / 2

    # rhs = copy(v_n)
    Landau_rhs_2!(rhs, v_midpoint, params )

    explicit_update!(view(F,1,:), v_n[1,:], Δt, rhs[1,:])
    explicit_update!(view(F,2,:), v_n[2,:], Δt, rhs[2,:])

    F .-= v_nplus1

end


function Picard_iterate_Landau_nls!(dist, sdist, tol, ftol, β, Δt, ti, t, v_prev, v_prev_2, rhs_prev, sdist2, max_iters, m)
    # β is the damping parameter for damped Picard iterations, with β = 1 yielding regular Picard iterations
    # ti is the time index at which v_new is being computed, i.e. for t = ti * Δt
    # v_prev is v at the previous timestep 
    # v_prev_2 is v at two timesteps prior to t
    # rhs_prev[:,:,1] is the rhs at t - Δt, and rhs_prev[:,:,2] is the rhs at t - 2Δt
    
    # creating this to store the guess for the moment, for diagnostic purposes
    v_guess = copy(dist.particles.v) 

    params = (dist = dist, sdist = sdist)
    
    # use Hermite extrapolation to get an initial guess
    if ti ≥ 4
        Extrapolators.extrapolate!(t - 2Δt, v_prev_2, rhs_prev[:,:,2], t - Δt, v_prev, rhs_prev[:,:,1], t, v_guess, Extrapolators.HermiteExtrapolation())
    end

    rhs_prev[:,:,2] .= rhs_prev[:,:,1]

    g!(F, v) = f!(F, v, v_prev, rhs_prev[:,:,1], params, Δt) 

    println("nlsolve")
    @time sol = nlsolve(g!, v_guess, method=:anderson, iterations = max_iters, m = m, beta = β, xtol = tol, ftol = ftol, show_trace = true)

    # v_new = sol.zero

    # update solution array
    dist.particles.v .= sol.zero

    # update rhs storage
    rhs_prev[:,:,1] .= Landau_rhs_2!(rhs_prev[:,:,1], dist.particles.v, params)

    # compute some diagnostics

    # @show v_guess[:,1:5]
    # @show dist.particles.v[:,1:5]

    # S_guess = projection(v_guess, dist, sdist)
    # coefs_guess = copy(sdist.coefficients)
    # S_final = projection(dist.particles.v, dist, sdist2)
    # coefs_final = copy(sdist2.coefficients)

    # diff(v) = (S_final(v) .- S_guess(v)).^2

    # L2_err_f, _ = hcubature(diff, [-3., -3.], [3., 3.])

    # println("L2 error in f_final - f_initial: ", L2_err_f)

    # max_err_v = norm(dist.particles.v .- v_guess, Inf)
    # max_err_coefs = norm(coefs_final .- coefs_guess, Inf)

    # println("L_inf error in v: ", max_err_v)
    # println("L_inf error in coefficients: ", max_err_coefs)
    # println("")

    # return solution at t
    return dist.particles.v
end


function Picard_iterate_Landau!(dist, sdist, tol, β, Δt, ti, t, v_prev_2, rhs_prev, sdist2, max_iters, m)
    # β is the damping parameter for damped Picard iterations, with β = 1 yielding regular Picard iterations
    # ti is the time index at which v_new is being computed, i.e. for t = ti * Δt
    # dist.particles.v is v at the previous timestep 
    # v_prev_2 is v at two timesteps prior to t
    # rhs_prev[:,:,1] is the rhs at t - Δt, and rhs_prev[:,:,2] is the rhs at t - 2Δt
    
    # create vectors to store the current and previous iterations
    v_new = zero(dist.particles.v)
    v_prev = copy(dist.particles.v) #TODO: this should be updated to computing an initial guess using Hermite extrapolation
    
    # creating this to store the guess for the moment, for diagnostic purposes
    v_guess = copy(dist.particles.v) 

    params = (dist = dist, sdist = sdist)
    
    # use Hermite extrapolation to get an initial guess
    if ti ≥ 4
        Extrapolators.extrapolate!(t - 2Δt, v_prev_2, rhs_prev[:,:,2], t - Δt, dist.particles.v, rhs_prev[:,:,1], t, v_guess, Extrapolators.HermiteExtrapolation())
    end

    rhs_prev[:,:,2] .= rhs_prev[:,:,1]

    v_prev .= v_guess

    # j is the Picard iteration index
    j = 0
    # println("entering Picard loop")
    err = 1.
    # start with no damping
    beta = 1.
    @time while err > tol
        println("j=",j)

        # compute midpoint either using initial guess or the previous Picard iteration
        v_midpoint = (v_prev .+ dist.particles.v) / 2

        # compute rhs using the approximated v_midpoint
        rhs_prev[:,:,1] .= Landau_rhs_2!(rhs_prev[:,:,1], v_midpoint, params )
        
        # perform explicit time update
        explicit_update!(view(v_new,1,:), dist.particles.v[1,:], Δt, rhs_prev[1,:,1])
        explicit_update!(view(v_new,2,:), dist.particles.v[2,:], Δt, rhs_prev[2,:,1])

        # apply damping
        v_new .*= beta 
        v_new .+= (1 - beta)*v_prev

        # update rhs storage
        rhs_prev[:,:,1] .= Landau_rhs_2!(rhs_prev[:,:,1], v_new, params )

        # compute error measures for iteration control
        err = norm(v_new .- v_prev)
        err_max = norm(v_new .- v_prev, Inf)

        println("L2 norm of (v_j+1 - v_j) = ", err)
        println("Linf norm of (v_j+1 - v_j) = ", err_max)
        println("")

        

        if j - 1 == max_iters && beta == 1.
            @warn " Maximum iterations reached, aborting regular Picard iterations and restarting with supplied damping parameter"
            v_prev .= v_guess
            beta = β
            j = 0
        end

    end

    # update solution array
    dist.particles.v .= v_new

    # compute some diagnostics

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
