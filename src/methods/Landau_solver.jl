# function IM_rule!(g, x, y, f, Δt)
#     g .= x .- y .- Δt * f(0.5 * (y .+ x))
# end


function IM_rule!(g, x, y, f, Δt, params)
    # g is the function to find the zero of 
    # x is the root 
    # y is the previous timestep for IM 
    # f is the vector field 
    # 1/2 * (x + y) is the midpoint 
    g .= x .- y .- Δt .* f((y .+ x) ./ 2, params)
end

# function IM_update!(y_new, y_new_guess, yₙ, f, Δt, params, ::Landau)
#     # println("define F!")
#     F!(g, x) = IM_rule!(g, x, yₙ, f, Δt, params)
#     obj = zero(yₙ)
#     opt = Options(x_reltol = 1e-8, f_reltol = 1e-8, x_abstol = 1e-12, f_abstol = 1e-12)
#     # println("construct NewtonSolver")
#     nl = NewtonSolver(y_new_guess, obj, config = opt)
#     # println("Newton solve")
#     SimpleSolvers.solve!(y_new, F!, nl)

#     return y_new
# end


# function explicit_update!(rhs, v, Δt)
#     @. rhs = v + Δt * rhs
# end

# function Picard_iterate_over_particles(f::Function, dist, sdist, tol, Δt)
#     # j is the Picard iteration index
#     # create vectors to store the current iteration and the two previous iterations
#     v_new = zero(dist.particles.v)
#     v_prev = copy(dist.particles.v) #TODO: this should be updated to computing an initial guess using Hermite extrapolation
    
#     err = 1.
#     j = 0
#     println("entering Picard loop")
#     while err > tol
#         println("j=", j)

#         S = projection(v_prev, dist, sdist)
#         L = VlasovMethods.compute_J(sdist)

#         params = (dist = dist, sdist = sdist, B = sdist.basis, L = L, v_array = v_prev)

#         @time for i in axes(v_prev,2)
#         # @time for i in axes(v_prev,2)
#             # @show i
#             IM_update!(view(v_new, :, i), view(v_prev, :, i), view(dist.particles.v, :, i), f, Δt, params)
#         end

#         @show err = sqeuclidean(v_new, v_prev)
#         v_prev .= v_new
#         j += 1
#     end

#     dist.particles.v .= v_new
    
#     return v_new
# end

function f!(f::AbstractArray{T}, vn::AbstractArray{T}, vp, params, Δt, landau) where T
    v_midpoint  = landau.cache[T].v
    v_midpoint .= (vn .+ vp) ./ 2

    collisional_vectorfield!(f, v_midpoint, params, landau)

    f .*= Δt
    f .-= (vn .- vp)

    return nothing
end

# function f(dv, v_nplus1::AbstractArray{T}, v_n, params, Δt) where T

#     # F = copy(v_nplus1)
#     rhs = copy(v_nplus1)

#     v_midpoint = (v_nplus1 + v_n) / 2

#     # rhs = copy(v_n)
#     Landau_rhs_2!(rhs, v_midpoint, params)

#     explicit_update!(view(rhs, 1, :), v_n[1,:], Δt)
#     explicit_update!(view(rhs, 2, :), v_n[2,:], Δt)

#     rhs .-= v_nplus1

#     return rhs::AbstractArray{T}

# end


function Picard_iterate_Landau_nls!(landau, tol, ftol, β, Δt, ti, t, v_prev, v_prev_2, rhs_prev, m, n, chunksize)
    # β is the damping parameter for damped Picard iterations, with β = 1 yielding regular Picard iterations
    # ti is the time index at which v_new is being computed, i.e. for t = ti * Δt
    # v_prev is v at the previous timestep 
    # v_prev_2 is v at two timesteps prior to t
    # rhs_prev[:,:,1] is the rhs at t - Δt, and rhs_prev[:,:,2] is the rhs at t - 2Δt

    dist = landau.dist
    ent = landau.entropy
    
    # creating this to store the guess for the moment, for diagnostic purposes
    v_guess = copy(dist.particles.v)
    
    params = (dist = dist, ent = ent, n = n)
    
    # use Hermite extrapolation to get an initial guess
    if ti ≥ 4
        Extrapolators.extrapolate!(t - 2Δt, v_prev_2, view(rhs_prev, :, :, 2), t - Δt, v_prev, view(rhs_prev, :, :, 1), t, v_guess, Extrapolators.HermiteExtrapolation())
    else
        problemGNI = GeometricEquations.ODEProblem((v̇,t,v,params) -> collisional_vectorfield!(v̇,v,params,landau), (t, t+Δt), Δt, v_prev; parameters = params)
        Extrapolators.extrapolate!(t - Δt, v_prev, t, v_guess, problemGNI, Extrapolators.MidpointExtrapolation(5))
    end

    probN = NonlinearProblem{true}((f, v, p) -> f!(f, v, v_prev, params, Δt, landau), v_guess)

    # println("nlsolve")
    # NonlinearSolve.jl using NewtonRaphson
    # @time sol = NonlinearSolve.solve(probN, 
    #     NewtonRaphson(linsolve = AppleAccelerateLUFactorization(), autodiff = AutoForwardDiff(; chunksize = chunksize)); 
    #     reltol = 5e-3, show_trace=Val(true), trace_level = TraceWithJacobianConditionNumber())
    
    # NonlinearSolve.jl using Picard w/ anderson acceleration
    @time sol = NonlinearSolve.solve(probN, 
        NonlinearSolve.NLsolveJL(; method=:anderson, m = m, beta = β); 
        reltol = 5e-3, show_trace=Val(true))

    # @time sol = nlsolve(g!, v_guess, method=:anderson, iterations = max_iters, m = m, beta = β, xtol = tol, ftol = ftol, show_trace = true)

    # v_new = sol.u

    # update solution array
    dist.particles.v .= sol.u

    # update rhs storage
    # rhs_prev[:,:,1] .= Landau_rhs_2!(view(rhs_prev, :, :, 1), dist.particles.v, params)
    rhs_prev[:,:,2] .= rhs_prev[:,:,1]
    collisions_rhs!(view(rhs_prev, :, :, 1), dist.particles.v, params, landau)

    # return solution at t
    return sol
end


# function Picard_iterate_Landau!(dist, sdist, tol, β, Δt, ti, t, v_prev_2, rhs_prev, sdist2, max_iters, m)
#     # β is the damping parameter for damped Picard iterations, with β = 1 yielding regular Picard iterations
#     # ti is the time index at which v_new is being computed, i.e. for t = ti * Δt
#     # dist.particles.v is v at the previous timestep 
#     # v_prev_2 is v at two timesteps prior to t
#     # rhs_prev[:,:,1] is the rhs at t - Δt, and rhs_prev[:,:,2] is the rhs at t - 2Δt
    
#     # create vectors to store the current and previous iterations
#     v_new = zero(dist.particles.v)
#     v_prev = copy(dist.particles.v) #TODO: this should be updated to computing an initial guess using Hermite extrapolation
    
#     # creating this to store the guess for the moment, for diagnostic purposes
#     v_guess = copy(dist.particles.v) 

#     params = (dist = dist, sdist = sdist)
    
#     # use Hermite extrapolation to get an initial guess
#     if ti ≥ 4
#         Extrapolators.extrapolate!(t - 2Δt, v_prev_2, rhs_prev[:,:,2], t - Δt, dist.particles.v, rhs_prev[:,:,1], t, v_guess, Extrapolators.HermiteExtrapolation())
#     end

#     rhs_prev[:,:,2] .= rhs_prev[:,:,1]

#     v_prev .= v_guess

#     # j is the Picard iteration index
#     j = 0
#     # println("entering Picard loop")
#     err = 1.
#     # start with no damping
#     beta = 1.
#     @time while err > tol
#         println("j=",j)

#         # compute midpoint either using initial guess or the previous Picard iteration
#         v_midpoint = (v_prev .+ dist.particles.v) / 2

#         # compute rhs using the approximated v_midpoint
#         rhs_prev[:,:,1] .= Landau_rhs_2!(rhs_prev[:,:,1], v_midpoint, params )
        
#         # perform explicit time update
#         explicit_update!(view(v_new,1,:), dist.particles.v[1,:], Δt, rhs_prev[1,:,1])
#         explicit_update!(view(v_new,2,:), dist.particles.v[2,:], Δt, rhs_prev[2,:,1])

#         # apply damping
#         v_new .*= beta 
#         v_new .+= (1 - beta)*v_prev

#         # update rhs storage
#         rhs_prev[:,:,1] .= Landau_rhs_2!(rhs_prev[:,:,1], v_new, params )

#         # compute error measures for iteration control
#         err = norm(v_new .- v_prev)
#         err_max = norm(v_new .- v_prev, Inf)

#         println("L2 norm of (v_j+1 - v_j) = ", err)
#         println("Linf norm of (v_j+1 - v_j) = ", err_max)
#         println("")

        

#         if j - 1 == max_iters && beta == 1.
#             @warn " Maximum iterations reached, aborting regular Picard iterations and restarting with supplied damping parameter"
#             v_prev .= v_guess
#             beta = β
#             j = 0
#         end

#     end

#     # update solution array
#     dist.particles.v .= v_new

#     # compute some diagnostics

#     @show v_guess[:,1:5]
#     @show v_new[:,1:5]

#     S_guess = projection(v_guess, dist, sdist)
#     coefs_guess = copy(sdist.coefficients)
#     S_final = projection(v_new, dist, sdist2)
#     coefs_final = copy(sdist2.coefficients)

#     diff(v) = (S_final(v) .- S_guess(v)).^2

#     L2_err_f, _ = hcubature(diff, [-3., -3.], [3., 3.])

#     println("L2 error in f_final - f_initial: ", L2_err_f)

#     max_err_v = norm(v_new .- v_guess, Inf)
#     max_err_coefs = norm(coefs_final .- coefs_guess, Inf)

#     println("L_inf error in v: ", max_err_v)
#     println("L_inf error in coefficients: ", max_err_coefs)
#     println("")

#     # return solution at t
#     return v_new
# end
