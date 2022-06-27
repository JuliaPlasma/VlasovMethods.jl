module BumpOnTailSimulation

using LaTeXStrings
using LinearAlgebra
using ParticleMethods
using ParticleMethods.BumpOnTail
using Plots


# simulation parameters

const dt = 1e-1             # timestep
const T = 50                # final time
const nₜ = Int(div(T, dt))  # nb. of timesteps
const nₕ = 16               # nb. of elements
const p  = 3                # spline degree

const nₚ = Int(5e4)         # nb. of particles

const vmin = -10
const vmax = +10

# bump-on tail instability parameters: κ, ε, a, v₀, σ
# in order: spartial perturbation wave number and amplitude of s.p., fast particle fraction, speed, and temperature
const params = (κ = 0.3, ε = 0.03, a = 0.1, v₀= 4.5, σ = 0.5, χ = 1.0)

const L = 2π / params.κ     # domain length
const h = L/nₕ              # element width

const IP = VPIntegratorParameters(dt, nₜ, nₜ+1, nₕ, nₚ)


function run()

    # B-spline Poisson solver
    poisson = PoissonSolverPBSplines(p, IP.nₕ, L)

    # initial data
    P = BumpOnTail.draw_accept_reject(nₚ, params)
    # P = BumpOnTaildraw_importance_sampling(nₚ, params)

    # initial potential
    solve!(poisson, P.x, P.w)

    # plot initial particle distribution
    F = plot_particles(P.x, P.v, P.w, 0, L, vmin, vmax)
    savefig("bump-on-tail-init.png")

    F = plot_distribution(P.x, P.v, P.w, 0, L, vmin, vmax; nx=48, nv=32)
    savefig("bump-on-tail-init-f.png")

    # create integrate cache
    IC = VPIntegratorCache(IP)

    # integrate all time steps
    # for t = 1:nₜ
        integrate_vp!(P, poisson, params, IP, IC)
    # end

    # plot diagnostics
    plot(0:dt:dt*nₜ, IC.W, yaxis = :log, legend = :none)
    savefig("bump-on-tail-W.png")

    plot(0:dt:dt*nₜ, IC.K, yaxis = :log, legend = :none)
    savefig("bump-on-tail-K.png")
    
    plot(0:dt:dt*nₜ, IC.W .+ IC.K, yaxis = :log, legend = :none)
    savefig("bump-on-tail-E.png")
    
    # create animation and save to file
    anim = @animate for t in 0:IP.nₜ
        plot_particles(IC.X[:,t+1], IC.V[:,t+1], P.w, 0, L, vmin, vmax)
    end
    gif(anim, "bump-on-tail.gif", fps=10)
end

end # module


using .BumpOnTailSimulation

BumpOnTailSimulation.run()
