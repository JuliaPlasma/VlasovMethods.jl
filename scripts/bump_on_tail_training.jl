# Bump-on-tail, step 1 of 3: the training sweep of the full-order model.
#
# Draws the particles once at the reference wave number, marches them with the full-order leapfrog
# for each of the ten training wave numbers κ in [0.1, 0.5], and evaluates the field of each saved
# state. The defaults are the paper's parameters: dt = 0.1 to T = 25, 10⁴ particles, a potential
# of 16 cubic spline cells, and the seed 1234. The paper command is the default:
#
#     julia --project=scripts scripts/bump_on_tail_training.jl
#
# Writes `runs_path("<run id>.h5")`: the snapshots of the sweep, with their accelerations, and the
# integrator parameters in the group `integrator`. `bump_on_tail_projections.jl` reads it, with the
# same options.

include("options.jl")
include("bump_on_tail_common.jl")

const OPT = parse_arguments(ARGS; BUMP_ON_TAIL...)

function main(o)
    particles = initial_particles(o)
    poisson = potential(o)
    pspace = parameter_space(training_kappas())

    nt = step_count(o.final_time, o.dt)
    nh = length(PoissonSolvers.coefficients(poisson))
    IP = IntegratorParameters(o.dt, nt, nt + 1, nh, o.particles, length(pspace))
    SS = Snapshots(Float64, 1, o.particles, nh, nt + 1, pspace)

    for p in 1:length(pspace)
        params = pspace(p)
        println("training sample $p of $(length(pspace)), χ = $(params.χ)")
        integrate_vp(particles, params, vp_field(poisson, params), SS, IP, p)
    end
    accelerations!(SS, poisson, vec(particles.w))

    file = runs_path(run_id(o) * ".h5")
    h5open(file, "w") do io
        h5save(io, SS)
        h5save(io, IP; path = "integrator")
    end
    println("wrote ", file)
    return nothing
end

main(OPT)
