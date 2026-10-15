# Measure the allocation of one evaluation of the Vlasov–Poisson right-hand side.
#
#   julia --project=scripts scripts/vlasov_poisson_allocations.jl \
#       [--npart 10000] [--ncells 16] [--order 3] [--runs-dir dir] [--results-dir dir]
#
# Run it in a fresh process, without code-coverage instrumentation. Each function is called once
# to compile it, and the second call is measured. The numbers go to `--runs-dir`, the report with
# the verdict to `--results-dir`, and the script exits with status 1 if any evaluation allocates.

include(joinpath(@__DIR__, "options.jl"))

using Printf

using PoissonSolvers
using Random
using VlasovMethods
using VlasovMethods: SAcceleration, VAcceleration

const OPT = parse_arguments(ARGS; npart = 10000, ncells = 16, order = 3)
const NPART = OPT.npart
const NCELLS = OPT.ncells
const ORDER = OPT.order

function measure(::Type{T}; npart = NPART, ncells = NCELLS, order = ORDER) where {T}
    # `ParticleDistribution(1, 1, npart)` always fills `Float64` zeros, so the particles of
    # element type `T` are built here.
    z = zeros(T, 3, npart)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    dist = ParticleDistribution(1, 1, VlasovMethods.ParticleList(z; variables = vars))
    Random.seed!(1234)
    initialize!(dist, NormalDistribution((T(0), T(1))))

    potential = Potential(PeriodicBasisSpline((T(0), T(1)), order, ncells))
    method = SplittingMethod(VlasovPoisson(dist, potential), (T(0), T(1)), T(0.1))
    params = method.equation.parameters

    z̄ = copy(dist.particles.z)
    ż = similar(z̄)
    t, t̄ = T(0.1), zero(T)

    s_acc = SAcceleration(method.model)
    v_acc = VAcceleration(method.model)
    s_acc(ż, t, z̄, t̄, params)
    v_acc(ż, t, z̄, params)

    return (SAcceleration = @allocated(s_acc(ż, t, z̄, t̄, params)),
        VAcceleration = @allocated(v_acc(ż, t, z̄, params)))
end

results = [T => measure(T) for T in (Float64, Float32)]

lines = [@sprintf("%-8s SAcceleration %d bytes, VAcceleration %d bytes",
             T, bytes.SAcceleration, bytes.VAcceleration) for (T, bytes) in results]

numbers = runs_path("vlasov_poisson_allocations.txt")
report = results_path("vlasov_poisson_allocations.md")
clean = all(r -> all(iszero, values(r[2])), results)

open(numbers, "w") do io
    for line in lines
        println(io, line)
    end
end

open(report, "w") do io
    println(io, "# Vlasov–Poisson right-hand-side allocations")
    println(io)
    @printf(io,
        "`--npart %d --ncells %d --order %d`, one warmed call of each acceleration.\n",
        NPART, NCELLS, ORDER)
    println(io)
    for line in lines
        println(io, "    ", line)
    end
    println(io)
    println(io,
        clean ?
        "Both precisions allocate zero bytes on both calls, as the scheme requires." :
        "**At least one call allocates.** The right-hand side is on the hot path of " *
        "every step and must allocate nothing.")
end

for line in lines
    println(line)
end
@printf("  wrote %s and %s\n", numbers, report)

exit(clean ? 0 : 1)
