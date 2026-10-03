# Measure the allocation of one evaluation of the Vlasov–Poisson right-hand side.
#
# Run it in a fresh process, without code-coverage instrumentation:
#
#     julia --startup-file=no --project=scripts scripts/vlasov_poisson_allocations.jl
#
# Each function is called once to compile it, and the second call is measured. The script exits
# with status 1 if any evaluation allocates.

using PoissonSolvers
using Random
using VlasovMethods
using VlasovMethods: s_acceleration!, v_acceleration!

function measure(::Type{T}; npart = 10000, ncells = 16, order = 3) where {T}
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

    s_acceleration!(ż, t, z̄, t̄, params)
    v_acceleration!(ż, t, z̄, params)

    return (s_acceleration! = @allocated(s_acceleration!(ż, t, z̄, t̄, params)),
        v_acceleration! = @allocated(v_acceleration!(ż, t, z̄, params)))
end

results = [T => measure(T) for T in (Float64, Float32)]
for (T, bytes) in results
    println(rpad(string(T), 8), " s_acceleration! ", bytes.s_acceleration!, " bytes, ",
        "v_acceleration! ", bytes.v_acceleration!, " bytes")
end

exit(all(r -> all(iszero, values(r[2])), results) ? 0 : 1)
