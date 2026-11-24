using GeometricEquations: ntime
using GeometricIntegrators: GeometricIntegrators
using HDF5: h5read
using StaticArrays
using Test
using VlasovMethods

include("../helpers/external_fields.jl")

@testset "run!(::SplittingMethod, h5file), $T" for T in (Float64, Float32)
    # The gyro benchmark's grid: ten periods at a step of `T_g/32`, 320 steps in exact arithmetic,
    # ending one ulp below `10 * T_g`.
    m, B₀ = T(3), T(1.7)
    q = T(2)
    T_g = 2 * T(π) * m / (abs(q) * B₀)
    x₀ = SVector{3, T}(0.3, -0.2, 0.1)
    v₀ = SVector{3, T}(0.7, 0.4, 0.25)
    field = UniformField(zero(SVector{3, T}), SVector{3, T}(0, 0, B₀))
    model = VlasovWithExternalFields(one_particle(T, x₀, v₀), field; charge = q, mass = m)
    method = SplittingMethod(model, (zero(T), prevfloat(10 * T_g)), T_g / 32)

    # `T_g/32` is exact, and `10 * T_g` is half an ulp or less off `320Δt`. The span ends one ulp
    # below that, so the last step of `320Δt` ends past `t₁`: the equation's `ntime` counts it,
    # 320 steps, and the solution's grid `t₀:Δt:t₁` drops it, 319 steps on 320 points, one state
    # per point. `run!` writes one slice per point of the grid, not `ntime + 1` (K93).
    @test ntime(method.equation) == 320
    sol = GeometricIntegrators.integrate(method.integrator)
    @test length(sol.t) == 320
    @test ntime(method.equation) == length(sol.t)

    h5file = joinpath(mktempdir(), "splitting.hdf5")
    @test run!(method, h5file) === model.distribution
    h5z = h5read(h5file, "z")

    # one slice per time point of the solution, the initial state at the first of them
    @test size(h5z, 3) == length(sol.t)
    @test h5z[:, :, 1] == sol.q[0]
    @test all(h5z[:, :, n + 1] == sol.q[n] for n in 1:(length(sol.t) - 1))
    # the last slice is the state the run leaves in the model's particles
    @test h5z[:, :, end] == sol.q[end]
    @test h5z[:, :, end] == model.distribution.particles.z
end
