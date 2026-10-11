using GeometricEquations: ntime
using GeometricIntegrators: GeometricIntegrators
using HDF5: h5read
using StaticArrays
using Test
using VlasovMethods

include("../helpers/external_fields.jl")

@testset "run!(::SplittingMethod, h5file), $T" for T in (Float64, Float32)
    # The gyro benchmark's grid: ten periods at a step that does not divide the span.
    m, B₀ = T(3), T(1.7)
    q = T(2)
    T_g = 2 * T(π) * m / (abs(q) * B₀)
    x₀ = SVector{3, T}(0.3, -0.2, 0.1)
    v₀ = SVector{3, T}(0.7, 0.4, 0.25)
    field = UniformField(zero(SVector{3, T}), SVector{3, T}(0, 0, B₀))
    model = VlasovWithExternalFields(one_particle(T, x₀, v₀), field; charge = q, mass = m)
    method = SplittingMethod(model, (zero(T), 10 * T_g), T_g / 32)

    # The equation's `ntime` rounds the step count up (`div(t₁ - t₀, Δt, RoundUp)`) while the
    # solution carries the range `t₀:Δt:t₁`, one state per point of it. The two differ by one
    # here, which is what the run below has to take its step count from.
    @test ntime(method.equation) == 321
    sol = GeometricIntegrators.integrate(method.integrator)
    @test length(sol.t) == 321
    @test ntime(method.equation) + 1 != length(sol.t)

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
