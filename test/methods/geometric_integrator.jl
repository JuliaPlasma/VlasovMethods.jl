using GeometricIntegrators: GeometricIntegrators
using GeometricEquations: ntime
using HDF5: h5read
using Random
using SimpleSplines: BSplineBasis, UniformMesh, Free, SplineQuadrature, nbasis, (..)
using Test
using VlasovMethods

# A uniform particle distribution inside the spline domain, so the projection stays positive and
# the conservative models do not stop on the positivity check (KNOWN_ISSUES.md, K8).
function lenard_bernstein(M; npart = 200)
    Random.seed!(1234)
    dist = initialize!(ParticleDistribution(1, 1, npart),
        UniformDistribution((0.0, 1.0), (-2.0, 2.0)))
    axis = BSplineBasis(UniformMesh(7, -2.5 .. 2.5), 2, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(nbasis(axis)))
    M(dist, CollisionEntropy(sdist))
end

@testset "run! of a GeometricIntegrator, $M" for M in (LenardBernstein,
    ConservativeLenardBernstein, RescaledConservativeLenardBernstein)
    tspan = (0.0, 5e-3)
    tstep = 1e-3
    times = tspan[1]:tstep:tspan[2]

    model = lenard_bernstein(M)
    v₀ = copy(model.dist.particles.v[1, :])
    method = VlasovMethods.GeometricIntegrator(model, tspan, tstep)
    h5file = joinpath(mktempdir(), "geometric_integrator.hdf5")
    run!(method, h5file)
    z = h5read(h5file, "z")
    t = h5read(h5file, "t")

    @test size(z) == (length(v₀), length(times))
    @test t ≈ times
    @test z[:, 1] == v₀
    @test z[:, 2] != z[:, 1]
    # the last slice is the final integration state, written back into the model's particles
    @test z[:, end] == model.dist.particles.v[1, :]
    # each slice is the state of the same step of `integrate` on the same problem
    sol = GeometricIntegrators.integrate(
        VlasovMethods.GeometricIntegrator(lenard_bernstein(M), tspan, tstep).integrator)
    @test ntime(sol) + 1 == length(times)
    @test all(z[:, n + 1] == sol.q[n] for n in 0:ntime(sol))
end

@testset "the method argument of GeometricIntegrator(::LenardBernstein, …)" begin
    tspan = (0.0, 5e-3)
    tstep = 1e-3
    build(method...) = VlasovMethods.GeometricIntegrator(
        lenard_bernstein(LenardBernstein), tspan, tstep, method...)

    default = build()
    midpoint = build(GeometricIntegrators.ImplicitMidpoint())
    @test GeometricIntegrators.Integrators.method(midpoint.integrator) isa
          GeometricIntegrators.ImplicitMidpoint

    # the default is `RK438`: its run is that of `RK438()` passed as the method, step for step
    rk = GeometricIntegrators.integrate(default.integrator)
    rk438 = GeometricIntegrators.integrate(build(GeometricIntegrators.RK438()).integrator)
    @test all(rk.q[n] == rk438.q[n] for n in 0:ntime(rk))

    # the first step of the given method solves the implicit midpoint rule, to a relative
    # residual of 4e-14, and the default's does not, at 6e-5
    sol = GeometricIntegrators.integrate(midpoint.integrator)
    problem = midpoint.equation
    residual(q) =
        let v̇ = similar(q[0])
            problem.equation.v(v̇, 0.0, (q[0] .+ q[1]) ./ 2, problem.parameters)
            maximum(abs, q[1] .- q[0] .- tstep .* v̇) / maximum(abs, q[1] .- q[0])
        end
    @test residual(sol.q) < 1e-10
    @test residual(rk.q) > 1e-6
end

# `ntime` rounds the step count of `0.05 / 8e-4` up to 63, and the range `0.0:8e-4:0.05` has 63
# points, so 62 steps; a loop over `1:ntime` reads past its end.
@testset "run! of a GeometricIntegrator over a span that is not a whole number of steps" begin
    tspan = (0.0, 0.05)
    tstep = 8e-4
    times = tspan[1]:tstep:tspan[2]

    model = lenard_bernstein(LenardBernstein)
    method = VlasovMethods.GeometricIntegrator(model, tspan, tstep)
    @test ntime(method.equation) == length(times) == 63
    h5file = joinpath(mktempdir(), "geometric_integrator.hdf5")
    @test run!(method, h5file) === model.dist
    z = h5read(h5file, "z")
    t = h5read(h5file, "t")

    @test size(z, 2) == length(t) == length(times)
    @test t == collect(times)
    @test z[:, end] == model.dist.particles.v[1, :]
    # the slices are the states of `integrate` on the same problem, which steps over `times`
    sol = GeometricIntegrators.integrate(
        VlasovMethods.GeometricIntegrator(lenard_bernstein(LenardBernstein), tspan, tstep).integrator)
    @test ntime(sol) + 1 == length(times)
    @test all(z[:, n + 1] == sol.q[n] for n in 0:ntime(sol))
end
