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

    model = lenard_bernstein(M)
    v₀ = copy(model.dist.particles.v[1, :])
    method = VlasovMethods.GeometricIntegrator(model, tspan, tstep)
    h5file = joinpath(mktempdir(), "geometric_integrator.hdf5")
    run!(method, h5file)
    z = h5read(h5file, "z")
    t = h5read(h5file, "t")

    @test size(z) == (length(v₀), ntime(method.equation) + 1)
    @test t ≈ 0.0:tstep:tspan[2]
    @test z[:, 1] == v₀
    @test z[:, 2] != z[:, 1]
    # the last slice is the final integration state, written back into the model's particles
    @test z[:, end] == model.dist.particles.v[1, :]
    # each slice is the state of the same step of `integrate` on the same problem
    sol = GeometricIntegrators.integrate(
        VlasovMethods.GeometricIntegrator(lenard_bernstein(M), tspan, tstep).integrator)
    @test all(z[:, n + 1] == sol.q[n] for n in 0:ntime(method.equation))
end
