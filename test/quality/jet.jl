using JET
using PoissonSolvers
using VlasovMethods
using VlasovMethods: v_acceleration!, s_acceleration!
using Test

# Static optimisation analysis of the hot paths: every function of `src/` that a test file
# asserts with `@allocated`, at the concrete argument types those tests pass. A runtime
# dispatch on one of these paths is what the allocation tests measure only indirectly.
# Each further element type that a test outside `test/quality/` passes directly to the same
# method gets one line too; every such call passes `Float64`, so each method has one line.

if isdefined(JET, :JET_AVAILABLE) ? JET.JET_AVAILABLE : JET.JET_LOADABLE
    # test/distributions/spline_distribution.jl: projection onto a 1-D velocity space
    pd1 = ParticleDistribution(1, 1, 10)
    sd1 = SplineDistribution(1, 1, 21, 4, (-6.0, 6.0), 0, Free())
    v1 = pd1.particles.v[1, :]
    @test isempty(JET.get_reports(JET.report_opt(projection,
        (typeof(v1), typeof(pd1), typeof(sd1)); target_modules = (VlasovMethods,))))

    # test/distributions/spline_distribution.jl: projection onto a 2-D velocity space
    pd2 = ParticleDistribution(1, 2, 10)
    sd2 = SplineDistribution(1, 2, 11, 4, (-6.0, 6.0), 0, Free())
    @test isempty(JET.get_reports(JET.report_opt(projection,
        (typeof(pd2.particles.v), typeof(pd2), typeof(sd2));
        target_modules = (VlasovMethods,))))

    # test/integration/vlasov_poisson.jl: the Vlasov–Poisson right-hand side, on its two
    # `@allocated` calls, at the `Float64` argument types of that test.
    pmodel = VlasovPoisson(ParticleDistribution(1, 1, 10),
        Potential(PeriodicBasisSpline((0.0, 1.0), 3, 16)))
    pmethod = SplittingMethod(pmodel, (0.0, 1.0), 0.1)
    pparams = pmethod.equation.parameters
    pz = pmodel.distribution.particles.z
    pż = similar(pz)
    @test isempty(JET.get_reports(JET.report_opt(v_acceleration!,
        (typeof(pż), Float64, typeof(pz), typeof(pparams));
        target_modules = (VlasovMethods,))))
    @test isempty(JET.get_reports(JET.report_opt(s_acceleration!,
        (typeof(pż), Float64, typeof(pz), Float64, typeof(pparams));
        target_modules = (VlasovMethods,))))

    # test/gridbased/collisions.jl: the three hot index paths, asserted `@allocated == 0`
    nx, nv = 3, 5
    vc = collect(range(-3.0, 3.0; length = nv))
    qc = VlasovMethods.QuadraticCollisions(
        nx, nv, 1.0 / nx, (vc[end] - vc[1]) / (nv - 1), vc)
    ct = VlasovMethods.CollisionTensor(Float64, nx, nv, qc)
    Nc = nx * nv
    rt = VlasovMethods.ReducedCollisionTensor(ct, ones(Nc, 2), ones(Nc, 3), ones(Nc, 2))
    @test isempty(JET.get_reports(JET.report_opt(qc,
        (CartesianIndex{2}, CartesianIndex{2}, CartesianIndex{2});
        target_modules = (VlasovMethods,))))
    @test isempty(JET.get_reports(JET.report_opt(getindex,
        (typeof(ct), Int, Int, Int); target_modules = (VlasovMethods,))))
    @test isempty(JET.get_reports(JET.report_opt(getindex,
        (typeof(rt), Int, Int, Int); target_modules = (VlasovMethods,))))
else
    @test_skip "JET does not work on this Julia version"  # aviatesk/JET.jl#681
end
