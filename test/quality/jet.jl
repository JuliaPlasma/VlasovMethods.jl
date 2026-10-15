using JET
using LinearAlgebra
using PoissonSolvers
using Random
using ReducedBasisMethods: get_DEIM_indices
using VlasovMethods
using VlasovMethods: DEIMElectricField, ReducedElectricField, SAcceleration,
                     VAcceleration, VlasovHamiltonian, efield!, energy, s_drift!,
                     s_lorentz!,
                     update!
using ElectromagneticFields
using Test

include("../helpers/external_fields.jl")
include("../helpers/poisson_field.jl")

# Static optimisation analysis of the hot paths: every function of `src/` that a test file
# asserts with `@allocated`, at the concrete argument types those tests pass. A runtime
# dispatch on one of these paths is what the allocation tests measure only indirectly.
# Each further element type that a test outside `test/quality/` passes directly to the same
# method is analysed at it too: the external-field tests run in `Float64` and `Float32`, so
# each of their methods is analysed once per precision.

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
    @test isempty(JET.get_reports(JET.report_opt(VAcceleration(pmodel),
        (typeof(pż), Float64, typeof(pz), typeof(pparams));
        target_modules = (VlasovMethods,))))
    @test isempty(JET.get_reports(JET.report_opt(SAcceleration(pmodel),
        (typeof(pż), Float64, typeof(pz), Float64, typeof(pparams));
        target_modules = (VlasovMethods,))))

    # test/models/vlasov_with_external_fields.jl and
    # test/integration/electromagnetic_fields_ext.jl: the external-field push, on the
    # `@allocated` calls of `test_push_does_not_allocate`, at the argument types those files pass
    # — a copied state matrix, in both precisions, with the constant `UniformField` and the ABC
    # `FieldFunctions` the extension adapts.
    for T in (Float64, Float32)
        for field in (UniformField(zeros(T, 3), zeros(T, 3)),
            FieldFunctions(ABCEquilibrium(T(1.3), T(0.7), T(0.45))))
            model = VlasovWithExternalFields(
                particle_distribution(T, 10), field; charge = 2, mass = 3)
            params = SplittingMethod(model, (zero(T), one(T)), T(0.1)).equation.parameters
            z = copy(model.distribution.particles.z)
            for f in (s_drift!, s_lorentz!)
                @test isempty(JET.get_reports(JET.report_opt(f,
                    (typeof(z), T, typeof(z), T, typeof(params));
                    target_modules = (VlasovMethods,))))
            end
        end
    end
    # test/integration/sodeensemble.jl: the Hamiltonian invariant, on its `@allocated` call, at
    # the types of an ensemble member's state and parameters, which the test builds as copies.
    pensemble = (w = copy(pmodel.distribution.particles.w),)
    @test isempty(JET.get_reports(JET.report_opt(VlasovHamiltonian(pmodel),
        (Float64, typeof(copy(pz)), typeof(pensemble));
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

    # test/electric_field.jl and test/particles/electric_field.jl: the full-order field and the
    # two reduced fields, on their `@allocated` calls, at the `Float64` and `Float32` argument
    # types those files pass. Each field is analysed once per precision, because each test file
    # runs both.
    for T in (Float64, Float32)
        # `poisson_field` sets the same seed; the `randn` draws below continue its stream
        Random.seed!(1234)
        efield, q, w = poisson_field(T; n = 200)
        Ψₓ = Matrix(qr(randn(T, 200, 20)).Q)
        Ψₑ = Matrix(qr(randn(T, 200, 8)).Q)
        p = get_DEIM_indices(Ψₑ)
        X = Ψₓ' * q
        E = similar(X)
        for f in (efield, ReducedElectricField(efield, Ψₓ, X),
            DEIMElectricField(efield, Ψₓ, Ψₑ, p, X))
            @test isempty(JET.get_reports(JET.report_opt(update!,
                (typeof(f), typeof(X), typeof(w), T); target_modules = (VlasovMethods,))))
            @test isempty(JET.get_reports(JET.report_opt(efield!,
                (typeof(f), typeof(E), typeof(X)); target_modules = (VlasovMethods,))))
            @test isempty(JET.get_reports(JET.report_opt(f,
                (typeof(E), typeof(X), typeof(w), T); target_modules = (VlasovMethods,))))
            @test isempty(JET.get_reports(JET.report_opt(energy,
                (typeof(f),); target_modules = (VlasovMethods,))))
        end
    end
else
    @test_skip "JET does not work on this Julia version"  # aviatesk/JET.jl#681
end
