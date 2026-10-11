using ElectromagneticFields
using GeometricIntegrators: GeometricIntegrators
using LinearAlgebra: ×, ⋅
using StaticArrays
using Test
using TOML
using VlasovMethods
using VlasovMethods: s_drift!, s_lorentz!, v_drift!, v_lorentz!

include("../helpers/external_fields.jl")

const EXTENSION = :VlasovMethodsElectromagneticFieldsExt
const EF_UUID = "d6c1ba6f-ee03-53af-b876-68cefeb88ec8"
const ELECTROMAGNETIC_FIELDS = Base.PkgId(Base.UUID(EF_UUID), "ElectromagneticFields")

# What a fresh process asserts. It is its own program: the parent compares the exit status and no
# output, so the `--color=yes` of CI cannot split a pattern.
const CHILD = """
    using VlasovMethods
    @assert Base.get_extension(VlasovMethods, :$EXTENSION) === nothing
    @assert !haskey(Base.loaded_modules, Base.PkgId(Base.UUID("$EF_UUID"), "ElectromagneticFields"))
    using ElectromagneticFields
    @assert Base.get_extension(VlasovMethods, :$EXTENSION) isa Module
    """

# The magnetic field of `ABCEquilibrium(a, b, c)`, which is equal to its own vector potential:
# `(a sin z + c cos y, b sin x + a cos z, c sin y + b cos x)`. Written out here rather than read
# from the package, so that the assertion is against the formula and not against the code that
# produced it.
function abc_field(::Type{T}, a, b, c, x) where {T}
    SVector{3, T}(
        a * sin(x[3]) + c * cos(x[2]),
        b * sin(x[1]) + a * cos(x[3]),
        c * sin(x[2]) + b * cos(x[1]))
end

# Five points with distinct non-power-of-two coordinates, so that a permuted component of the
# adapter is caught.
const ABC_POINTS = ((0.31, -0.52, 0.77), (-1.13, 0.63, 2.05), (2.71, -1.94, 0.37),
    (0.13, 0.29, -2.17), (-0.71, 2.31, 1.41))

@testset "the ElectromagneticFields extension" begin
    @testset "the extension loads only on demand" begin
        # `ElectromagneticFields` is a dependency of this test environment, so the extension is
        # loaded in this process; the clause is about a process that does not load the package.
        @test Base.get_extension(VlasovMethods, EXTENSION) isa Module
        @test ELECTROMAGNETIC_FIELDS in keys(Base.loaded_modules)

        # A fresh process, in the environment `Pkg.test` built.
        @test success(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $CHILD`)

        # `[deps]` and not `[weakdeps]` is what would make every user's environment carry the
        # package; the extension is declared on the weak one.
        project = TOML.parsefile(joinpath(pkgdir(VlasovMethods), "Project.toml"))
        @test haskey(project["weakdeps"], "ElectromagneticFields")
        @test !haskey(get(project, "deps", Dict{String, Any}()), "ElectromagneticFields")
        @test project["extensions"] ==
              Dict("VlasovMethodsElectromagneticFieldsExt" => "ElectromagneticFields")
    end

    @testset "the ABC field, $T" for T in (Float64, Float32)
        a, b, c = T(1.3), T(0.7), T(0.45)
        field = FieldFunctions(ABCEquilibrium(a, b, c))
        tolerance = 4 * eps(T) * (a + b + c)

        @testset "the adapter returns the field" begin
            for p in ABC_POINTS
                x = SVector{3, T}(p)
                @test magnetic_field(field, zero(T), x) ≈
                      abc_field(T, a, b, c, x) atol = tolerance
                # `ABCEquilibrium` defines no `φ`, so the default `φ = 0` applies and the
                # electric field is zero
                @test all(iszero, electric_field(field, zero(T), x))
            end
        end

        @testset "the perturbation" begin
            E₀ = T(0.3)
            perturbed = FieldFunctions(ABCEquilibrium(a, b, c), EzCosZPerturbation(E₀))
            for p in ABC_POINTS
                x = SVector{3, T}(p)
                E = electric_field(perturbed, zero(T), x)
                @test E ≈
                      SVector{3, T}(0, 0, -E₀ * cos(2 * T(π) * x[3])) atol = 4 * eps(T) * E₀
                # the perturbation leaves the magnetic field as the equilibrium alone had it
                @test magnetic_field(perturbed, zero(T), x) ≈
                      magnetic_field(field, zero(T), x)
            end
        end

        @testset "an orbit in the ABC field" begin
            q, m = T(2), T(3)
            x₀ = SVector{3, T}(0.3, -0.2, 0.1)
            v₀ = SVector{3, T}(0.7, 0.4, 0.25)
            tstep = T(5) / T(100)
            model = VlasovWithExternalFields(
                one_particle(T, x₀, v₀), field; charge = q, mass = m)
            method = SplittingMethod(model, (zero(T), T(5)), tstep)
            sol = GeometricIntegrators.integrate(method.integrator)
            # The solution's own step count: the equation's `ntime` rounds up, the time series is
            # the range, and the two disagree wherever the span is not an exact multiple.
            N = lastindex(sol.q)

            # `φ = 0`, so the Boris push is a rotation and keeps `|v|`; the bound is check 2's,
            # and no parameter of this field is a power of two.
            drift = zero(T)
            for n in 0:N
                v = sol.q[n]
                E = (v[4, 1]^2 + v[5, 1]^2 + v[6, 1]^2) / 2
                drift = max(drift, abs(E - v₀ ⋅ v₀ / 2))
            end
            @test drift ≤ N * eps(T) * (v₀ ⋅ v₀ / 2)
            @test all(isfinite, sol.q[N])
        end

        @testset "the push through the adapter does not allocate" begin
            q, m = T(2), T(3)
            model = VlasovWithExternalFields(
                particle_distribution(T, 3), field; charge = q, mass = m)
            method = SplittingMethod(model, (zero(T), T(1)), T(0.1))
            params = method.equation.parameters
            z̄ = copy(model.distribution.particles.z)
            z = similar(z̄)
            t̄ = zero(T)
            # the two maps are type-stable at the adapter's concrete field type too
            @test @inferred(s_drift!(z, t̄ + T(0.1), z̄, t̄, params)) === nothing
            @test @inferred(s_lorentz!(z, t̄ + T(0.1), z̄, t̄, params)) === nothing
            if Base.JLOptions().code_coverage == 0
                @test allocations_s(s_lorentz!, z, t̄ + T(0.1), z̄, t̄, params) == 0
            end
            @test allocations(collect, z) > 0
        end

        @testset "the concrete types of the adapter" begin
            model = VlasovWithExternalFields(
                particle_distribution(T, 3), field; charge = 2, mass = 3)
            @test all(isconcretetype, fieldtypes(typeof(model)))
            @test all(isconcretetype, fieldtypes(typeof(model.field)))
            # and one level further in: the generated functions are held by a concrete
            # `NamedTuple` of concrete `FieldFunction`s
            @test all(isconcretetype, fieldtypes(typeof((; model.field.functions...))))
            @test eltype(electric_field(field, zero(T), SVector{3, T}(0.1, 0.2, 0.3))) === T
        end

        @testset "the two vector fields of the adapter" begin
            # The right-hand side of each map, against values written out here, at the adapter's
            # own field type and with a nonunit `q/m`. The orbit and energy checks cannot see a
            # zero `ẋ` in the drift or a dropped `qm`: both leave them unchanged. The perturbation
            # makes `E` nonzero as well, which the equilibrium alone leaves at zero.
            q, m, E₀ = T(2), T(3), T(0.3)
            perturbed = FieldFunctions(ABCEquilibrium(a, b, c), EzCosZPerturbation(E₀))
            dist = particle_distribution(T, 2)
            dist.particles.x .= SVector{3, T}(0.31, -0.52, 0.77)
            dist.particles.v .= SVector{3, T}(0.7, 0.4, 0.25)
            model = VlasovWithExternalFields(dist, perturbed; charge = q, mass = m)
            method = SplittingMethod(model, (zero(T), T(1)), T(0.1))
            params = method.equation.parameters
            z = copy(model.distribution.particles.z)
            ż = similar(z)
            t = T(0.4)
            x = SVector{3, T}(0.31, -0.52, 0.77)
            v = SVector{3, T}(0.7, 0.4, 0.25)
            acc = (q / m) *
                  (electric_field(perturbed, t, x) + v × magnetic_field(perturbed, t, x))

            @test @inferred(v_drift!(ż, t, z, params)) === nothing
            @test @inferred(v_lorentz!(ż, t, z, params)) === nothing

            v_drift!(ż, t, z, params)
            for i in axes(z, 2)
                @test ż[1:3, i] == v
                @test all(iszero, ż[4:6, i])
            end

            v_lorentz!(ż, t, z, params)
            for i in axes(z, 2)
                @test all(iszero, ż[1:3, i])
                @test ż[4:6, i] == collect(acc)
            end
        end
    end

    @testset "the chart of the adapter" begin
        # The adapter is confined to Cartesian charts by dispatch: `B♯` and `E♯` are the
        # contravariant components, which are the Cartesian ones only where the chart is the
        # identity map, so the first kick of a curvilinear field raises its `MethodError`.
        curvilinear = FieldFunctions(
            AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 2.0))
        for T in (Float64, Float32)
            x = SVector{3, T}(0.31, -0.52, 0.77)
            @test_throws MethodError magnetic_field(curvilinear, zero(T), x)
            @test_throws MethodError electric_field(curvilinear, zero(T), x)
        end
    end
end
