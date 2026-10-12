using ElectromagneticFields
using GeometricIntegrators: GeometricIntegrators
using StaticArrays
using Test
using VlasovMethods

include("../helpers/external_fields.jl")

const EXTENSION = :VlasovMethodsElectromagneticFieldsExt
const EF_UUID = "d6c1ba6f-ee03-53af-b876-68cefeb88ec8"

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
        # `ElectromagneticFields` is a dependency of this test environment, so this process has
        # loaded the extension; a fresh process, in the environment `Pkg.test` built, shows that
        # `using VlasovMethods` alone loads neither the package nor the extension.
        @test success(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $CHILD`)
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
            # `φ = 0`, so the Boris push is a rotation and keeps `|v|`; no parameter of this field
            # is a power of two
            test_energy(sol, v₀)
            @test all(isfinite, sol.q[end])
        end

        @testset "the push through the adapter does not allocate" begin
            test_push_does_not_allocate(T, field)
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
            # at the adapter's own field type; the perturbation makes `E` nonzero, which the
            # equilibrium alone leaves at zero
            perturbed = FieldFunctions(ABCEquilibrium(a, b, c), EzCosZPerturbation(T(0.3)))
            test_vector_fields(T, perturbed, SVector{3, T}(0.31, -0.52, 0.77),
                SVector{3, T}(0.7, 0.4, 0.25))
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
