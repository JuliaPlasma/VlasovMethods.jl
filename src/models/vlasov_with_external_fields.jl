"""
    electric_field(field, t, x)
    magnetic_field(field, t, x)

The electric and magnetic field of `field` at the time `t` and the position `x`, each an
`SVector{3,T}` in Cartesian components, with `T` the element type of `x`.

The two functions are the interface between a [`VlasovWithExternalFields`](@ref) model and the
field its particles move in. The model calls them once per particle and per kick, with `x` the
particle's position and `t` the time the kick starts at, and a field type defines one method of
each. There is no abstract field type and no fallback method, so a field that defines neither
raises its `MethodError` at the first kick inside `integrate`, which is where a position exists.

[`UniformField`](@ref) is the one field type in this package. The `ElectromagneticFields`
extension adds methods for that package's `FieldFunctions` of a Cartesian chart; see
[`VlasovWithExternalFields`](@ref).
"""
function electric_field end

"""
    magnetic_field(field, t, x)

The magnetic field of `field` at the time `t` and the position `x`. See [`electric_field`](@ref)
for the interface the two functions share.
"""
function magnetic_field end

"""
    UniformField(E, B)

A constant electric and magnetic field, `E` and `B` in Cartesian components, stored as
`SVector{3,T}` with `T` the promoted element type of the two arguments. It is constant in `x` and
in `t`, so it is the field of the analytic cases: a gyro-orbit at `E = 0`, and a constant
acceleration along `B` at `E ∥ B`.

`E` and `B` are two vectors of length three, or two `SVector{3}`s; a length other than three
throws the `DimensionMismatch` of `SVector{3}`. There is no method taking a scalar, so a constant
field of `Int` coordinates is written `UniformField(SVector(0, 0, 0), SVector(0, 0, 0))`.
"""
struct UniformField{T}
    E::SVector{3, T}
    B::SVector{3, T}
end

function UniformField(E, B)
    T = promote_type(eltype(E), eltype(B))
    UniformField(SVector{3, T}(E), SVector{3, T}(B))
end

# The field's own precision is not the particle's: the interface returns the element type of the
# coordinates, so both methods convert.
@inline function electric_field(field::UniformField, t, x::SVector{3, T}) where {T}
    return SVector{3, T}(field.E)
end

@inline function magnetic_field(field::UniformField, t, x::SVector{3, T}) where {T}
    return SVector{3, T}(field.B)
end

"""
    VlasovWithExternalFields(dist, field; charge = one(T), mass = one(T))

A `VlasovModel` of one particle species in an external electromagnetic field: the
particles of `dist`, a `ParticleDistribution{T,3,3}`, pushed by `field`. `charge` and `mass` are
the species' charge and mass, in units in which the equation of motion reads
`ẍ = (q/m)(E + ẋ × B)`, and both are stored as `T`. The constructor throws an `ArgumentError`
unless the converted mass is positive.

The field is any value that defines [`electric_field`](@ref) and [`magnetic_field`](@ref) for an
`SVector{3}` coordinate; the model puts no bound on the type of its field parameter, so the two
functions are the whole contract. `UniformField` is the field type this package ships, and the
`ElectromagneticFields` extension adds the two methods for that package's `FieldFunctions` of a
Cartesian chart — `FieldFunctions(ABCEquilibrium(a, b, c))` is passed to the constructor as it
is. In a curvilinear chart the components are not Cartesian, so the extension defines no method
and the first kick raises its `MethodError`.

The method is `SplittingMethod(model, tspan, tstep)`: a Strang splitting of the drift `ẋ = v` and
the Lorentz force `v̇ = (q/m)(E + v × B)`, with the rotational part of the force a Boris push. One
step is a half drift, two half Boris pushes with the field evaluated at the step's start, and a
half drift. `GeometricIntegrators.integrate(method.integrator)` and
`run!(method, h5file)` serve it as they serve the Vlasov–Poisson model.
"""
struct VlasovWithExternalFields{T, DT <: ParticleDistribution{T, 3, 3}, FT} <: VlasovModel
    distribution::DT
    field::FT
    charge::T
    mass::T

    function VlasovWithExternalFields(dist::ParticleDistribution{T, 3, 3}, field,
            charge::Real, mass::Real) where {T}
        # The guard reads the converted mass: `1e-50` is not positive in `Float32`, where it
        # converts to `0f0`, although it is in `Float64`. `NaN` fails the comparison and is
        # rejected with it; `Inf` passes.
        m = T(mass)
        m > 0 || throw(ArgumentError("the particle mass must be positive, got $mass"))
        new{T, typeof(dist), typeof(field)}(dist, field, T(charge), m)
    end
end

function VlasovWithExternalFields(dist::ParticleDistribution{T, 3, 3}, field;
        charge = one(T), mass = one(T)) where {T}
    VlasovWithExternalFields(dist, field, charge, mass)
end

#################################################################
# Define Splitting Method for the external-field Vlasov model   #
#################################################################

# The state is `model.distribution.particles.z`, the positions in its first three rows and the
# velocities in its last three. The fields below write all six rows, so no substep leaves a row
# of the previous state behind.

# Vector field for the drift, `(ẋ, v̇) = (v, 0)`.
function v_drift!(ż, t, z, params)
    for i in axes(ż, 2)
        ż[1, i] = z[4, i]
        ż[2, i] = z[5, i]
        ż[3, i] = z[6, i]
        ż[4, i] = 0
        ż[5, i] = 0
        ż[6, i] = 0
    end
end

# Vector field for the Lorentz force, `(ẋ, v̇) = (0, (q/m)(E + v × B))`.
function v_lorentz!(ż, t, z, params)
    qm = params.qm
    for i in axes(ż, 2)
        x = SVector{3, eltype(z)}(z[1, i], z[2, i], z[3, i])
        v = SVector{3, eltype(z)}(z[4, i], z[5, i], z[6, i])
        a = qm *
            (electric_field(params.field, t, x) + v × magnetic_field(params.field, t, x))
        ż[1, i] = 0
        ż[2, i] = 0
        ż[3, i] = 0
        ż[4, i] = a[1]
        ż[5, i] = a[2]
        ż[6, i] = a[3]
    end
end

# Solution for the drift: `x ← x̄ + h v̄`, `v` unchanged.
function s_drift!(z, t, z̄, t̄, params)
    h = t - t̄
    for i in axes(z, 2)
        z[1, i] = z̄[1, i] + h * z̄[4, i]
        z[2, i] = z̄[2, i] + h * z̄[5, i]
        z[3, i] = z̄[3, i] + h * z̄[6, i]
        z[4, i] = z̄[4, i]
        z[5, i] = z̄[5, i]
        z[6, i] = z̄[6, i]
    end
end

# Solution for the Lorentz force: one Boris push with `h = t − t̄`, the field read at `t̄` and at
# the substep's own position. `v⁻`, `τ`, `s`, `v′` and `v⁺` are the standard names of the push.
# With `E = 0` it is a rotation of `v` about `B`, which keeps `|v|` for any `B`, constant or not.
function s_lorentz!(z, t, z̄, t̄, params)
    h = t - t̄
    qh = params.qm * h / 2
    for i in axes(z, 2)
        x = SVector{3, eltype(z)}(z̄[1, i], z̄[2, i], z̄[3, i])
        E = electric_field(params.field, t̄, x)
        B = magnetic_field(params.field, t̄, x)
        v⁻ = SVector{3, eltype(z)}(z̄[4, i], z̄[5, i], z̄[6, i]) + qh * E
        τ = qh * B
        s = 2 * τ / (1 + τ ⋅ τ)
        v⁺ = v⁻ + (v⁻ + v⁻ × τ) × s
        v = v⁺ + qh * E
        z[1, i] = z̄[1, i]
        z[2, i] = z̄[2, i]
        z[3, i] = z̄[3, i]
        z[4, i] = v[1]
        z[5, i] = v[2]
        z[6, i] = v[3]
    end
end

# Constructor for a splitting method from GeometricIntegrators, as for the Vlasov–Poisson model:
# one solution step pushes all particles.
function SplittingMethod(
        model::VlasovWithExternalFields, tspan::Tuple, tstep::Real)
    # collect parameters
    params = (field = model.field, qm = model.charge / model.mass)

    # create geometric problem
    equ = GeometricEquations.SODEProblem(
        (v_drift!, v_lorentz!),
        (s_drift!, s_lorentz!),
        tspan, tstep, copy(model.distribution.particles.z);
        parameters = params)

    # create integrator
    int = GeometricIntegrators.GeometricIntegrator(equ, GeometricIntegrators.Strang())

    # put together splitting method
    SplittingMethod(model, equ, int)
end
