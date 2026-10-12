module VlasovMethodsElectromagneticFieldsExt

import ElectromagneticFields
using StaticArrays: SVector

import VlasovMethods: electric_field, magnetic_field

# The adapter from a `FieldFunctions` to the two functions of VlasovMethods' field interface. It
# is confined to Cartesian charts: `B♯` and `E♯` are the contravariant components, which are the
# Cartesian ones only where the chart is the identity map. In a curvilinear chart they are the
# components in that chart's basis, and a kick built from them would integrate another model —
# the same reason `sodeproblem` of ChargedParticleDynamics refuses those charts. The bound is
# therefore by dispatch, and a curvilinear `FieldFunctions` has no method of either function, so
# the first kick raises its `MethodError`.
#
# The alias binds the first two of the six parameters of `FieldFunctions`: a Cartesian
# equilibrium, and either the default `ZeroPerturbation`, which is not a `CartesianPerturbation`,
# or a Cartesian perturbation such as `EzCosZPerturbation`.
const CartesianFieldFunctions = ElectromagneticFields.FieldFunctions{
    <:ElectromagneticFields.CartesianEquilibrium,
    <:Union{ElectromagneticFields.ZeroPerturbation,
        ElectromagneticFields.CartesianPerturbation}}

# The accessors already return the element type of the coordinates — `FieldFunction` converts to
# `float(eltype(ξ))` — so the two methods return what they get, and a `Float32` particle sees a
# `Float32` field whatever precision the equilibrium's parameters were written in.
@inline function electric_field(field::CartesianFieldFunctions, t, x::SVector{3})
    return ElectromagneticFields.E♯(field, t, x)
end

@inline function magnetic_field(field::CartesianFieldFunctions, t, x::SVector{3})
    return ElectromagneticFields.B♯(field, t, x)
end

end
