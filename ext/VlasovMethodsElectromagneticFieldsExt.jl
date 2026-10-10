module VlasovMethodsElectromagneticFieldsExt

using ElectromagneticFields
using StaticArrays

import VlasovMethods: electric_field, magnetic_field

# The adapter from a `FieldFunctions` to the two functions of VlasovMethods' field interface. It
# is confined to Cartesian charts: `B♯` and `E♯` are the contravariant components, which are the
# Cartesian ones only where the chart is the identity map. In a curvilinear chart they are the
# components in that chart's basis, and a kick built from them would integrate another model —
# the same reason `sodeproblem` of ChargedParticleDynamics refuses those charts. The bound is
# therefore by dispatch, and a curvilinear `FieldFunctions` has no method of either function, so
# the first kick raises its `MethodError`.
#
# `FieldFunctions` is a concrete type with six parameters, and the two important ones are the
# first two; the rest are left free. Both perturbations that a Cartesian chart admits are in the
# union: `ZeroPerturbation`, which the constructor takes as its default and which is not a
# `CartesianPerturbation`, and the Cartesian perturbations such as `EzCosZPerturbation`.
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
