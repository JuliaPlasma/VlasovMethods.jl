using Aqua
using VlasovMethods
using Test

# The two `Vararg` call operators of `DistributionFunction` are ambiguous for `{T, 0, 0}`.
Aqua.test_all(VlasovMethods; ambiguities = (; broken = true))  # issue #57
