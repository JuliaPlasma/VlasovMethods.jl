using Aqua
using VlasovMethods
using Test

# The two `Vararg` call operators of `DistributionFunction` in `src/distributions/distribution.jl`
# are ambiguous for a `DistributionFunction{T, 0, 0}`: issue #57.
Aqua.test_all(VlasovMethods; ambiguities = (; broken = true))
