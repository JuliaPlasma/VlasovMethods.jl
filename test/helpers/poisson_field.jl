# The Poisson-field fixture the electric-field test files share. It is included by
# `test/electric_field.jl`, `test/particles/electric_field.jl` and `test/quality/jet.jl`, and is
# listed in no `runtests.jl`.

using PoissonSolvers
using Random
using VlasovMethods: PoissonField

# The field of `n` seeded positions in `(0.2, 0.8)` with equal weights, on `ncells` periodic
# spline cells of degree `order`.
function poisson_field(::Type{T}; n = 40, ncells = 16, order = 3) where {T}
    Random.seed!(1234)
    q = T(0.2) .+ T(0.6) .* rand(T, n)
    w = fill(T(1) / n, n)
    potential = Potential(PeriodicBasisSpline((zero(T), one(T)), order, ncells))
    return PoissonField(potential), q, w
end
