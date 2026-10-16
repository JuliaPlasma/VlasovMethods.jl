# The reduced equations of the analytic quadratic model: the Hamiltonian that the forms of
# `ReducedBasisMethods.reduced_ode` and `reduced_pode` take as `∇H!`, its value, and the lift of a
# GeometricIntegrators solution back to `(q, p)`. The two forms themselves are that package's; this
# file supplies `H` and the lifting, and no integrator.
#
# The self-consistent `VlasovHamiltonian` of `src/models/vlasov_poisson.jl` is deliberately not
# wired to `reduced_pode`: with weights its canonical momentum is `w·v` and not `v`.

"""
    AnalyticHamiltonian(β)

The Hamiltonian of the analytic quadratic-potential model, ``H(q, p) = ½|p|² + ½β²|q|²``, for a
state `x = [q; p]` of `2N` entries, and the `∇H!` of the reduced forms of `ReducedBasisMethods`.
The two halves are the two halves of `x`.

The instance is both:

- the gradient, called as `∇H!(g, t, x, params)`, which writes `∇H(x)` — the first block `β²q` and
  the second `p` — into `g` and returns it. It is what `ReducedBasisMethods.reduced_ode` and
  `ReducedBasisMethods.reduced_pode` take as their `∇H!`, and it is not RBM's
  `hamiltonian_gradient!`, which is the reduced gradient of a `SymplecticDEIM`;
- the Hamiltonian, called as `h(t, x, params)`, which returns `H(x)`.

Neither call reads `t` or `params`: the model is autonomous and takes no equation parameters. The
gradient accepts any element type, which Newton's `JacobianAutodiff` needs — the Lobatto IIIA–IIIB(2)
run of the PODE form calls the vector field with `ForwardDiff.Dual` states.

# Example

```julia
h = AnalyticHamiltonian(6.0)
z₀ = project(b, x₀)
prob = reduced_pode(b, h, (0.0, 1000.0), 1e-4, z₀)
```
"""
struct AnalyticHamiltonian{T}
    β::T
end

function (h::AnalyticHamiltonian)(t, x, params)
    N = length(x) ÷ 2
    return @views (h.β^2 * sum(abs2, x[1:N]) + sum(abs2, x[(N + 1):(2N)])) / 2
end

function (h::AnalyticHamiltonian)(g, t, x, params)
    N = length(x) ÷ 2
    β² = h.β^2
    @inbounds for i in 1:N
        g[i] = β² * x[i]
        g[N + i] = x[N + i]
    end
    return g
end

"""
    lift_solution(b::ReducedBasis, sol) -> (q, p)

The full-order positions and velocities of the GeometricIntegrators solution `sol` of the reduced
PODE form in the basis `b`: the two blocks of the lifted trajectory, each of `N` rows and one column
per stored time.

The solution's own halves `sol.q` and `sol.p` are the reduced coordinates, and
`ReducedBasisMethods.lift(b, [q; p])` is `[Φq; Φp]` for a cotangent lift and
`[Φa - Ψb; Ψa + Φb]` for a complex SVD. The two blocks are views of that one lifted array. A POD
basis has no PODE solution and no `q` and `p` halves to split, so it is a `MethodError`.
"""
function lift_solution(
        b::ReducedBasis{<:Any, <:Union{CotangentLiftSVD, CotangentLiftEVD, ComplexSVD}},
        sol::GeometricIntegrators.AbstractSolution)
    X = lift(b, vcat(Array(sol.q), Array(sol.p)))
    N = size(X, 1) ÷ 2
    return @views X[1:N, :], X[(N + 1):end, :]
end
