
"""
    projection!(potential::PoissonSolvers.Potential, distribution::ParticleDistribution)

Deposit the particle weights of `distribution` onto the basis of `potential`, writing the charge
density into the potential's right-hand side buffer.

Each particle contributes to the `local_width(basis)` basis functions that do not vanish at its
position. `evaluate_all!` writes those values into a buffer and returns the index of the first,
before wrapping; `basis_index` is what wraps it, and is the identity where the basis is not
periodic.

A recombined basis has no single width: a function near an end spans the union of two parent
supports, so `local_width` is the largest block and a cell with fewer nonzero functions leaves
the tail of the buffer zero. Those padding entries carry indices past `nbasis(basis)`, which is
why the zeros are skipped rather than added.
"""
function projection!(potential::PoissonSolvers.Potential,
        distribution::ParticleDistribution)
    _deposit!(potential, _work(potential), distribution.particles.x, distribution.particles.w)
end

# A work vector of length `local_width(basis)` for the deposit and the field evaluation.
# `projection!` allocates one per call; `VlasovPoisson` holds one, so its RHS does not.
function _work(potential::PoissonSolvers.Potential)
    zeros(
        eltype(PoissonSolvers.rhs(potential)), local_width(basis(potential)))
end

# `x` reduced into the domain of a periodic basis `b`. The particle state keeps the unwrapped
# position, so every evaluation on the basis reduces it. Any other basis takes `x` unchanged.
function _reduce_into_domain(b::PeriodicBSplineBasis, x::Number)
    a = minimum(domain(b))
    return a + mod(x - a, maximum(domain(b)) - a)
end
_reduce_into_domain(b, x::Number) = x

# The deposit of `projection!`, from the positions in the first row of `x` and the weights in
# the first row of `w`, with `vals` as the buffer of basis values. The first row is read by
# index, so the positions of a state matrix `z` are deposited without a slice.
function _deposit!(potential::PoissonSolvers.Potential, vals::AbstractVector,
        x::AbstractMatrix, w::AbstractMatrix)
    b = basis(potential)
    ρ = PoissonSolvers.rhs(potential)
    ρ .= 0

    for i in axes(x, 2)
        j₀ = evaluate_all!(vals, b, _reduce_into_domain(b, x[1, i]))
        for (t, value) in pairs(vals)
            # Skipping the zeros drops the padding described in `projection!`. It is not a
            # bounds guard: a nonzero value at an out-of-range index still raises, rather than
            # being folded silently onto another basis function.
            iszero(value) && continue
            ρ[basis_index(b, j₀ + t - 1)] += w[1, i] * value
        end
    end

    return potential
end
