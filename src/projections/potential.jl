
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
    b = basis(potential)
    ρ = PoissonSolvers.rhs(potential)
    ρ .= 0

    vals = zeros(eltype(ρ), local_width(b))
    points = distribution.particles.x
    weights = distribution.particles.w

    for (x, w) in zip(points, weights)
        j₀ = evaluate_all!(vals, b, x)
        for (t, value) in pairs(vals)
            # Skipping the zeros drops the padding described above. It is not a bounds guard:
            # a nonzero value at an out-of-range index still raises, rather than being folded
            # silently onto another basis function.
            iszero(value) && continue
            ρ[basis_index(b, j₀ + t - 1)] += w * value
        end
    end

    return potential
end
