
"""
    projection!(potential::PoissonSolvers.Potential, distribution::ParticleDistribution)

Deposit the particle weights of `distribution` onto the basis of `potential`, writing the charge
density into the potential's right-hand side buffer.

Each particle contributes to the `local_width(basis)` basis functions that do not vanish at its
position. `evaluate_all!` writes those values into a buffer and returns the index of the first,
before wrapping; `basis_index` is what wraps it, and is the identity where the basis is not
periodic.
"""
function projection!(potential::PoissonSolvers.Potential,
        distribution::ParticleDistribution)
    b = basis(potential)
    ρ = PoissonSolvers.rhs(potential)
    ρ .= 0

    values = zeros(eltype(ρ), local_width(b))
    points = distribution.particles.x
    weights = distribution.particles.w

    for (x, w) in zip(points, weights)
        first = evaluate_all!(values, b, x, 0)
        for (t, value) in pairs(values)
            ρ[basis_index(b, first + t - 1)] += w * value
        end
    end

    return potential
end
