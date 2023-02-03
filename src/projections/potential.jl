
function projection!(out, potential::Potential{<:PeriodicBSplineBasis}, distribution::ParticleDistribution)
    b = Splines.PeriodicVector(out)
    b .= 0

    basis = potential.basis
    points = distribution.particles.x
    weights = distribution.particles.w

    for (x, w) in zip(points, weights)
        ilast, bs = basis(x)  # same as `evaluate_all`
        # Iterate over evaluated basis functions.
        # The indices of the evaluated basis functions are ilast:-1:(ilast - k + 1),
        # where k is the spline order.
        for (δi, bi) ∈ pairs(bs)
            i = ilast + 1 - δi
            b[i] += w * bi
        end
    end

    return out
end
