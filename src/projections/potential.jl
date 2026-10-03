
"""
    projection!(potential::PoissonSolvers.Potential, distribution::ParticleDistribution)

Deposit the particle weights of `distribution` onto the basis of `potential`, writing the charge
density into the potential's right-hand side buffer.

Each particle contributes to the `local_width(basis)` basis functions that do not vanish at its
position. `evaluate_all!` reduces a periodic argument into the domain and writes those values
into a buffer, returning the index of the first; `basis_index` is what wraps it, and is the
identity where the basis is not periodic. The particle state keeps the unwrapped position, so
the reduction belongs to the evaluation.

A recombined basis has no single width: a function near an end spans the union of two parent
supports, so `local_width` is the largest block and a cell with fewer nonzero functions leaves
the tail of the buffer zero. Those padding entries carry indices past `nbasis(basis)`, which is
why the zeros are skipped rather than added.
"""
function projection!(potential::PoissonSolvers.Potential,
        distribution::ParticleDistribution)
    # the same `local_width(basis)` buffer that `VlasovPoisson` holds, allocated here per call
    work = _local_buffers(basis(potential), eltype(PoissonSolvers.rhs(potential)))
    _deposit!(potential, work, distribution.particles.x, distribution.particles.w)
end

# The deposit of `projection!`, from the positions in the first row of `x` and the weights in
# the first row of `w`, with `vals` as the buffer of basis values. The first row is read by
# index, so the positions of a state matrix `z` are deposited without a slice. The reduction
# into the domain happens in `evaluate_all!`, so an unwrapped position needs no wrap here.
function _deposit!(potential::PoissonSolvers.Potential, vals::AbstractVector,
        x::AbstractMatrix, w::AbstractMatrix)
    b = basis(potential)
    ρ = PoissonSolvers.rhs(potential)
    ρ .= 0

    for i in axes(x, 2)
        j₀ = evaluate_all!(vals, b, x[1, i])
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
