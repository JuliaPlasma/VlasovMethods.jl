@doc raw"""
Projection of a particle distribution onto a spline basis, and of a prescribed function onto
the same basis.

The particle projection is the ``L^2`` projection of the Klimontovich distribution
``f_p(v) = \sum_\alpha w_\alpha \, \delta(v - v_\alpha)``, i.e. the coefficients

```math
\hat{f}_I = \sum_K \mathbb{M}^{-1}_{IK} \sum_{\alpha \,:\, v_\alpha \in \operatorname{supp} \varphi_K}
            w_\alpha \, \varphi_K(v_\alpha) ,
```

which is `eq:projection_coeffs` of the Landau manuscript and `eq:projection_f` of the
Lenard-Bernstein one. Pairing a delta with a basis function turns the integral into an
evaluation, so the load vector is a deposition: each particle contributes to the
``\prod_d (p_d+1)`` basis functions whose support contains it, and to no others.
"""

# The block of basis functions overlapping one point, deposited with weight `w`. Written
# against `evaluate_all!`, so it touches the (p+1) functions that are actually nonzero rather
# than looping over the whole basis: that is the difference between O(N_p p) and O(N_p N) for
# the whole deposition, and it is the innermost loop of every step.
function _deposit!(rhs::AbstractVector, b::AbstractBSplineBasis, v::Number, w,
        buf::AbstractVector)
    N = nbasis(b)
    j₀ = evaluate_all!(buf, b, v)
    for t in eachindex(buf)
        # `basis_index` wraps where the basis is periodic and is the identity where it is not.
        # Writing `mod1` at the call site instead is what made the earlier implementation
        # inconsistent: the spline was evaluated with the wrap and its gradient without it.
        j = basis_index(b, j₀ + t - 1)
        1 ≤ j ≤ N || continue
        rhs[j] += w * buf[t]
    end
    return rhs
end

function _deposit!(rhs::AbstractArray{RT, D}, B::TensorProductBasis{T, D}, v, w,
        bufs::NTuple{D, <:AbstractVector}) where {RT, T, D}
    j₀ = evaluate_all!(bufs, B, v)
    for t in CartesianIndices(ntuple(k -> length(bufs[k]), D))
        val = w
        inrange = true
        idx = ntuple(D) do k
            val *= bufs[k][t[k]]
            j = basis_index(B.bases[k], j₀[k] + t[k] - 1)
            inrange &= (1 ≤ j ≤ size(B, k))
            j
        end
        # The bound is checked on every component separately. Testing only the flattened
        # index, as the earlier implementation did, does not bound the components: with
        # i = 0, j = 3 and M = 10 the flat index (j-1)M + i = 20 passes a `1 ≤ k ≤ M²` test
        # and decodes to (10, 2) — an out-of-support contribution aliased onto a valid and
        # entirely unrelated basis function rather than discarded.
        inrange || continue
        rhs[idx...] += val
    end
    return rhs
end

_local_buffers(b::AbstractBSplineBasis, ::Type{VT}) where {VT} = zeros(VT, local_width(b))
function _local_buffers(B::TensorProductBasis{T, D}, ::Type{VT}) where {T, D, VT}
    ntuple(k -> zeros(VT, local_width(B.bases[k])), D)
end

# Whether a point lies in the domain the basis resolves. A periodic axis has no outside: an
# argument is reduced onto the domain before evaluation, so every particle contributes.
_axis_in_domain(b::PeriodicBSplineBasis, x) = true
_axis_in_domain(b::AbstractBSplineBasis, x) = x ∈ domain(b)

_in_domain(b::AbstractBSplineBasis, v::Number) = _axis_in_domain(b, v)
function _in_domain(B::TensorProductBasis{T, D}, v) where {T, D}
    all(k -> _axis_in_domain(B.bases[k], v[k]), 1:D)
end

@doc raw"""
    projection(velocities, dist::ParticleDistribution, final_dist::SplineDistribution)

Project the particles onto the basis of `final_dist`, writing the coefficients into it and
returning its `Spline`.

`velocities` is a vector for a one-dimensional velocity space and a `VD × N` matrix otherwise;
it is passed separately from `dist` because the implicit solvers evaluate the right-hand side
at a stage value rather than at the particles' stored positions.

!!! warning "A particle outside the domain is an error, not a warning"
    The velocity domain is truncated and the true one is not, so a particle can leave. When it
    does there is nothing sensible to do with it: it contributes nothing to the projection, so
    its weight is lost from ``f_s`` while it keeps moving, and `f_s` then evaluates to zero at
    its position — which makes the ``f_s'/f_s`` of every collision operator here a division by
    zero. The Lenard-Bernstein manuscript says as much in its closing remarks.

    The earlier implementation printed a warning per offending particle in the two-dimensional
    branch and dropped it silently in the one-dimensional one. This throws, naming the count
    and the domain, because a run that has lost mass is not a run whose output means anything.
    Enlarge the domain, or give the mesh oversized end cells with `length_big_cell`.
"""
function projection(velocities::AbstractVector{VT},
        dist::ParticleDistribution{PT, XD, 1},
        final_dist::SplineDistribution{ST, XD, 1}) where {VT, PT, ST, XD}
    b = final_dist.basis
    rhs = zeros(VT, nbasis(b))
    buf = _local_buffers(b, VT)
    w = dist.particles.w

    escaped = 0
    for p in eachindex(velocities)
        if !_in_domain(b, velocities[p])
            escaped += 1
            continue
        end
        _deposit!(rhs, b, velocities[p], w[1, p], buf)
    end
    _check_escaped(escaped, length(velocities), b)

    mass_solve!(final_dist.coefficients, mass_operator(final_dist), rhs)
    return final_dist.spline
end

function projection(velocities::AbstractMatrix{VT},
        dist::ParticleDistribution{PT, XD, VD},
        final_dist::SplineDistribution{ST, XD, VD}) where {VT, PT, ST, XD, VD}
    B = final_dist.basis
    rhs = zeros(VT, size(B)...)
    bufs = _local_buffers(B, VT)
    w = dist.particles.w

    escaped = 0
    for p in axes(velocities, 2)
        v = view(velocities, :, p)
        if !_in_domain(B, v)
            escaped += 1
            continue
        end
        _deposit!(rhs, B, v, w[1, p], bufs)
    end
    _check_escaped(escaped, size(velocities, 2), B)

    mass_solve!(final_dist.coefficients, mass_operator(final_dist), rhs)
    return final_dist.spline
end

function _check_escaped(escaped::Integer, total::Integer, basis)
    escaped == 0 && return nothing
    throw(ErrorException(
        "$(escaped) of $(total) particles left the velocity domain $(domain(basis)). " *
        "Their weight would be dropped from the spline projection while they keep moving, " *
        "and f_s would then be zero at their position, making f_s'/f_s a division by zero. " *
        "Enlarge the domain, or give the mesh oversized end cells with `length_big_cell`."))
end

@doc raw"""
    project_function(f, sdist::SplineDistribution)

The ``L^2`` projection of the function `f` onto the basis of `sdist`, written into its
coefficients.

`f` is called with a scalar for a one-dimensional velocity space and with a `VD`-tuple
otherwise. This is `eq:numerical-example-1-coefficients` of the Lenard-Bernstein manuscript,
``\hat{f} = \mathbb{M}^{-1} \int f \varphi_i \, dv``, and it works in any number of velocity
dimensions rather than only in two.
"""
function project_function(f, sdist::SplineDistribution)
    l2_projection!(sdist.coefficients, sdist.quadrature, f)
    return sdist.spline
end

"""
    project_Maxwellian(sdist::SplineDistribution)

Project the normalised Maxwellian onto the basis of `sdist`.

The projection goes through the Gauß-Legendre quadrature the basis already carries, in any
number of velocity dimensions. The previous implementation set up an `IntegralProblem` per
coefficient and called `Integrals.solve` with `HCubatureJL` — neither of which was imported,
so it raised `UndefVarError` whenever it was reached — and existed only for `VD == 2`.
"""
function project_Maxwellian(sdist::SplineDistribution{T, XD, 1}) where {T, XD}
    project_function(v -> MaxwellianDistribution(SVector(v)), sdist)
end

function project_Maxwellian(sdist::SplineDistribution{T, XD, VD}) where {T, XD, VD}
    project_function(v -> MaxwellianDistribution(SVector(v)), sdist)
end

function projection!(init_dist::SplineDistribution, final_dist::ParticleDistribution)
    error("projecting a spline distribution back onto particles is not implemented")
end
