struct CollisionEntropy{XD, VD, DT <: DistributionFunction{XD, VD}} <: Entropy
    dist::DT
    # entropy::ET 

    function CollisionEntropy(dist::DistributionFunction{XD, VD}) where {XD, VD}
        new{XD, VD, typeof(dist)}(dist)
    end
end

@doc raw"""
    (ent::CollisionEntropy)()

The discrete entropy

```math
S_h = \int_\Omega f_s \, \log f_s \, dv ,
```

evaluated on the Gauß-Legendre grid the spline basis already carries, in any number of velocity
dimensions.

!!! note "The sign is the manuscripts' `S`, and this used not to work at all"
    Both manuscripts define ``S = -\int f \log f`` in the continuum and then report
    ``\int f_s \log f_s`` in their figures; this returns the latter, so it *decreases*
    monotonically. Read as "entropy" with the manuscripts' continuum sign it is ``-S``, and a
    plot of it labelled "entropy" and expected to grow is inverted.

    The earlier implementation routed a one-dimensional spline through a hard-coded
    two-dimensional quadrature that built `SVector{2}` sample points, so a `VD == 1`
    distribution could not be evaluated at them at all; it also took `nquad` as a positional
    default rather than a keyword and applied `log` without checking positivity. There was
    consequently no working entropy diagnostic, which is why every entropy computation in
    `scripts/` is commented out.
"""
# One method per velocity dimensionality, because the two quadratures report their grid
# differently and must not be conflated: a `SplineQuadrature` returns the nodes as a plain
# vector, a `TensorProductQuadrature` as a tuple of per-axis vectors. Splatting the former into
# `Iterators.product` produces one point of as many components as there are nodes.
function (ent::CollisionEntropy{XD, 1})() where {XD}
    sdist = ent.dist
    q = sdist.quadrature
    fs = sdist.spline
    x = quadrature_nodes(q)
    w = quadrature_weights(q)

    S = zero(eltype(sdist))
    for r in eachindex(x)
        f = fs(x[r])
        _check_positive(f, x[r])
        S += w[r] * f * log(f)
    end
    return S
end

function (ent::CollisionEntropy{XD, VD})() where {XD, VD}
    sdist = ent.dist
    q = sdist.quadrature
    fs = sdist.spline
    X = quadrature_nodes(q)
    W = quadrature_weights(q)

    S = zero(eltype(sdist))
    for idx in CartesianIndices(map(length, X))
        pt = ntuple(k -> X[k][idx[k]], VD)
        f = fs(pt)
        _check_positive(f, pt)
        S += prod(ntuple(k -> W[k][idx[k]], VD)) * f * log(f)
    end
    return S
end

function _check_positive(f, v)
    f > 0 || throw(ErrorException(
        "the projected distribution is non-positive, f_s = $(f) at v = $(v), so the " *
        "entropy ∫ f log f dv is undefined there"))
    return nothing
end

# ## TODO: add functions for computing the entropy given a distribution
# function compute_entropy!(entropy, dist <: DistributionFunction{XD, VD}) where {XD, VD}

# end
