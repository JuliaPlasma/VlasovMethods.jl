struct CollisionEntropy{XD, VD, DT <: DistributionFunction{<:Any, XD, VD}} <: Entropy
    dist::DT
    # entropy::ET 

    function CollisionEntropy(dist::DistributionFunction{<:Any, XD, VD}) where {XD, VD}
        new{XD, VD, typeof(dist)}(dist)
    end
end

# One method per velocity dimensionality, because the two quadratures report their grid
# differently and must not be conflated: a `SplineQuadrature` returns the nodes as a plain
# vector, a `TensorProductQuadrature` as a tuple of per-axis vectors. Splatting the former into
# `Iterators.product` produces one point of as many components as there are nodes.
@doc raw"""
    (ent::CollisionEntropy)()

The discrete entropy

```math
S_h = \int_\Omega f_s \, \log |f_s| \, dv = \frac{1}{2} \int_\Omega f_s \, \log (f_s^2) \, dv ,
```

`eq:disc_entropy` of the Landau and Lenard-Bernstein metriplectic manuscripts, evaluated on the
Gauß-Legendre grid the spline basis already carries, in any number of velocity dimensions.

``S_h`` is the H-functional of the manuscripts' ``S = \int f \log f``: it decreases under the
collision operators. ``\log |f_s|`` equals ``\log (f_s)`` where ``f_s > 0``, and stays
defined where the ``L^2`` projection of a particle distribution undershoots below zero. A node
where ``f_s`` is zero or not finite throws a `DomainError`, because ``\log |f_s|`` is undefined
there.
"""
function (ent::CollisionEntropy{XD, 1})() where {XD}
    sdist = ent.dist
    q = sdist.quadrature
    fs = sdist.spline
    x = quadrature_nodes(q)
    w = quadrature_weights(q)

    S = zero(eltype(sdist))
    for r in eachindex(x)
        f = fs(x[r])
        _check_log_defined(f, x[r])
        S += w[r] * f * log(abs(f))
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
        _check_log_defined(f, pt)
        S += prod(ntuple(k -> W[k][idx[k]], VD)) * f * log(abs(f))
    end
    return S
end

# log|f| is -Inf where f is zero and undefined where f is not finite; every other value, a
# negative one included, has a finite logarithm. The entropy, 𝕃 of both metriplectic models and
# `compute_entropy` sample log|f_s| at the quadrature nodes and refuse such a node.
function _check_log_defined(f, v)
    (f != 0 && isfinite(f)) || throw(DomainError(f,
        "the projected distribution is zero or not finite at v = $(v), so " *
        "log|f_s| is undefined there"))
    return nothing
end
