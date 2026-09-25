@doc raw"""
    GridDistribution{DT, XT, VT, FT}

The distribution function of a 1D1V phase space, stored as its values at the nodes of an
`nx × nv` grid. It is a `DistributionFunction{DT, 1, 1}`, beside `ParticleDistribution` and
`SplineDistribution`.

# The grid

  - `x` is **periodic**: ``x_i = x_a + (i-1) h_x`` for ``i = 1, \dots, n_x``, with
    ``h_x = (x_b - x_a) / n_x``. The node ``x_b`` is the node ``x_a`` again and is not stored.
  - `v` is **bounded**: ``v_j = v_a + (j-1) h_v`` for ``j = 1, \dots, n_v``, with
    ``h_v = (v_b - v_a) / (n_v - 1)``. Both ends of the range are nodes.

`values` is the `nx × nv` matrix ``f_{ij} = f(x_i, v_j)``. Its column-major linear order, `x`
running fastest, is the order the grid stencils such as `_apply_∫dv!` use, so
`vec(dist.values)` is their state vector.

# Evaluation

`dist(x, v)` interpolates bilinearly between the four surrounding nodes. It is periodic in `x`
and zero outside ``[v_a, v_b]``.

# Moments

[`velocity_moments`](@ref) integrates over `v` by the rectangle rule,
``\int f(x_i, v) \, dv \approx \sum_j f_{ij} h_v``. This is the integral of the interpolant plus
``\tfrac{h_v}{2} (f_{i1} + f_{i n_v})``, so it is exact up to the values at the two ends of the
velocity range and first order in ``h_v`` where those do not vanish.
"""
struct GridDistribution{DT, XT <: AbstractRange{DT}, VT <: AbstractRange{DT},
    FT <: AbstractMatrix{DT}} <: DistributionFunction{DT, 1, 1}
    x::XT
    v::VT
    values::FT

    function GridDistribution(x::XT, v::VT,
            values::FT) where {DT, XT <: AbstractRange{DT}, VT <: AbstractRange{DT},
            FT <: AbstractMatrix{DT}}
        length(v) ≥ 2 || throw(ArgumentError(
            "a bounded velocity grid needs at least two nodes, got $(length(v))"))
        size(values) == (length(x), length(v)) || throw(DimensionMismatch(
            "values of size $(size(values)) on a grid of $(length(x)) × $(length(v)) nodes"))
        new{DT, XT, VT, FT}(x, v, values)
    end
end

"""
    GridDistribution(values::AbstractMatrix, xdomain, vdomain)

Wrap the node values `values` of size `nx × nv` on the periodic interval `xdomain` and the
bounded interval `vdomain`. Each domain is anything with a `first` and a `last`.
"""
function GridDistribution(values::AbstractMatrix{DT}, xdomain, vdomain) where {DT}
    nx, nv = size(values)
    x = range(DT(first(xdomain)), DT(last(xdomain)); length = nx + 1)[1:nx]
    v = range(DT(first(vdomain)), DT(last(vdomain)); length = nv)
    GridDistribution(x, v, values)
end

"""
    GridDistribution(nx, nv, xdomain, vdomain)

A `GridDistribution` of `Float64` zeros on `nx × nv` nodes.
"""
function GridDistribution(nx::Integer, nv::Integer, xdomain, vdomain)
    GridDistribution(zeros(nx, nv), xdomain, vdomain)
end

Base.eltype(::GridDistribution{DT}) where {DT} = DT
Base.size(dist::GridDistribution) = size(dist.values)
Base.length(dist::GridDistribution) = length(dist.values)

function Base.show(io::IO, dist::GridDistribution{DT}) where {DT}
    print(io, "GridDistribution{", DT, "}(", size(dist.values), ")")
end

## Evaluation

function evaluate(dist::GridDistribution{DT}, x, v) where {DT}
    nx, nv = size(dist.values)
    first(dist.v) ≤ only(v) ≤ last(dist.v) || return zero(DT)
    η = (only(v) - first(dist.v)) / step(dist.v)
    ξ = (only(x) - first(dist.x)) / step(dist.x)

    ξ₀ = floor(ξ)
    t = ξ - ξ₀
    i = mod(Int(ξ₀), nx) + 1
    i⁺ = mod(i, nx) + 1

    # `v = v_b` falls in the last cell, at its upper end
    j = min(floor(Int, η), nv - 2) + 1
    s = η - (j - 1)

    f = dist.values
    (1 - t) * (1 - s) * f[i, j] + t * (1 - s) * f[i⁺, j] +
    (1 - t) * s * f[i, j + 1] + t * s * f[i⁺, j + 1]
end

## Moments

@doc raw"""
    velocity_moments(dist::GridDistribution)

The velocity moments of `dist` at each `x`-node, as the named tuple `(density, momentum,
energy)` of three vectors of length `nx`:

```math
n_i = \sum_j f_{ij} h_v , \qquad
p_i = \sum_j v_j f_{ij} h_v , \qquad
\varepsilon_i = \frac{1}{2} \sum_j v_j^2 f_{ij} h_v .
```

`energy` is the kinetic energy density, with the factor ``\tfrac{1}{2}``. The sums are the
rectangle rule of `_apply_∫dv!` and its moment stencils.
"""
function velocity_moments(dist::GridDistribution{DT}) where {DT}
    nx, nv = size(dist.values)
    ci = CartesianIndices((nx, nv))
    li = LinearIndices((nx, nv))
    f = vec(dist.values)
    hx, hv = step(dist.x), step(dist.v)

    n = zeros(DT, nx)
    j = zeros(DT, nx)
    ε = zeros(DT, nx)
    _apply_∫dv!(n, f, ci, li, hx, hv)
    _apply_∫vdv!(j, f, dist.v, ci, li, hx, hv)
    _apply_∫v²dv!(ε, f, dist.v, ci, li, hx, hv)
    ε ./= 2

    return (density = n, momentum = j, energy = ε)
end
