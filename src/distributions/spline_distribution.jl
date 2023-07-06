struct SplineDistribution{XD, VD, DT, BT, MT, FT} <: DistributionFunction{XD,VD}
    spline::Spline{DT, BT, Vector{DT}}
    basis::BT
    coefficients::Vector{DT}
    mass_matrix::MT
    mass_fact::FT

    function SplineDistribution(xdim, vdim, basis::BT, coefficients::Vector{DT}) where {BT, DT}
        spline = Spline(basis, coefficients)
        mass_matrix = galerkin_matrix(basis)
        mass_fact = cholesky(mass_matrix)
        new{xdim, vdim, DT, typeof(basis), typeof(mass_matrix), typeof(mass_fact)}(spline, basis, coefficients, mass_matrix, mass_fact)
    end
end

Base.size(dist::SplineDistribution) = size(dist.coefficients)
Base.eltype(::SplineDistribution{XD, VD, DT}) where {XD, VD, DT} = DT

Cache(ST, s::SplineDistribution{XD, VD, DT, BT, MT, FT}) where {XD, VD, DT, BT, MT, FT} = SplineDistribution(XD, VD, s.basis, zeros(ST, axes(s.coefficients)))
CacheType(ST, ::SplineDistribution{XD, VD, DT, BT, MT, FT}) where {XD, VD, DT, BT, MT, FT} = SplineDistribution{XD, VD, ST, BT, MT, FT}


function SplineDistribution(xdim, vdim, nknots::KT, order::OT, domain::Tuple, bc::Symbol=:Dirichlet) where {KT, OT}
    knots = collect(LinRange(domain..., nknots))
    b = BSplineBasis(BSplineOrder(order), knots)

    if bc == :Dirichlet
        basis = RecombinedBSplineBasis(Derivative(0), b)
    else #TODO: add other boundary condition options here
        basis = b
    end

    coefficients = zeros(length(basis))

    return SplineDistribution(xdim, vdim, basis, coefficients)
end


_cachehash(ST) = hash(Threads.threadid(), hash(ST))

struct SplineDistributionCache{disttype}
    s::disttype
    caches::Dict{UInt64, SplineDistribution}

    function SplineDistributionCache(s::SplineDistribution)
        caches = Dict{UInt64, SplineDistribution}()
        caches[_cachehash(eltype(s))] = s
        new{typeof(s)}(s, caches)
    end
end

@inline function Base.getindex(c::SplineDistributionCache, ST::DataType)
    key = _cachehash(ST)
    if haskey(c.caches, key)
        c.caches[key]
    else
        c.caches[key] = Cache(ST, c.s)
    end::CacheType(ST, c.s)
end
