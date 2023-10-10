struct SplineDistribution{XD, VD, ST, DT, BT, MT, FT} <: DistributionFunction{XD,VD}
    spline::ST
    # spline::Spline{DT, BT, Vector{DT}}
    basis::BT
    coefficients::AbstractArray{DT}
    # coefficients::Vector{DT}
    mass_matrix::MT
    mass_fact::FT

    function SplineDistribution(xdim, vdim, basis::BT, coefficients::AbstractArray{DT}) where {BT, DT}
        if vdim == 1
            spline = Spline(basis, coefficients)
        elseif vdim == 2
            spline = TwoDSpline(basis, coefficients)
        end
        mass_matrix = galerkin_matrix(basis)
        mass_fact = cholesky(mass_matrix)
        new{xdim, vdim, typeof(spline), DT, typeof(basis), typeof(mass_matrix), typeof(mass_fact)}(spline, basis, coefficients, mass_matrix, mass_fact)
    end
end

Base.size(dist::SplineDistribution) = size(dist.coefficients)
Base.eltype(::SplineDistribution{XD, VD, ST, DT, BT, MT, FT}) where {XD, VD, ST, DT, BT, MT, FT} = DT

Cache(AT, s::SplineDistribution{XD, VD, ST, DT, BT, MT, FT}) where {XD, VD, ST, DT, BT, MT, FT} = SplineDistribution(XD, VD, s.basis, zeros(AT, axes(s.coefficients)))
CacheType(AT, ::SplineDistribution{XD, VD, TwoDSpline{DT, BT, BT2}, DT, BT, MT, FT}) where {XD, VD, ST, DT, BT, MT, FT, BT2} = SplineDistribution{XD, VD, TwoDSpline{AT, BT, BT2}, AT, BT, MT, FT}
CacheType(AT, ::SplineDistribution{XD, VD, Spline{DT, BT, Vector{DT}}, DT, BT, MT, FT}) where {XD, VD, ST, DT, BT, MT, FT} = SplineDistribution{XD, VD, Spline{AT, BT, Vector{AT}}, AT, BT, MT, FT}


function SplineDistribution(xdim, vdim, nknots::KT, s_order::OT, domain::Tuple, bc::Symbol=:Dirichlet) where {KT, OT}
    knots = collect(LinRange(domain..., nknots))
    b = BSplineBasis(BSplineOrder(s_order), knots)

    if bc == :Dirichlet
        basis = RecombinedBSplineBasis(Derivative(0), b)
    else #TODO: add other boundary condition options here
        basis = b
    end
    if vdim == 1
        coefficients = zeros(length(basis))
    elseif vdim == 2
        @show basis
        @show length(basis)
        coefficients = zeros(length(basis), length(basis))
    end
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
