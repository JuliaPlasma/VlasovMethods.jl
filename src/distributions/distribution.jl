
abstract type DistributionFunction{DT, XD, VD} end

xdim(::DistributionFunction{DT, XD, VD}) where {DT, XD, VD} = XD
vdim(::DistributionFunction{DT, XD, VD}) where {DT, XD, VD} = VD

function (d::DistributionFunction{T})(x::AbstractVector{T}, v::AbstractVector{T}) where {T}
    evaluate(d, x, v)
end
function (d::DistributionFunction{T, XD, 0})(x::Vararg{T, XD}) where {T, XD}
    evaluate(d, SVector{XD}(x...), missing)
end
function (d::DistributionFunction{T, 0, VD})(v::Vararg{T, VD}) where {T, VD}
    evaluate(d, missing, SVector{VD}(v...))
end

function (d::DistributionFunction{T, XD, VD})(z::Vararg{T, ZD}) where {T, XD, VD, ZD}
    @assert ZD == XD+VD
    d(SVector{XD}(ntuple(i -> z[i], XD)), SVector{VD}(ntuple(i -> z[XD + i], VD)))
end
