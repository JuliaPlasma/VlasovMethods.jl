
abstract type DistributionFunction{DT,XD,VD} end

(d::DistributionFunction{T})(x::AbstractVector{T}, v::AbstractVector{T}) where {T} = evaluate(d, x, v)
(d::DistributionFunction{T,XD,0})(x::Vararg{T,XD}) where {T,XD} = evaluate(d, SVector{XD}(x...), missing)
(d::DistributionFunction{T,0,VD})(v::Vararg{T,VD}) where {T,VD} = evaluate(d, missing, SVector{VD}(v...))

function (d::DistributionFunction{T,XD,VD})(z::Vararg{T,ZD}) where {T,XD,VD,ZD}
    @assert ZD == XD+VD
    x = @view z[1:XD]
    v = @view z[XD+1:XD+VD]
    d(SVector{XD}(x...), SVector{VD}(v...))
end
