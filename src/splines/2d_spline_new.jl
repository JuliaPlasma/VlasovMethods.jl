struct TwoDSpline{DT, BT, BT2} <: NDSpline{DT, BT}
    basis::BT
    coefficients::Vector{DT}
    basis_der::BT2
    mass_matrix::Matrix{DT}

    function TwoDSpline(basis::BT, coefficients::Vector{DT}) where {BT, DT}
        if BSplineKit.BSplines.has_parent_basis(basis)
            basis_der = BSplineKit.BSplines.basis_derivative(parent(basis), Derivative(1))
        else
            basis_der = BSplineKit.BSplines.basis_derivative(basis, Derivative(1))
        end
        new{DT, typeof(basis), typeof(basis_der)}(basis, coefficients, basis_der)
    end
end

# function TwoDSpline(basis::BT, coefficients::Matrix{DT}) where {BT, DT}
#     new{DT, typeof(basis), typeof(basis)}(basis, basis, coefficients)
# end


knots(S::TwoDSpline) = BSplineKit.knots(S.basis)
basis(S::TwoDSpline) = S.basis
# order(S::TwoDSpline) = order(S.basis)

order(::Type{<:TwoDSpline{T,Basis}}) where {T,Basis} = BSplineKit.order(Basis)
order(S::TwoDSpline) = order(typeof(S))

(S::TwoDSpline)(x::Float64,y::Float64) = (S::TwoDSpline)([x, y])
(S::TwoDSpline)(x) = evaluate(S, x)

function evaluate(S::TwoDSpline, x::AbstractVector)
    B = basis(S)
    M = length(B)
    k = order(S)
    result = zero(eltype(S.coefficients))
    
    ilast1, bs1 = evaluate_all(B, x[1])
    ilast2, bs2 = evaluate_all(B, x[2])

    for (δi, bi) ∈ pairs(bs1)
        for (δi2, b2) ∈ pairs(bs2)
            i = ilast1 + 1 - δi
            j = ilast2 + 1 - δi2
            if i > 0 && i <= M && j > 0 && j <= M
                result += S.coefficients[(j - 1)* M + i] * bi * b2
            end
        end
    end
    return result
end

function eval_2d_basis_func(B::AbstractBSplineBasis, x::AbstractVector, index::Int)
    i,j = ij_from_k(index, length(B))

    return B[i,T](x[1]) * B[j,T](x[2])
end

function ij_from_k(k::Int, M::Int)
    i = mod1(k, M)
    j = Int((k - i)/M) + 1

    return i, j
end

function eval_bfd(B::AbstractBSplineBasis, i, j, v::Vector{T}) where T
    a = B[i,T]
    b = B[j,T]

    return [a(v[1], Derivative(1)) * b(v[2]), a(v[1]) * b(v[2], Derivative(1))]
end

function eval_bfd(B::AbstractBSplineBasis, k, v::Vector{T}) where T
    # M = length(B)
    i,j = ij_from_k(k, length(B))

    a = B[i,T]
    b = B[j,T]

    return [a(v[1], Derivative(1)) * b(v[2]), a(v[1]) * b(v[2], Derivative(1))]
end

function evaluate_der_2d(B::AbstractBSplineBasis, v::Vector{T}) where T
    M = length(B)

    i1, bs1 = evaluate_all(B, v[1])
    i2, bs2 = evaluate_all(B, v[2])

    #derivatives 
    i1_der, bs1_der = evaluate_all(B, v[1], Derivative(1))
    i2_der, bs2_der = evaluate_all(B, v[2], Derivative(1))

    index_list = zeros(Int, length(bs1) * length(bs2))
    result = zeros(T,(2,length(bs1) * length(bs2)))

    count = 1 #TODO: should make this indexing cleaner
    for (δi, bi) ∈ pairs(bs1)
        for (δi2, b2) ∈ pairs(bs2)

            i = i1 + 1 - δi
            j = i2 + 1 - δi2
            k = (j - 1)*M + i

            index_list[count] = Int(k)
            result[:, count] .= [bs1_der[δi] * b2, bi * bs2_der[δi2]]

            count += 1
        end
    end

    return index_list, result
end