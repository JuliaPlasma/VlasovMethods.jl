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

(S::TwoDSpline)(x,y) = (S::TwoDSpline)([x, y])
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
            result += S.coefficients[(j - 1)* M + i] * bi * b2
        end
    end
    return result
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
    M = length(B)
    i = mod(k, M)
    j = Int((k - i)/M)

    a = B[i,T]
    b = B[j,T]

    return [a(v[1], Derivative(1)) * b(v[2]), a(v[1]) * b(v[2], Derivative(1))]
end