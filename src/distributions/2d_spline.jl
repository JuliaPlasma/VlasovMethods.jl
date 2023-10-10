abstract type NDSpline{DT, BT} end

struct TwoDSpline{DT, BT, BT2} <: NDSpline{DT, BT}
    basis::BT
    coefficients::Matrix{DT}
    basis_der::BT2

    function TwoDSpline(basis::BT, coefficients::Matrix{DT}) where {BT, DT}
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
    k = order(S)
    
    i1, bs1 = evaluate_all(B, x[1])
    i2, bs2 = evaluate_all(B, x[2])

    eval_b1 = collect(bs1[end:-1:1])
    eval_b2 = collect(bs2[end:-1:1])
    coefs = S.coefficients[i1-k+1:i1,i2-k+1:i2]
    
    return transpose(eval_b1)*coefs*eval_b2
end


function _first_derivative(S::TwoDSpline, component::Int)
    coefs = S.coefficients
    t = knots(S)
    k = order(S)

    d_coefs = similar(coefs)
    component == 1 ? copy!(d_coefs, transpose(coefs)) : copy!(d_coefs, coefs) # make sure axes of d_coefs matches with what we assume in the for loop below. 
    # du_i = similar(coefs[1,:])
    
    # assumes the first axis of coefs is invariant, and the second is being differentiated.
    for i in axes(coefs, 1)
        # du_i .= coefs[i,:] # get vector in coefficient matrix corresponding to the i-th spline in the non-differentiated direction 
        for j in Iterators.Reverse(axes(coefs, 2)) # now loop over other direction (in which we are differentiating)
            dt = t[j + k - 1] - t[j] 
            if iszero(dt) || j == firstindex(coefs[i,:]) # check for zeros in denom. 
                d_coefs[i,j] = zero(eltype(d_coefs))
            else
                d_coefs[i,j] = (k - 1) * (d_coefs[i,j] - d_coefs[i,j-1]) / dt
            end
        end
    end

    d_coefs = d_coefs[:,2:end]

    component == 1 ? d_coefs = transpose(d_coefs) : nothing

    return d_coefs
end

function evaluate_all_derivative(B::AbstractBSplineBasis, v::Vector{T}) where T


end

function evaluate_first_derivative(S::TwoDSpline, v::Vector{DT}) where {DT}
    # construct first derivative basis
    B = basis(S)
    k = order(S)
    # B_der = BSplineKit.BSplines.basis_derivative(B, Derivative(1))
    B_der = S.basis_der

    i1, bs1 = evaluate_all(B, v[1])
    eval_b1 = collect(bs1[end:-1:1])
    i2, bs2 = evaluate_all(B, v[2])
    eval_b2 = collect(bs2[end:-1:1])

    # compute ∂S/∂x
    coefs_dx = _first_derivative(S, 1)
    j1, ds1 = evaluate_all(B_der, v[1])
    eval_ds1 = collect(ds1[end:-1:1])

    # if j1-k+2 == 1
    #     m1 = 2 
    #     eval_ds1 = eval_ds1[2:end]
    # else
    #     m1 = j1-k+2
    # end
    m1 = j1-k+2
    dSdx = transpose(eval_ds1)*coefs_dx[m1:j1,i2-k+1:i2]*eval_b2

    # compute ∂S/∂y
    coefs_dy = _first_derivative(S, 2)
    j2, ds2 = evaluate_all(B_der, v[2])
    eval_ds2 = collect(ds2[end:-1:1])

    # if j2-k+2 == 1
    #     m2 = 2 
    #     eval_ds2 = eval_ds2[2:end]
    # else
    #     m2 = j2-k+2
    # end
    m2 = j2-k+2
    dSdy = transpose(eval_b1)*coefs_dy[i1-k+1:i1,m2:j2]*eval_ds2

    # return evaluated gradient vector
    return dSdx, dSdy, coefs_dx
end


# TODO: Can we delete this??
# function DB_ik(B_d, t, i, k, v::T) where T
#     # @show T
#     x = zero(T)
#     local Bi::T
#     local Bj::T

#     if i != 1
#         # println("before eval")
#         Bi = B_d[i](v)
#         # println("eval1")
#         # println(Bi, typeof(Bi))
#         Bj = B_d[i-1](v)
#         # println(Bj, typeof(Bj))
#         # println("eval2")
#         x += ((k-1)*(-Bi/(t[i+k] - t[i+1]) + Bj/(t[i+k-1] - t[i])))
#         # x += ((k-1)*(-B_d[i](v)/(t[i+k] - t[i+1]) + B_d[i-1](v)/(t[i+k-1] - t[i])))
#     end
#     return x
# end

# function basis_function_derivatives(S::TwoDSpline, i::Int, j::Int)
#     B = basis(S)
#     B_d = S.basis_der
#     k = order(S)
#     t = knots(S)

#     f_i(v) = DB_ik(B_d, t, i, k, v)* B[j](v)
#     f_j(v) = B[i](v)*DB_ik(B_d, t, j, k, v)

#     return f_i, f_j
# end

# # evaluate (i,j)-th tensor product basis function derivatives at v
# function eval_bfd(S::TwoDSpline, i::Int, j::Int, v::Vector{T}) where {T}
#     f, g = basis_function_derivatives(S, i, j)
#     return [f(v[1]), g(v[2])]::Vector{T}
# end

function eval_bfd(B::AbstractBSplineBasis, i, j, v::Vector{T}) where T
    a = B[i,T]
    b = B[j,T]

    return [a(v[1], Derivative(1)) * b(v[2]), a(v[1]) * b(v[2], Derivative(1))]
end