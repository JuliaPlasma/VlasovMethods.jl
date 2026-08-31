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

order(::Type{<:TwoDSpline{T, Basis}}) where {T, Basis} = BSplineKit.order(Basis)
order(S::TwoDSpline) = order(typeof(S))

# (S::TwoDSpline)(x, y) = (S::TwoDSpline)([x, y])
(S::TwoDSpline)(x) = evaluate(S, x)
(S::TwoDSpline)(x, y) = evaluate(S, [x, y])
# (S::TwoDSpline)(x, y) = evaluate(S, x, y)

function evaluate(S::TwoDSpline{T}, x::AbstractVector) where {T}
    B = basis(S)
    M = length(B)
    # k = order(S)
    result::T = 0

    ilast1, bs1 = BSplineKit.evaluate_all(B, x[1])
    ilast2, bs2 = BSplineKit.evaluate_all(B, x[2])

    for (δi, bi) in pairs(bs1)
        for (δi2, b2) in pairs(bs2)
            if typeof(B) <: PeriodicBSplineBasis
                i = mod1(ilast1 + 1 - δi, M)
                j = mod1(ilast2 + 1 - δi2, M)
                result += S.coefficients[(j - 1) * M + i] * bi * b2
            else
                i = ilast1 + 1 - δi
                j = ilast2 + 1 - δi2
                if i > 0 && i <= M && j > 0 && j <= M
                    result += S.coefficients[(j - 1) * M + i] * bi * b2
                end
            end
        end
    end

    return result
end

# function evaluate(S::TwoDSpline, x1::T, x2::T) where {T}
#     B = basis(S)
#     M = length(B)
#     k = order(S)
#     result = zero(eltype(S.coefficients))

#     ilast1, bs1 = evaluate_all(B, x1)
#     ilast2, bs2 = evaluate_all(B, x2)

#     for (δi, bi) ∈ pairs(bs1)
#         for (δi2, b2) ∈ pairs(bs2)
#             i = ilast1 + 1 - δi
#             j = ilast2 + 1 - δi2
#             if i > 0 && i <= M && j > 0 && j <= M
#                 result += S.coefficients[(j - 1)* M + i] * bi * b2
#             end
#         end
#     end
#     return result
# end

function eval_2d_basis_func(B::AbstractBSplineBasis, x::AbstractVector, index::Int)
    i, j = ij_from_k(index, length(B))

    return B[i, T](x[1]) * B[j, T](x[2])
end

function ij_from_k(k::Int, M::Int)
    i = mod1(k, M)
    j = div(k - i, M) + 1
    return i, j
end

# function eval_bfd(B::AbstractBSplineBasis, i, j, v::Vector{T}) where T
#     a = B[i,T]
#     b = B[j,T]

#     return [a(v[1], Derivative(1)) * b(v[2]), a(v[1]) * b(v[2], Derivative(1))]
# end

function eval_bfd(B::AbstractBSplineBasis, k, v::AbstractVector{T}) where {T}
    # M = length(B)
    i, j = ij_from_k(k, length(B))

    a = B[i, T]
    b = B[j, T]

    @SVector [a(v[1], BSplineKit.Derivative(1)) * b(v[2]),
        a(v[1]) * b(v[2], BSplineKit.Derivative(1))]
end

function eval_bfd(B::AbstractBSplineBasis, k, v1::T, v2::T) where {T}
    i, j = ij_from_k(k, length(B))

    a = B[i, T]
    b = B[j, T]

    @SVector [
        a(v1, BSplineKit.Derivative(1)) * b(v2), a(v1) * b(v2, BSplineKit.Derivative(1))]
end

# function eval_bfd!(derivative::AbstractVector{T}, B::AbstractBSplineBasis, k, v, α, β) where T
#     # M = length(B)
#     i, j = ij_from_k(k, length(B))

#     a = B[i,T]
#     b = B[j,T]

#     derivative[1] = α * derivative[1] + β * a(v[1], Derivative(1)) * b(v[2])
#     derivative[2] = α * derivative[2] + β * a(v[1]) * b(v[2], Derivative(1))
# end

# function eval_bfd!(derivative::AbstractVector{T}, B::AbstractBSplineBasis, k, v1::T, v2::T, α, β) where T
#     # M = length(B)
#     i, j = ij_from_k(k, length(B))

#     a = B[i,T]
#     b = B[j,T]

#     derivative[1] = α * derivative[1] + β * a(v1, Derivative(1)) * b(v2)
#     derivative[2] = α * derivative[2] + β * a(v1) * b(v2, Derivative(1))
# end

function evaluate_der_2d(B::AbstractBSplineBasis, v::AbstractVector{T}) where {T}
    M = length(B)

    i1, bs1 = BSplineKit.evaluate_all(B, v[1])
    i2, bs2 = BSplineKit.evaluate_all(B, v[2])

    #derivatives 
    i1_der, bs1_der = BSplineKit.evaluate_all(B, v[1], BSplineKit.Derivative(1))
    i2_der, bs2_der = BSplineKit.evaluate_all(B, v[2], BSplineKit.Derivative(1))

    index_list = zeros(Int, length(bs1) * length(bs2))
    result = zeros(T, (2, length(bs1) * length(bs2)))

    count = 1 #TODO: should make this indexing cleaner
    for (δi, bi) in pairs(bs1)
        for (δi2, b2) in pairs(bs2)
            i = i1 + 1 - δi
            j = i2 + 1 - δi2
            k = (j - 1)*M + i

            index_list[count] = Int(k)

            result[1, count] = bs1_der[δi] * b2
            result[2, count] = bi * bs2_der[δi2]

            count += 1
        end
    end

    return index_list, result
end

function evaluate_der_2d_indices(B::AbstractBSplineBasis, v::AbstractVector{T}) where {T}
    # performes the following calculation
    # for δi1 ∈ eachindex(bs1)
    #     for δi2 ∈ eachindex(bs2)
    #         i = i1 + 1 - δi1
    #         j = i2 + 1 - δi2
    #         k = (j - 1) * length(B) + i
    #         index_list[count] = k
    #         count += 1
    #     end
    # end

    i1, bs1 = BSplineKit.evaluate_all(B, v[1])
    i2, bs2 = BSplineKit.evaluate_all(B, v[2])

    i = i1 .+ 1 .- SVector{length(bs1)}(eachindex(bs1))
    j = i2 .+ 1 .- SVector{length(bs2)}(eachindex(bs2))

    i_ones = @SVector ones(Int, length(i))
    j_ones = @SVector ones(Int, length(j))

    index_list = (i_ones * j' .- 1) .* length(B) .+ i * j_ones'

    return vec(index_list)
end

function evaluate_der_2d_indices(B::AbstractBSplineBasis, v1::T, v2::T) where {T}
    M = length(B)

    i1, bs1 = BSplineKit.evaluate_all(B, v1)
    i2, bs2 = BSplineKit.evaluate_all(B, v2)

    #derivatives 
    # i1_der, bs1_der = evaluate_all(B, v[1], Derivative(1))
    # i2_der, bs2_der = evaluate_all(B, v[2], Derivative(1))

    index_list = zeros(Int, length(bs1) * length(bs2))
    # result = zeros(T,(2,length(bs1) * length(bs2)))

    count = 1 #TODO: should make this indexing cleaner
    for (δi, bi) in pairs(bs1)
        for (δi2, b2) in pairs(bs2)
            i = i1 + 1 - δi
            j = i2 + 1 - δi2
            k = (j - 1)*M + i

            index_list[count] = Int(k)

            # result[1, count] = bs1_der[δi] * b2
            # result[2, count] = bi * bs2_der[δi2]

            count += 1
        end
    end

    return index_list
end

function shift_scale(i, knots)
    # @show 0.5 * (knots[i] + knots[i+1])
    # @show 0.5 * (knots[i + 1] - knots[i])
    return (knots[i] + knots[i + 1]) / 2, (knots[i + 1] - knots[i]) / 2
end

function linear_transform(x::T, shift::T, scale::T) where {T}
    return x .* scale .+ shift
end

function gauss_quad(f::Function, basis::AbstractBSplineBasis, n::Int, params)
    T = eltype(params.sdist)
    # k = BSplineKit.order(basis)
    # knots = BSplineKit.knots(basis)[k:end-k+1]
    knots = BSplineKit.knots(basis)[3:(end - 1)]
    M = length(knots) - 1
    x, w = gausslegendre(n)

    res::T = 0

    for i in 1:M, j in 1:M, l in 1:M, m in 1:M
        shift_i, scale_i = shift_scale(i, knots)
        shift_j, scale_j = shift_scale(j, knots)
        shift_l, scale_l = shift_scale(l, knots)
        shift_m, scale_m = shift_scale(m, knots)

        shift_ij = @SVector [shift_i, shift_j]
        shift_lm = @SVector [shift_l, shift_m]
        scale_ij = @SVector [scale_i, scale_j]
        scale_lm = @SVector [scale_l, scale_m]

        for i1 in 1:n, j1 in 1:n, l1 in 1:n, m1 in 1:n
            x_ij = @SVector [x[i1], x[j1]]
            x_lm = @SVector [x[l1], x[m1]]

            v1 = linear_transform(x_ij, shift_ij, scale_ij)
            v2 = linear_transform(x_lm, shift_lm, scale_lm)

            res += f(v1, v2, params) * w[i1] * w[j1] * w[l1] * w[m1] * scale_i * scale_j *
                   scale_l * scale_m
        end
    end

    return res
end

# function gauss_quad(f::Function, basis::AbstractBSplineBasis, iknots::AbstractRange, jknots::AbstractRange, n::Int, params)
#     T = eltype(params.sdist)
#     # M = length(basis) + 1
#     # k = BSplineKit.order(basis)
#     knots = BSplineKit.knots(basis)
#     M = length(knots) - 1
#     x, w = gausslegendre(n)

#     res::T = 0

#     for i in iknots[1:end-1], j in jknots[1:end-1], l in iknots[1:end-1], m in jknots[1:end-1]
#         shift_i, scale_i = shift_scale(i, knots)
#         shift_j, scale_j = shift_scale(j, knots)
#         shift_l, scale_l = shift_scale(l, knots)
#         shift_m, scale_m = shift_scale(m, knots)

#         for i1 in 1:n, j1 in 1:n, l1 in 1:n, m1 in 1:n
#             v1 = convert(Array{T}, linear_transform([x[i1], x[j1]], [shift_i, shift_j], [scale_i, scale_j]))
#             v2 = convert(Array{T}, linear_transform([x[l1], x[m1]], [shift_l, shift_m], [scale_l, scale_m]))
#             res += f(v1, v2, params) * w[i1] * w[j1] * w[l1] * w[m1] * scale_i * scale_j * scale_l * scale_m
#         end
#     end

#     return res
# end

function gauss_quad_2d(f::Function, basis::AbstractBSplineBasis, n::Int, params)
    T = eltype(params.sdist)
    # M = length(basis) + 1
    # k = BSplineKit.order(basis)
    # knots = BSplineKit.knots(basis)[k:end-k+1]
    knots = BSplineKit.knots(basis)[3:(end - 1)]
    M = length(knots) - 1
    x, w = gausslegendre(n)

    res::T = 0

    for i in 1:M, j in 1:M

        shift_i, scale_i = shift_scale(i, knots)
        shift_j, scale_j = shift_scale(j, knots)

        shift_ij = @SVector [shift_i, shift_j]
        scale_ij = @SVector [scale_i, scale_j]

        for i1 in 1:n, j1 in 1:n

            x_ij = @SVector [x[i1], x[j1]]

            v1 = linear_transform(x_ij, shift_ij, scale_ij)

            res += f(v1, params) * w[i1] * w[j1] * scale_i * scale_j
        end
    end

    return res
end

function gauss_quad_1d(f::Function, basis::AbstractBSplineBasis, n::Int, params)
    T = eltype(params.sdist)
    # M = length(basis) + 1
    # k = BSplineKit.order(basis)
    knots = BSplineKit.knots(basis)[1:end]
    M = length(knots) - 1
    # knots = BSplineKit.knots(basis)[k:end-k+1]
    x, w = gausslegendre(n)

    res::T = 0

    for i in 1:M
        shift_i, scale_i = shift_scale(i, knots)
        for i1 in 1:n
            v1 = convert(Array{T}, linear_transform(x[i1], shift_i, scale_i))
            res += f(v1, params) * w[i1] * scale_i
        end
    end

    return res
end
