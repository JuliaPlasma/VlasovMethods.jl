
remap_unit_interval(y, x₀, x₁) = x₀ + y * (x₁ - x₀)

function remap_unit_to_knot_interval(f::Callable, knots, k)
    v -> f(SA_F64[remap_unit_interval(v[1], knots[k], knots[k+1])])
end

function remap_unit_to_knot_interval(f::Callable, knots, k, l)
    v -> f(SA_F64[remap_unit_interval(v[1], knots[k], knots[k+1]), remap_unit_interval(v[2], knots[l], knots[l+1])])
end

function remap_basis_from_unit_to_knot_interval(i::Int, basis::AbstractBSplineBasis, knots, k)
    v -> basis[i](remap_unit_interval(v, knots[k], knots[k+1]))
end

function remap_basis_from_unit_to_knot_interval(i::Int, basis::PeriodicBSplineBasis, knots, k)
    v -> basis[i](remap_unit_interval(v, knots[k], knots[k+1])) +
         basis[i](remap_unit_interval(v, knots[k], knots[k+1]) - (knots[end] - knots[begin])) +
         basis[i](remap_unit_interval(v, knots[k], knots[k+1]) + (knots[end] - knots[begin]))
end

unique_knots(basis::AbstractBSplineBasis) = sort([Set(BSplineKit.knots(basis))...])
unique_knots(basis::PeriodicBSplineBasis) = sort([Set(BSplineKit.knots(basis))...])[begin+div(BSplineKit.order(basis), 2):end-div(BSplineKit.order(basis) - 1, 2)]

function mass_matrix_quadrature(basis::AbstractBSplineBasis, knots, quadrature::QuadratureRule, i, j)
    m = 0.0
    for k in eachindex(knots[begin:end-1])
        integrand = v -> remap_basis_from_unit_to_knot_interval(i, basis, knots, k)(v) *
                         remap_basis_from_unit_to_knot_interval(j, basis, knots, k)(v)
        m += quadrature(integrand) * (knots[k+1] - knots[k])
    end
    return m
end

function mass_matrix(basis::AbstractBSplineBasis, quadrature::QuadratureRule)
    M = zeros(length(basis), length(basis))
    knots = unique_knots(basis)

    for i in axes(M, 1)
        for j in axes(M, 2)
            M[i, j] = mass_matrix_quadrature(basis, knots, quadrature, i, j)
        end
    end

    return M
end


spline_basis(p::BSplineOrder, knots::AbstractVector, ::Val{:Natural}) =
    BSplineBasis(p, knots)

spline_basis(p::BSplineOrder, knots::AbstractVector, ::Val{:Dirichlet}) =
    RecombinedBSplineBasis(BSplineKit.Derivative(0), BSplineBasis(p, knots))

spline_basis(p::BSplineOrder, knots::AbstractVector, ::Val{:Periodic}) =
    PeriodicBSplineBasis(p, knots)

function spline_basis_derivative(basis::AbstractBSplineBasis)
    if BSplineKit.BSplines.has_parent_basis(basis)
        return BSplineKit.BSplines.basis_derivative(parent(basis), Derivative(1))
    else
        return BSplineKit.BSplines.basis_derivative(basis, Derivative(1))
    end
end

function _mass_kron(D, mass_mat)
    if D == 1
        return mass_mat
    else
        kron((mass_mat for _ in 1:D)...)
    end
end


struct SplineND{T,D,BT<:AbstractBSplineBasis,DT,CT1<:AbstractArray{T},CT2<:AbstractVector{T},MT<:AbstractMatrix{T},FT,QT<:QuadratureRule{T}}
    basis::BT
    derivative::DT
    coefficients::CT1
    coeff_vector::CT2
    mass_matrix::MT
    mass_factor::FT
    quadrature::QT

    function SplineND{T,D}(basis::AbstractBSplineBasis, quadrature::QuadratureRule; mass_quadrature=quadrature) where {T,D}
        basis_der = spline_basis_derivative(basis)

        mass_mat = _mass_kron(D, mass_matrix(basis, mass_quadrature))
        mass_fac = cholesky(mass_mat)

        coefficients = zeros(T, (length(basis) for _ in 1:D)...)
        coeff_vector = vec(coefficients)

        new{
            T,
            D,
            typeof(basis),
            typeof(basis_der),
            typeof(coefficients),
            typeof(coeff_vector),
            typeof(mass_mat),
            typeof(mass_fac),
            typeof(quadrature)
        }(
            basis,
            basis_der,
            coefficients,
            coeff_vector,
            mass_mat,
            mass_fac,
            quadrature
        )
    end
end

SplineND(ndims::Int, args...; kwargs...) =
    SplineND{Float64,ndims}(args...; kwargs...)

SplineND(ndims::Int, o::Int, knots::AbstractVector, bcs::Symbol, args...; kwargs...) =
    SplineND(ndims, spline_basis(BSplineOrder(o), copy(knots), Val(bcs)), args...; kwargs...)

SplineND(ndims::Int, o::Int, knots::AbstractVector, args...; kwargs...) =
    SplineND(ndims, o, knots, :Natural, args...; kwargs...)


basis(s::SplineND) = s.basis
knots(s::SplineND) = BSplineKit.knots(basis(s))
unique_knots(s::SplineND) = unique_knots(basis(s))

order(::SplineND{T,D,BT}) where {T,D,BT} = BSplineKit.order(BT)

Base.eltype(::SplineND{T,D}) where {T,D} = T
Base.ndims(::SplineND{T,D}) where {T,D} = D
Base.size(s::SplineND{T,D}) where {T,D} = Tuple(length(basis(s)) for _ in 1:D)
Base.length(s::SplineND{T,D}) where {T,D} = *((length(basis(s)) for _ in 1:D)...)


map_index(::AbstractBSplineBasis, i) = i
map_index(basis::PeriodicBSplineBasis, i) = mod1(i, length(basis))

function evaluate(s::SplineND{T,1}, x::AbstractVector{T}) where {T}
    result::T = 0

    ilast, bi = BSplineKit.evaluate_all(basis(s), x[1])

    for (δi, bi) ∈ pairs(bi)
        i = map_index(basis(s), ilast - δi + 1)
        if i ≥ 1 && i ≤ length(basis(s))
            result += s.coefficients[i] * bi
        end
    end

    return result
end

function evaluate(s::SplineND{T,2}, x::AbstractVector{T}) where {T}
    result::T = 0

    ilast, bi = BSplineKit.evaluate_all(basis(s), x[1])
    jlast, bj = BSplineKit.evaluate_all(basis(s), x[2])

    for (δi, bi) ∈ pairs(bi)
        for (δj, bj) ∈ pairs(bj)
            i = map_index(basis(s), ilast - δi + 1)
            j = map_index(basis(s), jlast - δj + 1)
            if i ≥ 1 && i ≤ length(basis(s)) &&
               j ≥ 1 && j ≤ length(basis(s))
                result += s.coefficients[i, j] * bi * bj
            end
        end
    end

    return result
end


(s::SplineND{T})(x::AbstractVector{T}) where {T} = evaluate(s, x)
(s::SplineND{T,D})(x::Vararg{T,D}) where {D,T<:Number} = evaluate(s, SVector{D}(x...))


function cartesian_index(::SplineND{T,1}, i::Int) where {T}
    return (i,)
end

function cartesian_index(s::SplineND{T,2}, inds::Int) where {T}
    L = length(basis(s))
    i = mod1(inds, L)
    j = div(inds - i, L) + 1
    return (i, j)
end

function linear_index(::SplineND{T,1}, i::Int) where {T}
    return i
end

function linear_index(s::SplineND{T,2}, i::Int, j::Int) where {T}
    return (j - 1) * length(basis(s)) + i
end

function linear_indices(::SplineND{T,1}, i::AbstractVector{Int}) where {T}
    SVector{length(i)}(i)
end

function linear_indices(s::SplineND{T,2}, i::AbstractVector{Int}, j::AbstractVector{Int}) where {T}
    i_stat = SVector{length(i)}(i)
    j_stat = SVector{length(j)}(j)
    i_ones = @SVector ones(Int, length(i))
    j_ones = @SVector ones(Int, length(j))

    index_list = (i_ones * j_stat' .- 1) .* length(basis(s)) .+ i_stat * j_ones'

    return vec(index_list)
end


function _evaluate_basis(B::AbstractBSplineBasis{O,T}, x::AbstractVector{T}, I::Tuple) where {O,T}
    result = one(T)

    for d in eachindex(x)
        result *= B[I[d]](x[d])
    end

    return result
end

function evaluate_basis(s::SplineND, x::AbstractVector, I::Tuple)
    @assert ndims(s) == length(I) == length(x)
    _evaluate_basis(basis(s), x, I)
end

function evaluate_basis(s::SplineND, x::AbstractVector, i::Int)
    evaluate_basis(s, x, cartesian_index(s, i))
end


function _evaluate_basis_derivative(B::AbstractBSplineBasis{O,T}, x::AbstractVector{T}, I::Tuple, j::Int) where {O,T}
    result = one(T)
    for d in eachindex(x)
        result *= d == j ? B[I[d], T](x[d], BSplineKit.Derivative(1)) : B[I[d], T](x[d])
    end
    return result
end

function evaluate_basis_derivative(s::SplineND{T,D}, x::AbstractVector{T}, I::Tuple) where {T,D}
    @assert D == length(I) == length(x)
    SVector{D,T}((_evaluate_basis_derivative(basis(s), x, I, d) for d in 1:D)...)
end

function evaluate_basis_derivative(s::SplineND, x::AbstractVector, i::Int)
    evaluate_basis_derivative(s, x, cartesian_index(s, i))
end


function evaluate_indices(s::SplineND{T,D}, x::AbstractVector{T}) where {T,D}
    offset = BSplineKit.Recombinations.num_constraints(basis(s))[1]
    lastind = (BSplineKit.find_knot_interval(knots(s), x[i])[1] - offset for i in 1:D)
    indices = (map_index.(basis(s), l .+ 1 .- SVector{order(s)}(1:order(s))) for l in lastind)
    linear_indices(s, indices...)
end

evaluate_indices(s::SplineND{T,D}, x::Vararg{T,D}) where {T,D} = evaluate_indices(s, SVector{D}(x...))


function L2product_quadrature(f::Callable, basis::AbstractBSplineBasis, knots, quadrature::QuadratureRule, i)
    m = 0.0

    for k in eachindex(knots[begin:end-1])
        g = v -> remap_basis_from_unit_to_knot_interval(i, basis, knots, k)(v) *
                 remap_unit_to_knot_interval(f, knots, k)(SA_F64[v])
        m += quadrature(g) * (knots[k+1] - knots[k])
    end

    return m
end

function L2product_quadrature(f::Callable, basis::AbstractBSplineBasis, knots, quadrature::QuadratureRule, i, j)
    m = 0.0

    for k in eachindex(knots[begin:end-1])
        for l in eachindex(knots[begin:end-1])
            integrand = (u, v) -> remap_basis_from_unit_to_knot_interval(i, basis, knots, k)(u) *
                                  remap_basis_from_unit_to_knot_interval(j, basis, knots, l)(v) *
                                  remap_unit_to_knot_interval(f, knots, k, l)(SA_F64[u, v])

            quad = (i, j) -> quadrature.weights[i] * quadrature.weights[j] * integrand(quadrature.nodes[i], quadrature.nodes[j])

            interval = (knots[k+1] - knots[k]) * (knots[l+1] - knots[l])

            m += mapreduce(ind -> quad(ind...), +, zip(eachindex(quadrature), eachindex(quadrature))) * interval
        end
    end

    return m
end

function L2projection!(f::Callable, s::SplineND{T}) where {T}
    L = [L2product_quadrature(f, basis(s), unique_knots(s), s.quadrature, cartesian_index(s, l)...) for l in 1:length(s)]
    ldiv!(s.coeff_vector, s.mass_factor, L)
    return s
end
