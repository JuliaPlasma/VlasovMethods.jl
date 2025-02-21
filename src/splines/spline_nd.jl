
remap_unit_interval(y, x₀, x₁) = x₀ + y * (x₁ - x₀)

unique_knots(basis::AbstractBSplineBasis) = sort([Set(BSplineKit.knots(basis))...])
unique_knots(basis::PeriodicBSplineBasis) = sort([Set(BSplineKit.knots(basis))...])[begin+div(BSplineKit.order(basis),2):end-div(BSplineKit.order(basis)-1,2)]


function mass_matrix_integrand(basis::AbstractBSplineBasis, knots, i, j, k)
    remap = y -> remap_unit_interval(y, knots[k], knots[k+1])
    v -> (basis[i] ∘ remap)(v) * (basis[j] ∘ remap)(v)
end

function mass_matrix_integrand(basis::PeriodicBSplineBasis, knots, i, j, k)
    remap  = y -> remap_unit_interval(y, knots[k], knots[k+1])
    remapm = y -> remap_unit_interval(y, knots[k], knots[k+1]) - (knots[end] - knots[begin])
    remapp = y -> remap_unit_interval(y, knots[k], knots[k+1]) + (knots[end] - knots[begin])
    v -> ((basis[i] ∘ remap)(v) + (basis[i] ∘ remapm)(v) + (basis[i] ∘ remapp)(v)) *
         ((basis[j] ∘ remap)(v) + (basis[j] ∘ remapm)(v) + (basis[j] ∘ remapp)(v))
end

function mass_matrix(basis::AbstractBSplineBasis, quadrature::QuadratureRule)
    M = zeros(length(basis), length(basis))
    knots = unique_knots(basis)

    for i in axes(M,1)
        for j in axes(M,2)
            for k in eachindex(knots[begin:end-1])
                M[i,j] += quadrature(mass_matrix_integrand(basis, knots, i, j, k)) * (knots[k+1] - knots[k])
            end
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


struct SplineND{T, D, BT <: AbstractBSplineBasis, DT, CT1 <: AbstractArray{T}, CT2 <: AbstractVector{T}, MT <: AbstractMatrix{T}, FT, QT <: QuadratureRule{T}}
    basis::BT
    derivative::DT
    coefficients::CT1
    coeff_vector::CT2
    mass_matrix::MT
    mass_factor::FT
    quadrature::QT

    function SplineND{T,D}(basis::AbstractBSplineBasis, quadrature::QuadratureRule; mass_quadrature = quadrature) where {T,D}
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
                result += s.coefficients[i,j] * bi * bj
            end
        end
    end
    
    return result
end


(s::SplineND{T})(x::AbstractVector{T}) where {T} = evaluate(s, x)
(s::SplineND{T,D})(x::Vararg{T,D}) where {D, T <: Number} = evaluate(s, SVector{D}(x...))


function cartesian_indices(::SplineND{T,1}, i::Int) where {T}
    return (i,)
end

function cartesian_indices(s::SplineND{T,2}, inds::Int) where {T}
    L = length(basis(s))
    i = mod1(inds, L)
    j = div(inds - i, L) + 1
    return (i, j)
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
    evaluate_basis(s, x, cartesian_indices(s, i))
end


function _evaluate_basis_derivative(B::AbstractBSplineBasis{O,T}, x::AbstractVector{T}, I::Tuple, j::Int) where {O,T}
    result = one(T)
    for d in eachindex(x)
        result *= d == j ? B[I[d],T](x[d], BSplineKit.Derivative(1)) : B[I[d],T](x[d])
    end
    return result
end

function evaluate_basis_derivative(s::SplineND{T,D}, x::AbstractVector{T}, I::Tuple) where {T,D}
    @assert D == length(I) == length(x)
    SVector{D,T}((_evaluate_basis_derivative(basis(s), x, I, d) for d in 1:D)...)
end

function evaluate_basis_derivative(s::SplineND, x::AbstractVector, i::Int)
    evaluate_basis_derivative(s, x, cartesian_indices(s, i))
end
