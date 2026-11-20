### integrate over velocity space
@doc raw"""
    _apply_∫dv!(y, f, ci, li, hx, hv)

Integrate the grid values `f` over the velocity by the rectangle rule, and write the density
``y_i = \sum_j f_{ij} h_v`` at each ``x``-node into `y`.

`f` is the state vector `vec(dist.values)` of a [`GridDistribution`](@ref) of `nx × nv` nodes,
`ci` and `li` are the `CartesianIndices` and the `LinearIndices` of `(nx, nv)`, and `hx` and `hv`
are the grid widths; `hx` is not read. `y` holds `nx` values. A `y` or an `f` of another length
throws a `DimensionMismatch`. [`velocity_moments`](@ref) calls it, and with it the stencils of the
first and the second moment.
"""
function _apply_∫dv!(y::AbstractVector, f::AbstractVector, ci, li, h₁, h₂)
    nx, nv = size(ci)
    (length(f) == nx*nv && length(y) == nx) || throw(DimensionMismatch())
    y .= 0
    @inbounds for ij in li
        i, j = Tuple(ci[ij])
        # (∫ f dv)(xᵢ) = ∑ⱼ f(xᵢ,vⱼ) hᵥ
        y[i] += f[li[i, j]] * h₂
    end
end

function _apply_∫dvᵀ!(y::AbstractVector, ψ::AbstractVector, ci, li, h₁, h₂)
    nx, nv = size(ci)
    (length(ψ) == nx && length(y) == nx*nv) || throw(DimensionMismatch())
    y .= 0
    @inbounds for ij in li
        i, _ = Tuple(ci[ij])
        # (∫ f dv)(xᵢ) = ∑ⱼ f(xᵢ,vⱼ) hᵥ
        y[ij] = ψ[i] * h₂
    end
end

### weighted integral over velocity space
function _apply_∫dv_μ!(y::AbstractVector, f::AbstractVector, μ, ci, li, h₁, h₂)
    nx, nv = size(ci)
    (length(f) == nx*nv && length(y) == nx) || throw(DimensionMismatch())
    y .= 0
    @inbounds for ij in li
        i, j = Tuple(ci[ij])
        # (∫ f dv)(xᵢ) = ∑ⱼ f(xᵢ,vⱼ) hᵥ
        y[i] += f[li[i, j]] * h₂ * μ[j]
    end
end

function _apply_∫dvᵀ_μ!(y::AbstractVector, ψ::AbstractVector, μ, ci, li, h₁, h₂)
    nx, nv = size(ci)
    (length(ψ) == nx && length(y) == nx*nv) || throw(DimensionMismatch())
    y .= 0
    @inbounds for ij in li
        i, j = Tuple(ci[ij])
        # (∫ f dv)(xᵢ) = ∑ⱼ f(xᵢ,vⱼ) hᵥ
        y[ij] = ψ[i] * h₂ * μ[j]
    end
end

### First moment over velocity space
function _apply_∫vdv!(
        y::AbstractVector, f::AbstractVector, v::AbstractVector, ci, li, h₁, h₂)
    nx, nv = size(ci)
    (length(f) == nx*nv && length(y) == nx && length(v) == nv) || throw(DimensionMismatch())
    y .= 0
    @inbounds for ij in li
        i, j = Tuple(ci[ij])
        # (∫ f dv)(xᵢ) = ∑ⱼ f(xᵢ,vⱼ) hᵥ
        y[i] += f[li[i, j]] * h₂ * v[j]
    end
end

### Second moment over velocity space
function _apply_∫v²dv!(
        y::AbstractVector, f::AbstractVector, v::AbstractVector, ci, li, h₁, h₂)
    nx, nv = size(ci)
    (length(f) == nx*nv && length(y) == nx && length(v) == nv) || throw(DimensionMismatch())
    y .= 0
    @inbounds for ij in li
        i, j = Tuple(ci[ij])
        # (∫ f dv)(xᵢ) = ∑ⱼ f(xᵢ,vⱼ) hᵥ
        y[i] += f[li[i, j]] * h₂ * v[j]^2
    end
end
