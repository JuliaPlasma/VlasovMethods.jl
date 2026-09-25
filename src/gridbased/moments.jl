### integrate over velocity space
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
        i, j = Tuple(ci[ij])
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
