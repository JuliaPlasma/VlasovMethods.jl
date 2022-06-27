module BumpOnTail

using Parameters
using ..ParticleMethods: draw_g_accept_reject, draw_g_importance_sampling

# x-part of distribution function
function fₓ(x::T, params::NamedTuple) where {T}
    @unpack κ, ε = params
    @. ( 1 - ε * cos(κ*x) )
end

# v-part of distribution function
function fᵥ(v::T, params::NamedTuple) where {T}
    @unpack a, v₀, σ = params
    @. ( (1-a) * exp( - v^2 / 2 ) + a / σ * exp( - (v-v₀)^2 / (2σ^2) ) ) / √(2π)
end

# distribution function
function f(x::T, v::T, params) where {T}
    fₓ(x, params) .* fᵥ(v, params)
end

draw_accept_reject(Nₚ, params) = draw_g_accept_reject(Nₚ, fₓ, params)
draw_importance_sampling(Nₚ, params) = draw_g_importance_sampling(Nₚ, fₓ, params)

end
