# The reference of the HNN trainer: its algorithm in plain Julia, with nested ForwardDiff for both
# of its derivatives. It shares no code with `ext/VlasovMethodsHNNExt.jl`, which writes the
# gradient out by hand: here the network is three matrix products over weights laid out
# `(out, in)`, its input gradient and its parameter gradient are both forward mode, and the step is
# `p .-= η .* dp`. `test/hnn.jl` includes it in a module of its own.

using ForwardDiff

"""
    Network

The network of the reference as the five arrays of its parameter set: `W1` and `W2` of size
`(m, m)` with their biases `b1` and `b2`, and `W3` of size `(1, m)` with no bias.
"""
struct Network{T}
    W1::Matrix{T}
    b1::Vector{T}
    W2::Matrix{T}
    b2::Vector{T}
    W3::Matrix{T}
end

"""
    network(p, m)

The network of the flat parameter vector `p`, `[vec(W1); b1; vec(W2); b2; vec(W3)]`.
"""
function network(p::AbstractVector, m::Integer)
    k = m * m
    W1 = reshape(p[1:k], m, m)
    b1 = p[k .+ (1:m)]
    W2 = reshape(p[k .+ m .+ (1:k)], m, m)
    b2 = p[2k .+ m .+ (1:m)]
    W3 = reshape(p[2k .+ 2m .+ (1:m)], 1, m)
    return Network(W1, b1, W2, b2, W3)
end

"""
    potential(net, x)

`H(x) = sum(model(x))`: the two hidden layers apply `tanh` to `W x + b`, the output layer is linear
with no bias, and its single value is the sum of the output.
"""
function potential(net::Network, x::AbstractVector)
    a1 = tanh.(net.W1 * x .+ net.b1)
    a2 = tanh.(net.W2 * a1 .+ net.b2)
    return (net.W3 * a2)[1]
end

"""
    input_gradient(net, x)

The gradient of `potential` with respect to its input, by forward mode. The closure is built at
every call, so that this differentiation carries a tag of its own: with parameters that are duals
themselves, this is the inner half of the nested gradient, and a tag shared with the outer level
would mix the two perturbations.
"""
function input_gradient(net::Network, x::AbstractVector)
    ForwardDiff.gradient(xi -> potential(net, xi), x)
end

"""
    data_loss(net, data, target)

`Σᵢ ‖∇ₓH(xᵢ) − yᵢ‖²` over the columns of `data` and `target`.
"""
function data_loss(net::Network{T}, data::AbstractMatrix, target::AbstractMatrix) where {T}
    s = zero(T)
    for j in axes(data, 2)
        g = input_gradient(net, @view data[:, j])
        for i in eachindex(g)
            d = g[i] - target[i, j]
            s += d * d
        end
    end
    return s
end

"""
    parameter_gradient(p, data, target)

`∇ₚ Σᵢ ‖∇ₓH(xᵢ) − yᵢ‖²` at the flat parameters `p`, by forward mode over the `2m² + 3m`
parameters.
"""
function parameter_gradient(p::AbstractVector, data::AbstractMatrix, target::AbstractMatrix)
    ForwardDiff.gradient(pv -> data_loss(network(pv, size(data, 1)), data, target), p)
end

"""
    train(data, target, p, η, batches) -> (p, losses)

The training run: for each epoch `j`, one step `p .-= η .* dp` of the parameter gradient of the loss
over the samples of `batches[j]`, and then `losses[j]`, the loss over all the samples of `data` and
`target` after the step.
"""
function train(data::AbstractMatrix, target::AbstractMatrix, p, η, batches)
    m = size(data, 1)
    losses = zeros(eltype(p), length(batches))
    for j in eachindex(batches)
        idx = batches[j]
        p = p .- η .* parameter_gradient(p, data[:, idx], target[:, idx])
        losses[j] = data_loss(network(p, m), data, target)
    end
    return p, losses
end
