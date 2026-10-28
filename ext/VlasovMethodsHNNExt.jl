module VlasovMethodsHNNExt

using AbstractNeuralNetworks: AbstractNeuralNetworks, Chain, Dense, IdentityActivation,
                              TanhActivation, input_dimension, output_dimension
using GeometricOptimizers: GradientMethod, NetworkParameters, TrainingOptimizer,
                           optimization_step!
using LinearAlgebra: mul!
using Random: rand!

import VlasovMethods: train_hnn

# The trainer of `src/hnn.jl`: the gradient of the loss is written out by hand, as reverse mode
# through the input gradient of the network, and the step comes from `GeometricOptimizers`. The two
# packages are weak dependencies, so this module loads only when both of them are loaded.

"The activations whose first and second derivatives the gradient below is written for."
const Identity = Union{
    IdentityActivation, AbstractNeuralNetworks.GenericActivation{typeof(identity)}}

"""
    Passes{T}

The buffers of the passes over one set of samples, one column per sample, for a network of `L`
layers. With `aₗ = σ(zₗ)`, `zₗ = Wₗ aₗ₋₁ + bₗ`, `sₗ = σ'(zₗ)` and the input gradient
`gₗ₋₁ = Wₗᵀ dₗ`, `dₗ = sₗ ⊙ gₗ`, from `g_L = 1`: `a` and `g` hold `L + 1` arrays, `a[1]` the
samples and `g[1]` the input gradient of the Hamiltonian, and `z`, `s` and `d` hold `L`. The
adjoints of the gradient pass are `ḡ`, of `L + 1` arrays like `g`, and `z̄` and `ā`, of `L` arrays
like `z` and the first `L` of `a`. `ā[1]`, the adjoint of the samples, keeps the indices of `a` and
is not written.
"""
struct Passes{T}
    z::Vector{Matrix{T}}
    a::Vector{Matrix{T}}
    s::Vector{Matrix{T}}
    g::Vector{Matrix{T}}
    d::Vector{Matrix{T}}
    ḡ::Vector{Matrix{T}}
    z̄::Vector{Matrix{T}}
    ā::Vector{Matrix{T}}
end

function Passes{T}(widths, samples) where {T}
    buffers(ws) = [Matrix{T}(undef, w, samples) for w in ws]
    layers = widths[2:end]
    return Passes{T}(buffers(layers), buffers(widths), buffers(layers), buffers(widths),
        buffers(layers), buffers(widths), buffers(layers), buffers(widths[1:(end - 1)]))
end

"""
    TrainingState

Everything an epoch reads and writes, built once for the run: the parameters and their gradient,
held as the weights and biases of each layer, the activation of each layer, the target of the whole
data set and of one batch, and the buffers of the passes over each. `W`, `b`, `dW` and `db` are the
arrays of `ps` and `dp`, not copies, so that a step of the optimizer on `ps` is the step the next
pass reads. A layer with no bias has an empty `b`.
"""
struct TrainingState{T, P}
    ps::P
    dp::P
    W::Vector{Matrix{T}}
    b::Vector{Vector{T}}
    dW::Vector{Matrix{T}}
    db::Vector{Vector{T}}
    tanh::Vector{Bool}
    target::Matrix{T}
    batch_target::Matrix{T}
    full::Passes{T}
    batch::Passes{T}
end

_bias(layer, ::Type{T}) where {T} = haskey(layer, :b) ? layer.b : T[]

"""
    TrainingState(model, ps, data, target, batchsize)

The state of a training run of the network `model` from the parameters `ps` over the columns of
`data` and `target`. The state holds `ps` itself and trains it in place; `data` and `target` are
copied into its buffers, which have one column per sample, and those of a batch one per sample of
the batch.
"""
function TrainingState(model::Chain, ps::NetworkParameters{T}, data::AbstractMatrix{T},
        target::AbstractMatrix{T}, batchsize::Integer) where {T}
    dp = deepcopy(ps)
    layers = keys(ps)
    W = Matrix{T}[ps[k].W for k in layers]
    widths = [size(W[1], 2); [size(w, 1) for w in W]]
    full = Passes{T}(widths, size(data, 2))
    copyto!(full.a[1], data)
    return TrainingState(ps, dp, W, Vector{T}[_bias(ps[k], T) for k in layers],
        Matrix{T}[dp[k].W for k in layers], Vector{T}[_bias(dp[k], T) for k in layers],
        Bool[layer.σ isa TanhActivation for layer in model.layers], copy(target),
        Matrix{T}(undef, size(target, 1), batchsize), full, Passes{T}(widths, batchsize))
end

"""
    loss!(passes, state, target)

`Σᵢ ‖∇ₓH(xᵢ) − yᵢ‖²` over the samples in `passes.a[1]` and the columns of `target`: the forward
pass, then the input gradient, both into `passes`.
"""
function loss!(p::Passes{T}, state::TrainingState{T}, target::AbstractMatrix{T}) where {T}
    L = length(state.W)
    for l in 1:L
        mul!(p.z[l], state.W[l], p.a[l])
        isempty(state.b[l]) || (p.z[l] .+= state.b[l])
        if state.tanh[l]
            p.a[l + 1] .= tanh.(p.z[l])
            p.s[l] .= 1 .- p.a[l + 1] .^ 2
        else
            p.a[l + 1] .= p.z[l]
            p.s[l] .= 1
        end
    end
    p.g[L + 1] .= 1
    for l in L:-1:1
        p.d[l] .= p.s[l] .* p.g[l + 1]
        mul!(p.g[l], state.W[l]', p.d[l])
    end
    loss = zero(T)
    for i in eachindex(target)
        loss += (p.g[1][i] - target[i])^2
    end
    return loss
end

"""
    gradient!(state)

The gradient of the loss over the batch with respect to the parameters, into `state.dp`. It is the
reverse of `loss!`: first back through the input gradient, from `ḡ₀ = 2(g₀ − y)`, where
`sₗ = 1 − aₗ²` of a `tanh` layer gives `z̄ₗ = −2 aₗ sₗ ⊙ gₗ ⊙ d̄ₗ`, and then back through the
forward pass, which adds `āₗ ⊙ sₗ` to `z̄ₗ`.
"""
function gradient!(state::TrainingState)
    p = state.batch
    L = length(state.W)
    loss!(p, state, state.batch_target)
    p.ḡ[1] .= 2 .* (p.g[1] .- state.batch_target)
    for l in 1:L
        # gₗ₋₁ = Wₗᵀ dₗ: the weights' share, and d̄ₗ in the buffer of z̄ₗ
        mul!(state.dW[l], p.d[l], p.ḡ[l]')
        mul!(p.z̄[l], state.W[l], p.ḡ[l])
        # dₗ = sₗ ⊙ gₗ, and sₗ = σ'(zₗ)
        p.ḡ[l + 1] .= p.z̄[l] .* p.s[l]
        if state.tanh[l]
            p.z̄[l] .*= -2 .* p.a[l + 1] .* p.s[l] .* p.g[l + 1]
        else
            p.z̄[l] .= 0
        end
    end
    # the loss reads the network's output through its input gradient only, so `ā_L = 0`
    for l in L:-1:1
        l < L && (p.z̄[l] .+= p.ā[l + 1] .* p.s[l])
        mul!(state.dW[l], p.z̄[l], p.a[l]', true, true)
        isempty(state.db[l]) || sum!(state.db[l], p.z̄[l])
        # the samples are not parameters, so the first layer passes no adjoint back
        l > 1 && mul!(p.ā[l], state.W[l]', p.z̄[l])
    end
    return state.dp
end

"""
    _fill_batch!(state, indices)

Copy the samples `indices` of the data and the target into the batch buffers of `state`, one
sample per column.
"""
function _fill_batch!(state::TrainingState, indices)
    data, batch_data = state.full.a[1], state.batch.a[1]
    for (column, sample) in enumerate(indices)
        for row in axes(batch_data, 1)
            batch_data[row, column] = data[row, sample]
            state.batch_target[row, column] = state.target[row, sample]
        end
    end
    return state
end

"""
    epoch!(state, opt, indices, j, losses)

Epoch `j` of a training run: copy the samples `indices` into the state's batch buffers, form the
gradient of the loss over that batch, take one step of `opt`, and record in `losses[j]` the loss
over all the samples of `state`. A loss that is not finite stops the run, naming its epoch.
"""
function epoch!(state::TrainingState, opt, indices, j::Integer, losses)
    _fill_batch!(state, indices)
    optimization_step!(state.ps, opt, gradient!(state))
    loss = loss!(state.full, state, state.target)
    isfinite(loss) || throw(ArgumentError("the loss of epoch $j is $loss"))
    losses[j] = loss
    return losses
end

"""
    _check_batches(batches, batchsize, epochs, samples)

Refuse a set of batches that the run cannot take: one batch per epoch, `batchsize` indices in each,
each of them a sample of the data. The message of a bad batch names its epoch.
"""
function _check_batches(batches, batchsize::Integer, epochs::Integer, samples::Integer)
    length(batches) == epochs || throw(DimensionMismatch(
        "batches hold $(length(batches)) epochs; the run takes $epochs"))
    for (j, indices) in enumerate(batches)
        length(indices) == batchsize || throw(ArgumentError(
            "the batch of epoch $j holds $(length(indices)) indices; a batch holds $batchsize"))
        all(i -> 1 ≤ i ≤ samples, indices) || throw(ArgumentError(
            "the batch of epoch $j holds an index outside 1:$samples"))
    end
    return batches
end

"""
    _epochs!(state, opt, batches, rng, losses)

The epochs of a run. The indices of a drawn batch are filled into one buffer, allocated once for
the run, so that an epoch that draws its batch allocates what an epoch that is given one allocates.
"""
function _epochs!(state::TrainingState, opt, ::Nothing, rng, losses::AbstractVector)
    indices = Vector{Int}(undef, size(state.batch_target, 2))
    for j in eachindex(losses)
        rand!(rng, indices, 1:size(state.target, 2))
        epoch!(state, opt, indices, j, losses)
    end
    return state.ps
end

function _epochs!(state::TrainingState, opt, batches::AbstractVector, rng,
        losses::AbstractVector)
    for j in eachindex(losses)
        epoch!(state, opt, batches[j], j, losses)
    end
    return state.ps
end

function train_hnn(model::Chain, ps₀::NetworkParameters{T}, data::AbstractMatrix,
        target::AbstractMatrix; η, epochs, batchsize, rng, batches = nothing) where {T}
    for (l, layer) in enumerate(model.layers)
        layer isa Dense && layer.σ isa Union{TanhActivation, Identity} ||
            throw(ArgumentError(
                "layer $l of the model is a $(typeof(layer)); the trainer takes `Dense` layers of " *
                "`tanh` or `identity`"))
    end
    inputs, outputs = input_dimension(first(model.layers)),
    output_dimension(last(model.layers))
    samples = size(data, 2)
    size(target) == size(data) || throw(DimensionMismatch(
        "the target is of size $(size(target)); the data are of size $(size(data))"))
    size(data, 1) == inputs || throw(DimensionMismatch(
        "the data hold $(size(data, 1)) reduced positions; the model takes $inputs"))
    outputs == 1 || throw(ArgumentError(
        "the model maps to $outputs outputs; the Hamiltonian is one"))
    samples > 0 || throw(ArgumentError("the data hold no sample"))
    batchsize ≥ 1 ||
        throw(ArgumentError("the batch size is $batchsize; a batch holds a sample"))
    epochs ≥ 1 ||
        throw(ArgumentError("the number of epochs is $epochs; a run has an epoch"))
    isfinite(η) && η > 0 ||
        throw(ArgumentError("the step size is $η; it is finite and positive"))
    eltype(data) === T && eltype(target) === T || throw(ArgumentError(
        "the data are of element type $(eltype(data)), the target $(eltype(target)) and the " *
        "parameters $T; one run is of one element type"))
    batches === nothing || _check_batches(batches, batchsize, epochs, samples)

    # the run trains its own copy, which the state and the optimizer both hold
    state = TrainingState(model, deepcopy(ps₀), data, target, batchsize)
    opt = TrainingOptimizer(state.ps; algorithm = GradientMethod(), linesearch = η)
    losses = zeros(T, epochs)
    return _epochs!(state, opt, batches, rng, losses), losses
end

end
