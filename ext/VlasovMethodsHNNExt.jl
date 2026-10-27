module VlasovMethodsHNNExt

using AbstractNeuralNetworks: AbstractNeuralNetworks, Chain, Dense, NeuralNetwork,
                              input_dimension, output_dimension, params
using GeometricOptimizers
using SymbolicNeuralNetworks: SymbolicNeuralNetwork, build_nn_function, derivative,
                              Jacobian,
                              symbolic_parameter_gradient, symbolic_variables

import VlasovMethods: train_hnn

# The trainer of `src/hnn.jl`: the gradient of the loss comes from `SymbolicNeuralNetworks`, which
# builds it symbolically once and compiles it, and the step from `GeometricOptimizers`. The three
# packages are weak dependencies, so this module loads only when all of them are loaded, and
# nothing of the symbolic engine reaches a caller who does not load them.

"""
    TrainingState

Everything an epoch reads and writes, built once for the run: the two compiled functions of the
network, the training data and its target, and the buffers of one batch. An epoch writes the
samples of its batch into the buffers and reads the whole data set for its loss, so it allocates
nothing but what those calls allocate themselves.
"""
struct TrainingState{T, G, L}
    gradient::G
    loss::L
    data::Matrix{T}
    target::Matrix{T}
    batch_data::Matrix{T}
    batch_target::Matrix{T}
end

"""
    TrainingState(model, data, target, batchsize)

The state of a training run of the network `model` over the columns of `data` and `target`: both
are copied, so that the state holds the run's own data for the whole of it, and the batch buffers
have one column per sample of a batch.
"""
function TrainingState(model::Chain, data::AbstractMatrix{T}, target::AbstractMatrix{T},
        batchsize::Integer) where {T}
    gradient, loss = _functions(model, T)
    batch_buffer() = Matrix{T}(undef, size(data, 1), batchsize)
    return TrainingState{T, typeof(gradient), typeof(loss)}(gradient, loss, copy(data),
        copy(target), batch_buffer(), batch_buffer())
end

"""
    _functions(model, T)

The two compiled functions of `model`: the batch gradient of the loss with respect to the
parameters, which sums over the samples of its batch, and the loss itself, a scalar. Both take the
batch as two matrices of samples in their columns — the inputs and their targets — and the
parameters.

The input gradient of the Hamiltonian is the network's Jacobian with respect to its input, whose
single row is the gradient of the single output; the target enters as the symbolic variables `y`,
so that the loss and its parameter gradient are one expression each.
"""
function _functions(model::Chain, ::Type{T}) where {T}
    snn = SymbolicNeuralNetwork(model)
    y = symbolic_variables(zeros(T, input_dimension(model)), :y)
    ∇H = vec(derivative(Jacobian(snn)))
    loss = sum((∇H .- y) .^ 2)
    batch_gradient = symbolic_parameter_gradient(loss, snn)
    gradient = build_nn_function(batch_gradient, params(snn), snn.input, y; reduce = +)
    value = build_nn_function(loss, params(snn), snn.input, y; reduce = +)
    return gradient, value
end

"""
    _fill_batch!(batch_data, batch_target, data, target, indices)

Copy the samples `indices` of `data` and `target` into the batch buffers, one sample per column.
"""
function _fill_batch!(batch_data::AbstractMatrix, batch_target::AbstractMatrix,
        data::AbstractMatrix, target::AbstractMatrix, indices)
    for (column, sample) in enumerate(indices)
        for row in axes(batch_data, 1)
            batch_data[row, column] = data[row, sample]
            batch_target[row, column] = target[row, sample]
        end
    end
    return batch_data, batch_target
end

"""
    epoch!(ps, opt, state, indices, j, losses)

Epoch `j` of a training run: copy the samples `indices` into the state's batch buffers, form the
gradient of the loss over that batch, take one step of `opt`, and record in `losses[j]` the loss
over all the samples of `state`. A loss that is not finite stops the run, naming its epoch.
"""
function epoch!(ps, opt, state::TrainingState, indices, j::Integer, losses)
    _fill_batch!(state.batch_data, state.batch_target, state.data, state.target, indices)
    dp = state.gradient(state.batch_data, state.batch_target, ps)
    optimization_step!(ps, opt, dp)
    loss = state.loss(state.data, state.target, ps)
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

function train_hnn(model::Chain, ps₀::NetworkParameters{T}, data::AbstractMatrix,
        target::AbstractMatrix; η, epochs, batchsize, rng, batches = nothing) where {T}
    samples = size(data, 2)
    size(target) == size(data) || throw(DimensionMismatch(
        "the target is of size $(size(target)); the data are of size $(size(data))"))
    size(data, 1) == input_dimension(model) || throw(DimensionMismatch(
        "the data hold $(size(data, 1)) reduced positions; the model takes " *
        "$(input_dimension(model))"))
    output_dimension(model) == 1 || throw(ArgumentError(
        "the model maps to $(output_dimension(model)) outputs; the Hamiltonian is one"))
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

    state = TrainingState(model, data, target, batchsize)
    # the optimizer keeps its own state of the parameters it is built on, so it is built on a copy
    ps = deepcopy(ps₀)
    opt = TrainingOptimizer(ps; algorithm = GradientMethod(), linesearch = η)
    losses = zeros(T, epochs)
    for j in 1:epochs
        indices = batches === nothing ? rand(rng, 1:samples, batchsize) : batches[j]
        epoch!(ps, opt, state, indices, j, losses)
    end
    return ps, losses
end

end
