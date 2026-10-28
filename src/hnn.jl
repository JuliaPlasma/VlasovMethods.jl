# The trainer of the reduced field is a package extension. This file declares the function and
# documents it; `ext/VlasovMethodsHNNExt.jl` holds its only method, so that loading the package
# loads none of the machine-learning packages, and `train_hnn` is a function with no methods until
# the extension loads.

@doc raw"""
    train_hnn(model, ps₀, data, target; η, epochs, batchsize, rng, batches = nothing)
        -> (ps, losses)

Train the Hamiltonian neural network `model` on the reduced field of a bump-on-tail sweep, and
return its parameters and the loss of every epoch.

Defined in a package extension: load `GeometricOptimizers` and `AbstractNeuralNetworks` to make it
available.

```julia
using VlasovMethods, GeometricOptimizers, AbstractNeuralNetworks
m = 10
model = Chain(Dense(m, m, tanh), Dense(m, m, tanh), Dense(m, 1, identity; use_bias = false))
ps, losses = train_hnn(model, params(NeuralNetwork(model)), data, target;
    η = 1e-4, epochs = 1000, batchsize = 100, rng = Xoshiro(42))
```

`model` is the Hamiltonian `H` of the reduced field, a function of the reduced positions, and
`model(x, ps)` its value at one sample. It is a `Chain` of `Dense` layers whose activations are
`tanh` or `identity`, with one output. `data` holds one sample per column, the leading
coefficients of the particle basis divided by the square root of their values of its spectrum, and
`target` the matching accelerations, of the same size. The loss of a sample is
``||\nabla_x H(x) - y||^2``: the gradient of the Hamiltonian with respect to its input against the
acceleration. The loss of an epoch is the sum of that over the samples of the batch it steps on,
and `losses[j]` is the loss over all `n` samples after the step of epoch `j`.

Epoch `j` forms the gradient of the loss over its batch with respect to the parameters, takes one
`GeometricOptimizers.optimization_step!` of a `TrainingOptimizer` built with
`GradientMethod()` and `linesearch = η`, and records the loss over all samples. `batches` holds
one batch per epoch as the indices of the samples it steps on, its batches drawing with
replacement; without it, epoch `j` draws its batch as `rand(rng, 1:n, batchsize)`.

`ps₀` is a `NetworkParameters` of `model`, and the trainer copies it: the returned parameters are
a new set, and `ps₀` holds the values it held at the call. `data`, `target` and `ps₀` are of one
element type, which the parameters and losses are returned in.

# Throws

- `DimensionMismatch` when `data` and `target` do not have the same size, when the first axis of
  `data` is not the input dimension of `model`, and when `batches` does not hold one batch per
  epoch.
- `ArgumentError` when a layer of `model` is not a `Dense` layer of `tanh` or `identity`, when
  `model` maps to more than one output, when `data` holds no sample, when
  `batchsize` or `epochs` is below one, when `η` is not finite and positive, when a batch of
  `batches` is not of `batchsize` indices in `1:n`, when the element types differ, and when a loss
  is not finite, which names the epoch.
"""
function train_hnn end
