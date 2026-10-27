# The reference run of the HNN record: the bump-on-tail trainer's algorithm in plain Julia, with
# nested ForwardDiff for both of its derivatives, and no machine-learning framework. It shares no
# code with the trainer it checks: the network is three matrix products over weights laid out
# `(out, in)`, its input gradient and its parameter gradient are both forward mode, and the step is
# `p .-= η .* dp`.
#
# It reads the two files that `scripts/bump_on_tail_projections.jl` writes at its defaults: the
# particle spectrum `Λ` of `<run id>_projections.h5`, and the reduced positions `X` and
# accelerations `A` of `<run id>_vectorfield.h5`. It trains on the leading `m` rows of the latter —
# one coefficient per row, in the column order of `reshape(X, k, :)` — and writes `record.h5`
# beside this file, holding the initial weights, the batch of every epoch and the loss of the whole
# data set after each step, so that the package's trainer can be checked against a run that reaches
# its numbers another way.
#
#     julia --project=scripts/hnn_reference scripts/hnn_reference/train_reference.jl --runs <dir>
#
# The defaults are those of the record: `m = 10`, 1000 epochs of 100 samples, `η = 1e-4` and the
# seed 42. The environment is the one beside this file, and it names ForwardDiff, HDF5 and the
# standard library only.

using ForwardDiff
using HDF5
using Random

"The defaults of the record, and the meaning of every option."
const DEFAULTS = (m = 10, epochs = 1000, batchsize = 100, eta = 1e-4, seed = 42)

const USAGE = "train_reference.jl --runs <dir> [--m 10] [--epochs 1000] [--batchsize 100] " *
              "[--eta 1.0e-4] [--seed 42] [--out <file>]"

"""
    Network

The network of the record as the five arrays of its parameter set: `W1` and `W2` of size `(m, m)`
with their biases `b1` and `b2`, and `W3` of size `(1, m)` with no bias.
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

The network of the flat parameter vector `p`, `[vec(W1); b1; vec(W2); b2; vec(W3)]`: the order of
`NeuralNetworkParameters.flatten` on the `Chain` of the three `Dense` layers, so that the record's
weights and the package's parameter set are one layout.
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
    glorot_uniform(rng, out, in)

The weights of one layer, drawn uniform in `±sqrt(6 / (out + in))`: the Glorot uniform
initialisation of a dense layer.
"""
function glorot_uniform(rng, out::Integer, in::Integer)
    limit = sqrt(6 / (out + in))
    return limit .* (2 .* rand(rng, out, in) .- 1)
end

"""
    initial_parameters(rng, m)

The flat parameter vector of the record's initial weights: the two weight matrices and the output
weights drawn Glorot-uniform from `rng`, and the two biases zero.
"""
function initial_parameters(rng, m::Integer)
    return [vec(glorot_uniform(rng, m, m)); zeros(m); vec(glorot_uniform(rng, m, m));
            zeros(m);
            vec(glorot_uniform(rng, 1, m))]
end

"""
    train(data, target, p, η, batches) -> (p, losses)

The training run: for each epoch `j`, one step `p .-= η .* dp` of the parameter gradient of the loss
over the samples of `batches[j]`, and then `losses[j]`, the loss over all the samples of `data` and
`target` after the step.
"""
function train(data::AbstractMatrix, target::AbstractMatrix, p, η, batches)
    m = size(data, 1)
    losses = zeros(length(batches))
    for j in eachindex(batches)
        idx = batches[j]
        p = p .- η .* parameter_gradient(p, data[:, idx], target[:, idx])
        losses[j] = data_loss(network(p, m), data, target)
    end
    return p, losses
end

"""
    read_data(runs, m)

The leading `m` rows of the reduced positions and accelerations of the one run under `runs`, the
leading `m` values of its spectrum, and the commit its files record.
"""
function read_data(runs::AbstractString, m::Integer)
    names = readdir(runs)
    vectorfields = filter(endswith("_vectorfield.h5"), names)
    length(vectorfields) == 1 || throw(ArgumentError(
        "the directory $runs holds $(length(vectorfields)) files named `_vectorfield.h5`; the " *
        "run of the record is one"))
    vectorfield = only(vectorfields)
    projections = replace(vectorfield, "_vectorfield.h5" => "_projections.h5")
    projections in names || throw(ArgumentError(
        "$(joinpath(runs, projections)) is not there; the spectrum of the run is in it"))

    X = h5open(io -> read(io["X"]), joinpath(runs, vectorfield), "r")
    A = h5open(io -> read(io["A"]), joinpath(runs, vectorfield), "r")
    λ = h5open(io -> read(io["Λ"]), joinpath(runs, projections), "r")
    commit = h5open(io -> read(attributes(io)["commit"]), joinpath(runs, vectorfield), "r")

    size(X, 1) ≥ m || throw(ArgumentError(
        "the run holds $(size(X, 1)) coefficients; the record trains on $m"))
    inputs = reshape(X[1:m, :, :], m, :)
    targets = reshape(A[1:m, :, :], m, :)
    return inputs, targets, λ[1:m], commit
end

"""
    write_record(file, o, inputs, targets, λ, initial, losses, batches, command)

Write the record: the raw inputs and targets of the training run, the leading spectrum, the initial
weights, the losses and the batches, and the attributes that name the run.
"""
function write_record(file, o, inputs, targets, λ, initial, losses, batches, command)
    m = o.m
    net = network(initial, m)
    h5open(file, "w") do io
        a = attributes(io)
        a["dimension"] = m
        a["epochs"] = o.epochs
        a["batchsize"] = o.batchsize
        a["eta"] = o.eta
        a["seed"] = o.seed
        a["command"] = command
        a["input_commit"] = o.commit
        a["forwarddiff_version"] = string(pkgversion(ForwardDiff))

        io["inputs"] = inputs
        io["targets"] = targets
        io["lambda"] = λ
        io["losses"] = losses
        io["batches"] = batches

        g = create_group(io, "initial")
        g["W1"] = net.W1
        g["b1"] = net.b1
        g["W2"] = net.W2
        g["b2"] = net.b2
        g["W3"] = net.W3
    end
    return file
end

"""
    parse_arguments(args)

The options of the run, from the `--name value` pairs of `args`, at the defaults of `DEFAULTS`.
"""
function parse_arguments(args)
    o = DEFAULTS
    runs = ""
    out = joinpath(@__DIR__, "record.h5")
    i = firstindex(args)
    while i ≤ lastindex(args)
        name = args[i]
        i < lastindex(args) || throw(ArgumentError("$name has no value; $USAGE"))
        value = args[i + 1]
        if name == "--runs"
            runs = value
        elseif name == "--out"
            out = value
        elseif name == "--m"
            o = merge(o, (m = parse(Int, value),))
        elseif name == "--epochs"
            o = merge(o, (epochs = parse(Int, value),))
        elseif name == "--batchsize"
            o = merge(o, (batchsize = parse(Int, value),))
        elseif name == "--eta"
            o = merge(o, (eta = parse(Float64, value),))
        elseif name == "--seed"
            o = merge(o, (seed = parse(Int, value),))
        else
            throw(ArgumentError("unknown option $name; $USAGE"))
        end
        i += 2
    end
    isempty(runs) && throw(ArgumentError("the runs directory is required; $USAGE"))
    return merge(o, (runs = runs, out = out))
end

function main(args)
    o = parse_arguments(args)
    inputs, targets, λ, commit = read_data(o.runs, o.m)
    o = merge(o, (commit = commit,))

    rng = Xoshiro(o.seed)
    initial = initial_parameters(rng, o.m)
    batches = [rand(rng, 1:size(inputs, 2), o.batchsize) for _ in 1:o.epochs]

    # the run trains on the inputs scaled by the square root of the spectrum, and on the raw
    # accelerations
    data = inputs ./ sqrt.(λ)
    _, losses = train(data, targets, initial, o.eta, batches)

    write_record(o.out, o, inputs, targets, λ, initial, losses, reduce(hcat, batches),
        join([basename(PROGRAM_FILE); args], " "))
    return println("losses: ", losses[1], " … ", losses[end], ", record: ", o.out)
end

# The run happens when this file is the program: `include`ing it defines its functions and runs
# nothing, so that a test or a probe can use the reference against the trainer it checks.
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
