# Bump-on-tail, the HNN of the reduced field: the loss curve of the Hamiltonian network trained on
# the projections.
#
# Reads the two files that `bump_on_tail_projections.jl` writes with the same options, takes the
# leading `--m` coefficients of every saved state as the training data, and trains the three-layer
# network on them with `train_hnn`, whose step comes from GeometricOptimizers. The default command:
#
#     julia --project=scripts scripts/bump_on_tail_train_hnn.jl
#
# Writes `results_path("hnn_loss.png")`: the loss over the whole training set after each epoch. The
# run draws its initial weights and its batches from `--seed`, so two runs at the same options draw
# one curve.

include("options.jl")
include("bump_on_tail_common.jl")

using CairoMakie
using AbstractNeuralNetworks: Chain, Dense, NeuralNetwork, params
using GeometricOptimizers

const OPT = parse_arguments(ARGS; BUMP_ON_TAIL..., m = 10, epochs = 1000, batchsize = 100,
    eta = 1e-4)

function main(o)
    id = run_id(o)
    Λ = h5open(io -> ReducedBasis(io).Λ, runs_path(id * "_projections.h5"), "r")
    X = h5open(io -> read(io["X"]), runs_path(id * "_vectorfield.h5"), "r")
    A = h5open(io -> read(io["A"]), runs_path(id * "_vectorfield.h5"), "r")
    size(X, 1) ≥ o.m || throw(ArgumentError(
        "the run holds $(size(X, 1)) coefficients; the network trains on $(o.m)"))

    # the columns of `reshape(X, k, :)`: one sample per saved state, in the order of the parameters
    # and of the times
    data, target = hnn_data(
        reshape(X[1:o.m, :, :], o.m, :), reshape(A[1:o.m, :, :], o.m, :),
        Λ[1:o.m])
    model = Chain(Dense(o.m, o.m, tanh), Dense(o.m, o.m, tanh),
        Dense(o.m, 1, identity; use_bias = false))

    _, losses = train_hnn(model, params(NeuralNetwork(model; rng = Xoshiro(o.seed))), data,
        target; η = o.eta, epochs = o.epochs, batchsize = o.batchsize, rng = Xoshiro(o.seed))
    println("loss: $(losses[1]) … $(losses[end])")

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "epoch", ylabel = "loss", yscale = log10)
    lines!(ax, 1:o.epochs, losses)
    save(results_path("hnn_loss.png"), fig)
    return nothing
end

main(OPT)
