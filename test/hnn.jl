using Test
using TOML
using ForwardDiff
using HDF5
using Random
using AbstractNeuralNetworks: Chain, Dense, NeuralNetwork, params
using GeometricOptimizers: GradientMethod, NetworkParameters, TrainingOptimizer,
                           optimization_step!
using SymbolicNeuralNetworks
using VlasovMethods

# `train_hnn` is declared in `src/hnn.jl` and its one method is in `ext/VlasovMethodsHNNExt.jl`,
# which loads with the three packages above. The method is checked against
# `scripts/hnn_reference/record.h5`, the run of a trainer in `scripts/hnn_reference/` that shares no
# code with it: at the record's initial weights, its data and its batches, the two must produce one
# loss curve.
#
# The other half of the split is the fresh-environment testset below, which runs in a subprocess
# whose own manifest names none of the three packages, so no check here can see it.

include(joinpath(@__DIR__, "..", "scripts", "bump_on_tail_common.jl"))

# the reference trainer, in a module of its own: it shares no name and no code with the package, and
# including it here runs nothing
module Reference
include(joinpath(@__DIR__, "..", "scripts", "hnn_reference", "train_reference.jl"))
end

const RECORD = joinpath(@__DIR__, "..", "scripts", "hnn_reference", "record.h5")

"The network of the trainer: two hidden layers of `tanh`, and a linear output with no bias."
function hnn_model(m)
    Chain(Dense(m, m, tanh), Dense(m, m, tanh), Dense(m, 1, identity; use_bias = false))
end

"The parameter set of `hnn_model(m)` holding the weights `W1`, `b1`, `W2`, `b2` and `W3`."
function network_parameters(W1, b1, W2, b2, W3)
    return NetworkParameters(
        (L1 = (W = W1, b = b1), L2 = (W = W2, b = b2), L3 = (W = W3,)))
end

"The refusal of a call, as its message, so that a check can name the guard that refused."
function refusal(f)
    try
        f()
    catch exception
        return exception
    end
    return nothing
end

@testset "the reference trainer's network" begin
    m = 3
    model = hnn_model(m)
    ps = params(NeuralNetwork(model))
    net = Reference.Network(ps.L1.W, ps.L1.b, ps.L2.W, ps.L2.b, ps.L3.W)
    x = randn(Xoshiro(2), m)

    # the two forward passes and the two input gradients are the same function, written apart
    @test Reference.potential(net, x) ≈ sum(model(x, ps))
    @test Reference.input_gradient(net, x) ≈
          ForwardDiff.gradient(xi -> sum(model(xi, ps)), x)

    # the options of the run, which the record's command line went through
    @test Reference.parse_arguments(["--runs", "here", "--m", "4"]).m == 4
    @test_throws ArgumentError Reference.parse_arguments(String[])
    @test_throws ArgumentError Reference.parse_arguments(["--runs"])
    @test_throws ArgumentError Reference.parse_arguments(["--unknown", "1"])
end

@testset "the reference trainer's data" begin
    mktempdir() do dir
        X = reshape(collect(1.0:24.0), 2, 3, 4)
        A = 2 .* X
        λ = [1.0, 2.0]
        h5open(joinpath(dir, "run_vectorfield.h5"), "w") do io
            io["X"] = X
            io["A"] = A
            attributes(io)["commit"] = "c"^40
        end
        h5open(joinpath(dir, "run_projections.h5"), "w") do io
            io["Λ"] = λ
        end

        inputs, targets, lambda, commit = Reference.read_data(dir, 2)
        # the columns of `reshape(X, k, :)`: one sample per saved state, in the order of the
        # parameters and of the times
        @test inputs == reshape(X, 2, :)
        @test targets == reshape(A, 2, :)
        @test lambda == λ && commit == "c"^40
        @test_throws ArgumentError Reference.read_data(dir, 3)
        @test_throws ArgumentError Reference.read_data(mktempdir(), 1)
    end
end

@testset "the curve of the record" begin
    h5open(RECORD, "r") do io
        attributes_ = attributes(io)
        m = read(attributes_["dimension"])
        # the width the clause and this test are written for: at a larger width the run's loss
        # leaves the `1e-8` bound of the record's curve in the last few hundred epochs
        @test m == 6
        epochs = read(attributes_["epochs"])
        batchsize = read(attributes_["batchsize"])
        η = read(attributes_["eta"])
        inputs, targets = read(io["inputs"]), read(io["targets"])
        λ, recorded = read(io["lambda"]), read(io["losses"])
        batches = read(io["batches"])
        weights = ("W1", "b1", "W2", "b2", "W3")
        W1, b1, W2, b2, W3 = map(k -> read(io["initial/$k"]), weights)

        model = hnn_model(m)
        ps₀ = network_parameters(W1, b1, W2, b2, W3)
        reference = params(NeuralNetwork(model))
        @test keys(ps₀) == keys(reference) == (:L1, :L2, :L3)
        @test keys(ps₀.L3) == (:W,)
        @test size(ps₀.L1.W) == size(reference.L1.W) == (m, m)
        @test size(ps₀.L3.W) == size(reference.L3.W) == (1, m)
        data, target = hnn_data(inputs, targets, λ)
        batches = [batches[:, j] for j in axes(batches, 2)]

        _, losses = train_hnn(model, ps₀, data, target;
            η = η, epochs = epochs, batchsize = batchsize, rng = Xoshiro(0), batches = batches)

        # one step of the trainer is one step of the reference: the loss of the first epoch agrees to
        # fourteen digits, which is the round-off of two summations of the same 2510 terms
        @test isapprox(losses[1], recorded[1]; rtol = 1e-12, atol = 0)

        # every loss agrees with the record to `1e-8`, the bound of the clause, which is stated over
        # the first 400. The margin is wide: the largest relative difference over the 1000 epochs is
        # about `2e-14`, which leaves room for another machine's summation order and for another BLAS
        @test all(isapprox.(losses, recorded; rtol = 1e-8, atol = 0))

        # the control: the same run at a step size that differs in the fourth digit, which the
        # bounds above see at the first epoch
        _, other = train_hnn(model, ps₀, data, target;
            η = 1.01e-4, epochs = epochs, batchsize = batchsize, rng = Xoshiro(0),
            batches = batches)
        @test !isapprox(other[1], recorded[1]; rtol = 1e-8, atol = 0)
    end
end

@testset "the edges of train_hnn" begin
    m, n, batchsize = 3, 20, 4
    model = hnn_model(m)
    rng = Xoshiro(3)
    data = randn(rng, m, n)
    target = randn(rng, m, n)
    ps₀ = params(NeuralNetwork(model))
    batches = [rand(Xoshiro(11), 1:n, batchsize) for _ in 1:3]
    keywords = (;
        η = 0.1, epochs = 3, batchsize = batchsize, rng = Xoshiro(5), batches = batches)

    @testset "the shapes" begin
        err = refusal() do
            train_hnn(model, ps₀, data, target[:, 1:(end - 1)]; keywords...)
        end
        @test err isa DimensionMismatch
        @test occursin("$(size(data))", err.msg) &&
              occursin("$(size(target[:, 1:(end - 1)]))", err.msg)

        err = refusal() do
            train_hnn(model, ps₀, randn(rng, m + 1, n), randn(rng, m + 1, n); keywords...)
        end
        @test err isa DimensionMismatch
        @test occursin("$(m + 1) reduced positions", err.msg) &&
              occursin("takes $m", err.msg)

        err = refusal() do
            train_hnn(model, ps₀, randn(rng, m - 1, n), randn(rng, m - 1, n); keywords...)
        end
        @test err isa DimensionMismatch
        @test occursin("$(m - 1) reduced positions", err.msg) &&
              occursin("takes $m", err.msg)

        wide = Chain(Dense(m, m, tanh), Dense(m, 2))
        err = refusal() do
            train_hnn(wide, ps₀, data, target; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("maps to 2 outputs", err.msg)
    end

    @testset "the sizes of a run" begin
        empty_data = zeros(m, 0)
        err = refusal() do
            train_hnn(model, ps₀, empty_data, empty_data; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("hold no sample", err.msg)

        for (keywords_, phrase) in (
            ((; keywords..., batchsize = 0), "the batch size is 0"),
            ((; keywords..., epochs = 0), "the number of epochs is 0"),
            ((; keywords..., η = Inf), "the step size is Inf"),
            ((; keywords..., η = 0.0), "the step size is 0.0"),
            ((; keywords..., η = -1e-4), "the step size is -0.0001"))
            err = refusal() do
                train_hnn(model, ps₀, data, target; keywords_...)
            end
            @test err isa ArgumentError
            @test occursin(phrase, err.msg)
        end
    end

    @testset "the batches" begin
        err = refusal() do
            train_hnn(model, ps₀, data, target; keywords..., batches = batches[1:2])
        end
        @test err isa DimensionMismatch
        @test occursin("hold 2 epochs", err.msg) && occursin("takes 3", err.msg)

        err = refusal() do
            train_hnn(model, ps₀, data, target;
                keywords..., batches = [
                    rand(Xoshiro(2), 1:n, batchsize + 1), batches[2:3]...])
        end
        @test err isa ArgumentError
        @test occursin("the batch of epoch 1 holds $(batchsize + 1) indices", err.msg)
        @test occursin("a batch holds $batchsize", err.msg)

        err = refusal() do
            train_hnn(model, ps₀, data, target;
                keywords..., batches = [batches[1], [1, 2, 3, n + 1], batches[3]])
        end
        @test err isa ArgumentError
        @test occursin("the batch of epoch 2", err.msg) && occursin("outside 1:$n", err.msg)
    end

    @testset "the element types" begin
        err = refusal() do
            train_hnn(model, ps₀, Float32.(data), target; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("Float32", err.msg) && occursin("Float64", err.msg)

        err = refusal() do
            train_hnn(model, ps₀, data, Float32.(target); keywords...)
        end
        @test err isa ArgumentError
        @test occursin("the target Float32", err.msg)
    end

    @testset "a loss that is not finite" begin
        target_nan = copy(target)
        target_nan[1, 1] = NaN
        err = refusal() do
            train_hnn(model, ps₀, data, target_nan; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("the loss of epoch 1 is", err.msg) && occursin("NaN", err.msg)
    end

    @testset "the parameters of the call are not the run's" begin
        before = deepcopy(ps₀)
        train_hnn(model, ps₀, data, target; keywords...)
        @test ps₀ == before
    end
end

@testset "the drawn batches and the other precision" begin
    m, n, epochs, batchsize = 3, 20, 3, 5
    model = hnn_model(m)
    rng = Xoshiro(9)
    data = randn(rng, m, n)
    target = randn(rng, m, n)
    ps₀ = params(NeuralNetwork(model))
    keywords = (; η = 1e-2, epochs = epochs, batchsize = batchsize)

    # each run draws from a generator of its own, so the three are the same run: the drawn batches
    # are the batches it is given
    draw = Xoshiro(7)
    drawn = [rand(draw, 1:n, batchsize) for _ in 1:epochs]
    first_run = train_hnn(model, ps₀, data, target; keywords..., rng = Xoshiro(7))
    second_run = train_hnn(model, ps₀, data, target; keywords..., rng = Xoshiro(7))
    given_run = train_hnn(model, ps₀, data, target; keywords..., rng = Xoshiro(7),
        batches = drawn)
    @test first_run == second_run == given_run
    @test length(first_run[2]) == epochs && all(>(0), first_run[2])

    data32, target32 = Float32.(data), Float32.(target)
    ps₀32 = network_parameters(Float32.(ps₀.L1.W), Float32.(ps₀.L1.b), Float32.(ps₀.L2.W),
        Float32.(ps₀.L2.b), Float32.(ps₀.L3.W))
    ps₃₂, losses32 = train_hnn(
        model, ps₀32, data32, target32; keywords..., rng = Xoshiro(7),
        batches = drawn)
    @test ps₃₂.L1.W isa Matrix{Float32} && losses32 isa Vector{Float32}
end

@testset "the allocations of a training step" begin
    # `--check-bounds=yes` frames a copy of the arrays, so the reading is taken elsewhere
    if Base.JLOptions().code_coverage == 0
        extension = Base.get_extension(VlasovMethods, :VlasovMethodsHNNExt)
        m, n, batchsize, epochs = 3, 20, 4, 3
        model = hnn_model(m)
        rng = Xoshiro(4)
        data = randn(rng, m, n)
        target = randn(rng, m, n)
        state = extension.TrainingState(model, data, target, batchsize)
        ps = deepcopy(params(NeuralNetwork(model)))
        opt = TrainingOptimizer(ps; algorithm = GradientMethod(), linesearch = 1e-2)
        indices = [1, 4, 7, 10]
        losses = zeros(epochs)

        # a fixed-arity barrier of concrete arguments, which calls once before it measures
        function epoch_allocations(f, ps, opt, state, indices, j, losses)
            f(ps, opt, state, indices, j, losses)
            return @allocated f(ps, opt, state, indices, j, losses)
        end
        function call_allocations(f, a, b, c)
            f(a, b, c)
            return @allocated f(a, b, c)
        end
        function step_allocations(f, ps, opt, dp)
            f(ps, opt, dp)
            return @allocated f(ps, opt, dp)
        end

        dp = state.gradient(state.batch_data, state.batch_target, ps)
        @test step_allocations(optimization_step!, ps, opt, dp) == 0

        epoch = epoch_allocations(extension.epoch!, ps, opt, state, indices, 1, losses)
        gradient = call_allocations(state.gradient, state.batch_data, state.batch_target, ps)
        loss = call_allocations(state.loss, state.data, state.target, ps)

        # the SNN calls allocate their results, and the epoch adds nothing to them
        @test gradient > 0 && loss > 0                       # the controls
        @test epoch > 0
        @test epoch ≤ gradient + loss
    end
end

@testset "the docstring names the extension" begin
    @test occursin("package extension", string(@doc VlasovMethods.train_hnn))
    @test occursin("load `SymbolicNeuralNetworks`", string(@doc VlasovMethods.train_hnn))
end

@testset "the three packages are weak dependencies" begin
    project = TOML.parsefile(joinpath(pkgdir(VlasovMethods), "Project.toml"))

    @test get(project["extensions"], "VlasovMethodsHNNExt", nothing) ==
          ["AbstractNeuralNetworks", "GeometricOptimizers", "SymbolicNeuralNetworks"]
    for name in ("AbstractNeuralNetworks", "GeometricOptimizers", "SymbolicNeuralNetworks")
        @test haskey(project["weakdeps"], name)
        @test haskey(project["compat"], name)
        @test !haskey(project["deps"], name)
    end

    # neither the symbolic engine nor the machine-learning stack is a dependency, so a caller who
    # loads none of the three loads none of them
    for name in ("Optimisers", "Zygote", "Lux", "Symbolics", "NeuralNetworkParameters")
        @test !haskey(project["deps"], name)
    end
end

@testset "a fresh environment without the three packages" begin
    program = """
        import Pkg
        Pkg.activate(; temp = true, io = devnull)
        Pkg.develop(path = $(repr(pkgdir(VlasovMethods))); io = devnull)

        # a weak dependency of the package is not in the manifest of an environment that names
        # neither it nor a package that depends on it
        dependencies = Pkg.dependencies()
        for (uuid, name) in (
            ("aed23131-dcd0-47ca-8090-d21e605652e3", "SymbolicNeuralNetworks"),
            ("fc236c15-5557-4942-aa65-b650f329279e", "GeometricOptimizers"),
            ("60874f82-5ada-4c70-bd1c-fa6be7711c8a", "AbstractNeuralNetworks"))
            @assert !haskey(dependencies, Base.UUID(uuid)) "\$(name) is in the environment"
            @assert !any(info.name == name for info in values(dependencies)) "\$(name) is in the environment"
        end

        using VlasovMethods

        @assert Base.get_extension(VlasovMethods, :VlasovMethodsHNNExt) === nothing
        @assert isempty(methods(VlasovMethods.train_hnn))

        exception = try
            VlasovMethods.train_hnn(nothing, nothing, zeros(2, 2), zeros(2, 2))
            nothing
        catch thrown
            thrown
        end
        @assert exception isa MethodError "an unloaded train_hnn threw \$(typeof(exception))"

        print("no network package in the environment: the extension is not loaded and train_hnn has no method")
        """

    # `JULIA_LOAD_PATH` without `@v1.x`, so that a package of the shared environment cannot satisfy
    # a stray `using` in the program
    load_path = join(("@", "@stdlib"), Sys.iswindows() ? ';' : ':')
    command = addenv(`$(Base.julia_cmd()) --startup-file=no -e $program`,
        "JULIA_LOAD_PATH" => load_path)

    output = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = output, stderr = output))
    success(process) || print(String(take!(output)))
    @test success(process)
end
