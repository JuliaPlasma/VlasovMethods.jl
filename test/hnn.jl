using Test
using TOML
using ForwardDiff
using Random
using AbstractNeuralNetworks: Chain, Dense, NeuralNetwork, params
using GeometricOptimizers: GradientMethod, NetworkParameters, TrainingOptimizer,
                           optimization_step!
using VlasovMethods

# `train_hnn` is declared in `src/hnn.jl` and its one method is in `ext/VlasovMethodsHNNExt.jl`,
# which loads with the two packages above. The method is checked against the reference of
# `test/helpers/hnn_reference.jl`, a trainer that shares no code with it: from the same initial
# weights, on the same data and the same batches, the two must produce one loss curve.
#
# The other half of the split is the fresh-environment testset below, which runs in a subprocess
# whose own manifest names neither of the two packages, so no check here can see it.

# the reference trainer, in a module of its own: it shares no name and no code with the package
module Reference
include(joinpath(@__DIR__, "helpers", "hnn_reference.jl"))
end

"The network of the trainer: two hidden layers of `tanh`, and a linear output with no bias."
function hnn_model(m)
    Chain(Dense(m, m, tanh), Dense(m, m, tanh), Dense(m, 1, identity; use_bias = false))
end

"The parameters of `hnn_model` as the reference's flat vector `[vec(W1); b1; vec(W2); b2; vec(W3)]`."
flat(ps) = [vec(ps.L1.W); ps.L1.b; vec(ps.L2.W); ps.L2.b; vec(ps.L3.W)]

"""
The refusal of a call, as its exception, so that a check can name the guard that refused. The call
is held to leaving `ps` as it held it: a refused call trains nothing, and the caller's parameters
are the ones it passed.
"""
function refusal(f, ps)
    before = deepcopy(ps)
    exception = try
        f()
    catch caught
        caught
    end
    @test ps == before
    return exception
end

@testset "the reference trainer's network" begin
    m = 3
    model = hnn_model(m)
    ps = params(NeuralNetwork(model; rng = Xoshiro(2)))
    net = Reference.network(flat(ps), m)
    x = randn(Xoshiro(2), m)

    # the flat layout is the parameter set's, and the two forward passes and the two input
    # gradients are the same function, written apart
    @test net.W1 == ps.L1.W && net.b2 == ps.L2.b && net.W3 == ps.L3.W
    @test Reference.potential(net, x) ≈ sum(model(x, ps))
    @test Reference.input_gradient(net, x) ≈
          ForwardDiff.gradient(xi -> sum(model(xi, ps)), x)
end

@testset "the curve of the reference" begin
    # the width of the driver's default, on data whose target is near the gradient of `|x|² / 2`,
    # at the driver's step size
    m, n, batchsize, epochs, η = 10, 200, 20, 1000, 1e-4
    model = hnn_model(m)
    rng = Xoshiro(1)
    data = randn(rng, m, n)
    target = data .+ 0.1 .* randn(rng, m, n)
    ps₀ = params(NeuralNetwork(model; rng = Xoshiro(2)))
    batches = [rand(Xoshiro(3), 1:n, batchsize) for _ in 1:epochs]

    _, losses = train_hnn(model, ps₀, data, target;
        η = η, epochs = epochs, batchsize = batchsize, rng = Xoshiro(0), batches = batches)
    _, reference = Reference.train(data, target, flat(ps₀), η, batches)

    # the run trains: the loss falls by a third over the epochs
    @test losses[end] < 0.75 * losses[1]

    # one step of the trainer is one step of the reference: the loss of the first epoch agrees to
    # the round-off of two summations of the same terms in another order
    @test isapprox(losses[1], reference[1]; rtol = 1e-12, atol = 0)

    # every loss agrees with the reference to `1e-8`. The margin is wide: the largest relative
    # difference over the 1000 epochs is below `1e-15`. At ten times this step size the run is
    # unstable, and the two curves part at the round-off of either
    @test all(isapprox.(losses, reference; rtol = 1e-8, atol = 0))

    # the control: the same run at a step size that differs in the third digit, which the bound
    # above sees at the first epoch
    _, other = train_hnn(model, ps₀, data, target;
        η = 1.01e-4, epochs = epochs, batchsize = batchsize, rng = Xoshiro(0),
        batches = batches)
    @test !isapprox(other[1], reference[1]; rtol = 1e-8, atol = 0)
end

@testset "the edges of train_hnn" begin
    m, n, batchsize = 3, 20, 4
    model = hnn_model(m)
    rng = Xoshiro(3)
    data = randn(rng, m, n)
    target = randn(rng, m, n)
    ps₀ = params(NeuralNetwork(model; rng = Xoshiro(3)))
    batches = [rand(Xoshiro(11), 1:n, batchsize) for _ in 1:3]
    keywords = (;
        η = 0.1, epochs = 3, batchsize = batchsize, rng = Xoshiro(5), batches = batches)

    @testset "the shapes" begin
        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, target[:, 1:(end - 1)]; keywords...)
        end
        @test err isa DimensionMismatch
        @test occursin("$(size(data))", err.msg) &&
              occursin("$(size(target[:, 1:(end - 1)]))", err.msg)

        err = refusal(ps₀) do
            train_hnn(model, ps₀, randn(rng, m + 1, n), randn(rng, m + 1, n); keywords...)
        end
        @test err isa DimensionMismatch
        @test occursin("$(m + 1) reduced positions", err.msg) &&
              occursin("takes $m", err.msg)

        err = refusal(ps₀) do
            train_hnn(model, ps₀, randn(rng, m - 1, n), randn(rng, m - 1, n); keywords...)
        end
        @test err isa DimensionMismatch
        @test occursin("$(m - 1) reduced positions", err.msg) &&
              occursin("takes $m", err.msg)

        wide = Chain(Dense(m, m, tanh), Dense(m, 2))
        err = refusal(ps₀) do
            train_hnn(wide, ps₀, data, target; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("maps to 2 outputs", err.msg)
    end

    @testset "the layers" begin
        # the gradient is written for `tanh` and `identity`, so another activation is refused
        # with the layer that has it
        other = Chain(Dense(m, m, tanh), Dense(m, m, sin), Dense(m, 1, identity; use_bias = false))
        err = refusal(ps₀) do
            train_hnn(other, ps₀, data, target; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("layer 2 of the model", err.msg) &&
              occursin("`tanh` or `identity`", err.msg)
    end

    @testset "the sizes of a run" begin
        empty_data = zeros(m, 0)
        err = refusal(ps₀) do
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
            err = refusal(ps₀) do
                train_hnn(model, ps₀, data, target; keywords_...)
            end
            @test err isa ArgumentError
            @test occursin(phrase, err.msg)
        end
    end

    @testset "a batch larger than the samples" begin
        # the batches draw with replacement, so a batch may hold more indices than the data hold
        # samples: such a run trains for the epochs it is asked for, whether it draws its batches
        # or is given them
        _, losses = train_hnn(model, ps₀, data, target;
            η = keywords.η, epochs = keywords.epochs, batchsize = 30, rng = Xoshiro(5))
        @test length(losses) == keywords.epochs
        @test all(isfinite, losses)

        large = [rand(Xoshiro(6), 1:n, 30) for _ in 1:keywords.epochs]
        _, given = train_hnn(model, ps₀, data, target;
            η = keywords.η, epochs = keywords.epochs, batchsize = 30, rng = Xoshiro(5),
            batches = large)
        @test length(given) == keywords.epochs
        @test all(isfinite, given)
    end

    @testset "the batches" begin
        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, target; keywords..., batches = batches[1:2])
        end
        @test err isa DimensionMismatch
        @test occursin("hold 2 epochs", err.msg) && occursin("takes 3", err.msg)

        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, target;
                keywords..., batches = [
                    rand(Xoshiro(2), 1:n, batchsize + 1), batches[2:3]...])
        end
        @test err isa ArgumentError
        @test occursin("the batch of epoch 1 holds $(batchsize + 1) indices", err.msg)
        @test occursin("a batch holds $batchsize", err.msg)

        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, target;
                keywords..., batches = [batches[1], [1, 2, 3, n + 1], batches[3]])
        end
        @test err isa ArgumentError
        @test occursin("the batch of epoch 2", err.msg) && occursin("outside 1:$n", err.msg)

        # the indices of a batch start at one, so a zero is outside them too
        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, target;
                keywords..., batches = [batches[1], [0, 2, 3, 4], batches[3]])
        end
        @test err isa ArgumentError
        @test occursin("the batch of epoch 2", err.msg) && occursin("outside 1:$n", err.msg)
    end

    @testset "the element types" begin
        err = refusal(ps₀) do
            train_hnn(model, ps₀, Float32.(data), target; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("Float32", err.msg) && occursin("Float64", err.msg)

        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, Float32.(target); keywords...)
        end
        @test err isa ArgumentError
        @test occursin("the target Float32", err.msg)
    end

    @testset "a loss that is not finite" begin
        target_nan = copy(target)
        target_nan[1, 1] = NaN
        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, target_nan; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("the loss of epoch 1 is", err.msg) && occursin("NaN", err.msg)

        # an infinite loss is not finite either: the target that squares out of range is refused
        # with the epoch that found it
        target_inf = copy(target)
        target_inf[1, 1] = 1e300
        err = refusal(ps₀) do
            train_hnn(model, ps₀, data, target_inf; keywords...)
        end
        @test err isa ArgumentError
        @test occursin("the loss of epoch 1 is Inf", err.msg)
    end

    @testset "the parameters of the call are not the run's" begin
        before = deepcopy(ps₀)
        ps, _ = train_hnn(model, ps₀, data, target; keywords...)
        @test ps₀ == before
        @test ps != before
    end
end

@testset "the drawn batches and the other precision" begin
    m, n, epochs, batchsize = 3, 20, 3, 5
    model = hnn_model(m)
    rng = Xoshiro(9)
    data = randn(rng, m, n)
    target = randn(rng, m, n)
    ps₀ = params(NeuralNetwork(model; rng = Xoshiro(9)))
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
    ps₀32 = NetworkParameters((L1 = (W = Float32.(ps₀.L1.W), b = Float32.(ps₀.L1.b)),
        L2 = (W = Float32.(ps₀.L2.W), b = Float32.(ps₀.L2.b)), L3 = (W = Float32.(ps₀.L3.W),)))
    ps₃₂, losses32 = train_hnn(
        model, ps₀32, data32, target32; keywords..., rng = Xoshiro(7),
        batches = drawn)
    @test ps₃₂.L1.W isa Matrix{Float32} && losses32 isa Vector{Float32}
    @test isapprox(losses32, given_run[2]; rtol = 1e-4)
end

@testset "the allocations of a training run" begin
    # a run under coverage instruments every call, so the readings are taken only where it is off
    if Base.JLOptions().code_coverage == 0
        extension = Base.get_extension(VlasovMethods, :VlasovMethodsHNNExt)
        m, n, batchsize = 10, 200, 20
        model = hnn_model(m)
        rng = Xoshiro(4)
        data = randn(rng, m, n)
        target = randn(rng, m, n)
        ps = deepcopy(params(NeuralNetwork(model; rng = Xoshiro(4))))
        state = extension.TrainingState(model, ps, data, target, batchsize)
        opt = TrainingOptimizer(state.ps; algorithm = GradientMethod(), linesearch = 1e-8)
        indices = rand(Xoshiro(5), 1:n, batchsize)

        # fixed-arity barriers of concrete arguments, which call once before they measure
        function epoch_allocations(f::F, state, opt, indices, losses) where {F}
            f(state, opt, indices, 1, losses)
            return @allocated f(state, opt, indices, 1, losses)
        end
        function run_allocations(f::F, state, opt, batches, rng, losses) where {F}
            f(state, opt, batches, rng, losses)
            return @allocated f(state, opt, batches, rng, losses)
        end

        # an epoch writes into the buffers of the state and allocates nothing
        @test epoch_allocations(extension.epoch!, state, opt, indices, zeros(3)) == 0

        # a run of drawn batches allocates the buffer of their indices once: its bytes do not
        # grow with the number of epochs, up to the 16 bytes Windows rounds a reading to
        drawn(E) = run_allocations(
            extension._epochs!, state, opt, nothing, Xoshiro(6), zeros(E))
        given(E) = run_allocations(
            extension._epochs!, state, opt, fill(indices, E), Xoshiro(6), zeros(E))
        @test abs(drawn(300) - drawn(100)) ≤ 16
        @test given(100) == 0
    end
end

@testset "the docstring documents the method's arguments" begin
    # the call of the docstring names the method's own positional arguments, in their order, so that
    # an argument the extension takes and `?train_hnn` does not show is a failing test rather than a
    # stale page; the prose around the call is free. The docsystem holds the docstring as written
    # until something renders it, and as Markdown afterwards, and `string` reaches both.
    method = only(methods(VlasovMethods.train_hnn))
    documented = string(Base.@doc(VlasovMethods.train_hnn))
    call = "train_hnn(" * join(string.(Base.method_argnames(method)[2:end]), ", ")
    @test occursin(call, documented)
end

@testset "the two packages are weak dependencies" begin
    project = TOML.parsefile(joinpath(pkgdir(VlasovMethods), "Project.toml"))

    @test get(project["extensions"], "VlasovMethodsHNNExt", nothing) ==
          ["AbstractNeuralNetworks", "GeometricOptimizers"]
    for name in ("AbstractNeuralNetworks", "GeometricOptimizers")
        @test haskey(project["weakdeps"], name)
        @test haskey(project["compat"], name)
        @test !haskey(project["deps"], name)
    end

    # neither a symbolic engine nor a machine-learning framework is a dependency
    for name in ("Optimisers", "Zygote", "Lux", "Symbolics", "SymbolicNeuralNetworks",
        "NeuralNetworkParameters")
        @test !haskey(project["deps"], name)
        @test !haskey(project["weakdeps"], name)
    end
end

@testset "a fresh environment without the two packages" begin
    program = """
        import Pkg
        Pkg.activate(; temp = true, io = devnull)
        Pkg.develop(path = $(repr(pkgdir(VlasovMethods))); io = devnull)

        # a weak dependency of the package is not in the manifest of an environment that names
        # neither it nor a package that depends on it
        dependencies = Pkg.dependencies()
        for (uuid, name) in (
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
