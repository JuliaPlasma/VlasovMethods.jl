# The functions of `scripts/geometric_integrators_scaling.jl`, which compares the RK4 step of
# GeometricIntegrators with a hand-written RK4 loop on the grid 1D1V Vlasov–Poisson state. The
# script is included in a module of its own, which runs no measurement. `scaling_record` runs at
# the smallest size only, and no assertion here bounds a time: the measurement itself is the run
# of the script as a program.

using GeometricIntegrators
using LinearAlgebra: BLAS
using TOML
using Test
using VlasovMethods

module Scaling
include(joinpath(@__DIR__, "..", "..", "scripts", "geometric_integrators_scaling.jl"))
end

const S = Scaling

const PACKAGES = ["VlasovMethods", "GeometricIntegrators", "GeometricIntegratorsBase",
    "GeometricBase", "GeometricEquations", "GeometricBrackets", "PoissonSolvers",
    "SimpleSolvers", "RungeKutta"]

@testset "the verdict at and around the ratio limit" begin
    @test S.RATIO_LIMIT == 1.2
    z = zeros(3)
    @test S.scaling_verdict([1.0, 1.2, 1.2], z) == "accepted"
    @test S.scaling_verdict([1.19, 1.0, 1.1], z) == "accepted"
    @test S.scaling_verdict([1.0, 1.21, 1.0], z) == "rejected"
    @test S.scaling_verdict([1.3, 1.0, 1.0], z) == "rejected"
    @test S.scaling_verdict([1.0, 1.0, 1.3], z) == "rejected"
    @test S.scaling_verdict([1.0, NaN, 1.0], z) == "rejected"
    @test S.scaling_verdict([1.0, Inf, 1.0], z) == "rejected"
    @test S.scaling_verdict([1.21], [0.0]) == "rejected"
end

@testset "the verdict on the bytes per step across the sizes" begin
    r = ones(3)
    @test S.scaling_verdict(r, [0.0, 0.0, 16.0]) == "rejected"
    @test S.scaling_verdict(r, [0.0, 16.0, 0.0]) == "rejected"
    @test S.scaling_verdict(r, [100.0, 100.0, 100.0]) == "accepted"
    @test S.scaling_verdict(r, [100.0, 50.0, 50.0]) == "accepted"
    @test S.scaling_verdict(r, [100.0, 50.0, 80.0]) == "rejected"
    @test S.scaling_verdict([1.0], [1e6]) == "accepted"
end

@testset "the verdict refuses no sizes and inputs of different lengths" begin
    @test_throws r"^ArgumentError: scaling_verdict needs at least one size$" S.scaling_verdict(
        Float64[], Float64[])
    @test_throws r"^ArgumentError: scaling_verdict takes one ratio and one byte count per size" S.scaling_verdict(
        [1.0, 1.0], [0.0])
end

@testset "the verdict of a record reads the integrator's bytes" begin
    entry(r, ib, lb) = Dict("ratio" => r, "integrator_bytes_per_step" => ib,
        "loop_bytes_per_step" => lb)
    @test S.rk4_verdict([entry(1.0, 0.0, 0.0), entry(1.0, 16.0, 0.0)]) == "rejected"
    @test S.rk4_verdict([entry(1.0, 16.0, 0.0), entry(1.0, 16.0, 32.0)]) == "accepted"
    @test S.rk4_verdict([entry(1.0, 0.0, 0.0), entry(1.3, 0.0, 0.0)]) == "rejected"
end

@testset "the clock grows by dt at every step" begin
    model, dt = S.landau_model(16, 32)
    options = (warmup = 1, samples = 2, steps = 3)
    int = VlasovMethods.GeometricIntegrator(model, (0.0, 10dt), dt, RK4())
    solstep, t₀ = S.solution_step(int)
    S.measure_integrator(solstep, int.integrator, t₀, dt, options)
    @test solstep.t[] ≈ t₀ + 7dt rtol = 8eps()

    int = VlasovMethods.GeometricIntegrator(model, (0.0, 10dt), dt, RK4())
    solstep, t₀ = S.solution_step(int)
    q = copy(vec(model.distribution.values))
    work = (
        k₁ = similar(q), k₂ = similar(q), k₃ = similar(q), k₄ = similar(q), y = similar(q))
    S.measure_rk4(
        solstep, int.integrator, q, work, int.equation.parameters, t₀, dt, options)
    @test solstep.t[] ≈ t₀ + 7dt rtol = 8eps()
end

@testset "the implicit midpoint rule at its bound" begin
    J = 536_870_912
    @test S.jacobian_bytes(64, 128) == J
    @test S.implicit_fits(64, 128, 4J)
    @test !S.implicit_fits(64, 128, 4J - 1)
    @test S.implicit_fits(64, 128, UInt64(64) * 2^30)
    @test !S.implicit_fits(256, 512, UInt64(64) * 2^30)
    @test !S.implicit_fits(1024, 1024, UInt64(1024) * 2^30)
end

@testset "the per-step statistics" begin
    @test S.per_step([0.3, 0.2, 0.4], [160, 320, 160], 16) ==
          (seconds = 0.2 / 16, bytes = 20.0)
    @test S.per_step([0.5], [0], 2) == (seconds = 0.25, bytes = 0.0)
end

"Whether `f()` throws the `ErrorException` of `check_agreement`'s bound."
function refused(f)
    try
        f()
    catch e
        return e isa ErrorException &&
               startswith(e.msg, "the integrator's and the loop's RK4 states differ by ")
    end
    return false
end

@testset "the agreement of the two RK4 states" begin
    q = collect(range(1.0, 2.0; length = 101))
    @test S.check_agreement(q, q, 10) == 0.0
    near = copy(q)
    near[end] += 16 * eps()     # 8 eps relative to the largest entry, 2: the bound of one step
    @test S.check_agreement(near, q, 1) == 8 * eps()
    above = copy(q)
    above[end] += 18 * eps()    # 9 eps relative: above the bound of one step, within that of two
    @test refused(() -> S.check_agreement(above, q, 1))
    @test S.check_agreement(above, q, 2) == 9 * eps()
    far = copy(q)
    far[7] += 1e-6
    @test refused(() -> S.check_agreement(far, q, 10))
    nan = copy(q)
    nan[7] = NaN
    @test refused(() -> S.check_agreement(nan, q, 10))
end

@testset "the uptime text" begin
    @test S.uptime_text(nothing) == "unavailable"
    command = Sys.which("uptime")
    command === nothing || @test !isempty(S.uptime_text(command))
end

@testset "scaling_record refuses no sizes and sizes out of order" begin
    @test_throws r"^ArgumentError: scaling_record needs at least one size$" S.scaling_record([])
    order = r"^ArgumentError: scaling_record takes sizes in strictly increasing nₓ n_v"
    @test_throws order S.scaling_record([(256, 512), (64, 128)])
    @test_throws order S.scaling_record([(64, 128), (128, 64)])
end

@testset "the record at the smallest size" begin
    dir = mktempdir()
    record = S.scaling_record([(64, 128)];
        rk4 = (warmup = 1, samples = 2, steps = 2),
        midpoint = (warmup = 1, samples = 2, steps = 1))
    path = S.write_record(dir, record)
    @test path == joinpath(dir, "geometric_integrators_scaling.toml")
    R = TOML.parsefile(path)

    @test R["julia_version"] == string(VERSION)
    @test R["threads"] == Threads.nthreads()
    @test R["blas_threads"] == BLAS.get_num_threads()
    @test R["uptime"] isa String && !isempty(R["uptime"])
    @test length(R["load_average"]) == 3
    @test R["total_memory"] == Sys.total_memory()
    @test sort(collect(keys(R["packages"]))) == sort(PACKAGES)
    @test R["packages"]["VlasovMethods"] == string(pkgversion(VlasovMethods))
    @test R["packages"]["GeometricIntegrators"] == string(pkgversion(GeometricIntegrators))
    @test R["sizes"] == [[64, 128]]

    e = only(R["rk4"])
    @test (e["nx"], e["nv"]) == (64, 128)
    @test e["dt"] ≈ 0.5 * (4π / 64) / 8.6 rtol = 1e-12
    @test (e["warmup_steps"], e["samples"], e["steps_per_sample"]) == (1, 2, 2)
    for k in ("integrator_seconds_per_step", "loop_seconds_per_step")
        @test isfinite(e[k]) && e[k] > 0
    end
    for k in ("integrator_bytes_per_step", "loop_bytes_per_step")
        @test isfinite(e[k]) && e[k] ≥ 0
    end
    @test e["ratio"] ≈ e["integrator_seconds_per_step"] / e["loop_seconds_per_step"] rtol = 4eps()
    @test 0 ≤ e["state_difference"] ≤ 8 * 5 * eps()

    m = only(R["implicit_midpoint"])
    @test (m["nx"], m["nv"], m["jacobian_bytes"]) == (64, 128, 536_870_912)
    @test m["measured"] == S.implicit_fits(64, 128, Sys.total_memory())
    if m["measured"]
        @test isfinite(m["integrator_seconds_per_step"]) &&
              m["integrator_seconds_per_step"] > 0
        @test isfinite(m["integrator_bytes_per_step"]) && m["integrator_bytes_per_step"] ≥ 0
    end

    @test R["ratio_limit"] == 1.2
    @test R["verdict"] == S.scaling_verdict([e["ratio"]], [e["integrator_bytes_per_step"]])
end

@testset "an unmeasured implicit midpoint entry" begin
    record = S.scaling_record([(64, 128)];
        memory = 4 * S.jacobian_bytes(64, 128) - 1,
        rk4 = (warmup = 1, samples = 1, steps = 1))
    m = only(record["implicit_midpoint"])
    @test !m["measured"]
    @test m["jacobian_bytes"] == 536_870_912
    @test !haskey(m, "integrator_seconds_per_step")
    @test !haskey(m, "integrator_bytes_per_step")
    @test record["total_memory"] == 4 * 536_870_912 - 1
end

@testset "the record file is overwritten" begin
    dir = mktempdir()
    @test S.write_record(dir, Dict{String, Any}("a" => 1)) == joinpath(dir, S.RECORD)
    @test S.write_record(dir, Dict{String, Any}("a" => 2)) == joinpath(dir, S.RECORD)
    @test TOML.parsefile(joinpath(dir, S.RECORD)) == Dict{String, Any}("a" => 2)
end
