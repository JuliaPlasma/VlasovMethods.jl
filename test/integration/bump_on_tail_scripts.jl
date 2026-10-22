# The functions of the bump-on-tail drivers in `scripts/`, which `scripts/bump_on_tail_common.jl`
# holds so that this file can include it: the field of a march, the sweeps and the step count, the
# growth-rate fit of the reference, and the writer of the two files that the HNN reads. It also
# reads the drivers and the docs as text, for what no function of theirs can show: the one place a
# scaled field is built, the output paths, and the docs page. No driver runs here; a driver writes
# its figures with CairoMakie, which the common file does not load.

using HDF5
using LinearAlgebra
using PoissonSolvers
using Random
using ReducedBasisMethods: ReducedBasis, CotangentLiftEVD, PODEVD
using ReducedComplexityModeling: Parameter, ParameterSpace
using Test
using VlasovMethods
using VlasovMethods: PoissonField, ScaledField, efield!, update!

const REPOSITORY = joinpath(@__DIR__, "..", "..")
const SCRIPTS = joinpath(REPOSITORY, "scripts")

# The include loads no Makie package: every module it loads is one this environment had loaded
# before or one the common file names, and none of them is a Makie.
loaded_before = Set(keys(Base.loaded_modules))
include(joinpath(SCRIPTS, "bump_on_tail_common.jl"))
loaded_by_include = setdiff(Set(keys(Base.loaded_modules)), loaded_before)

"The text of every Julia file in `scripts/`, by name."
function script_texts()
    Dict(f => read(joinpath(SCRIPTS, f), String)
    for f in readdir(SCRIPTS) if endswith(f, ".jl"))
end

@testset "the common file loads no Makie" begin
    @test !any(id -> occursin("Makie", id.name), loaded_by_include)
    @test !occursin(r"^\s*(using|import)\b.*Makie"m,
        read(joinpath(SCRIPTS, "bump_on_tail_common.jl"), String))
    for f in ("bump_on_tail_training.jl", "bump_on_tail_projections.jl",
        "bump_on_tail_testing.jl", "analytic_quadratic.jl")
        @test isfile(joinpath(SCRIPTS, f))
    end
end

@testset "vp_field" begin
    potential = Potential(PeriodicBasisSpline((0.0, 2π / KAPPA_REFERENCE), 4, 16))
    Random.seed!(1234)
    x = (2π / KAPPA_REFERENCE) .* rand(200)
    w = fill(1 / 200, 200)

    field = vp_field(potential, (χ = 2.0,))
    @test field isa ScaledField
    @test field.field isa PoissonField
    @test field.field.potential === potential
    @test field.χ == 2.0

    # the field of the march at χ = 2 is one quarter of the plain field, exactly
    update!(field, x, w, 0.0)
    plain = PoissonField(potential)
    @test efield!(field, zeros(200), x) == efield!(plain, zeros(200), x) ./ 4
    @test any(!iszero, efield!(plain, zeros(200), x))

    # the χ is the one of the parameters
    @test vp_field(potential, (χ = 0.5, κ = 0.3)).χ == 0.5

    # one construction of a scaled field in the drivers, in the common file
    counts = Dict(f => count("ScaledField(", text) for (f, text) in script_texts())
    @test sum(values(counts)) == 1
    @test counts["bump_on_tail_common.jl"] == 1
end

@testset "accelerations" begin
    # the accelerations of a saved state are the field of the march at that state, `E/χ²`
    potential = Potential(PeriodicBasisSpline((0.0, 2π / KAPPA_REFERENCE), 4, 16))
    Random.seed!(7)
    n, ns = 50, 3
    space = ParameterSpace(Parameter(:χ, [0.5, 2.0]))
    SS = Snapshots(Float64, 1, n, 16, ns, space)
    SS.X .= (2π / KAPPA_REFERENCE) .* rand(1, n, ns, 2)
    w = fill(1 / n, n)
    accelerations!(SS, potential, w)

    for p in 1:2, ts in 1:ns

        x = SS.X[1, :, ts, p]
        plain = PoissonField(potential)
        update!(plain, x, w, 0.0)
        @test SS.A[1, :, ts, p] ≈ efield!(plain, zeros(n), x) ./ space(p).χ^2
    end
    @test SS.A[1, :, 1, 1] != SS.A[1, :, 1, 2]
end

@testset "paths" begin
    # every file goes through `runs_path` or `results_path`: no driver names the old output paths
    for (f, text) in script_texts(), line in eachline(IOBuffer(text))

        @test !occursin("../runs", line)
        @test !occursin("/toks/", line)
    end
end

@testset "sweeps" begin
    # the step count rounds: `Int(div(25.0, 0.1))` is `249`
    @test step_count(25.0, 0.1) == 250
    @test step_count(1000.0, 0.001) == 10^6

    # the training sweep: the ten κ of the old training script, evenly spaced in [0.1, 0.5]
    κ = training_kappas()
    @test length(κ) == 10
    @test first(κ) == 0.1
    @test last(κ) ≈ 0.5
    @test all(diff(κ) .≈ 0.4 / 9)

    # the test sweep: ten κ drawn uniform in [0.1, 0.5] after the seed, sorted, with χ = κ/0.3
    κ₁, χ₁ = test_sweep(1234)
    @test length(κ₁) == 10
    @test issorted(κ₁)
    @test all(0.1 .≤ κ₁ .≤ 0.5)
    @test χ₁ == κ₁ ./ 0.3
    @test test_sweep(1234) == (κ₁, χ₁)
    κ₂, _ = test_sweep(1)
    @test κ₂ != κ₁
    @test test_sweep(1234; n = 3)[1] |> length == 3

    # the parameter space of a sweep holds `χ = κ/κ₀`, one sample per κ, in the order of a sorted κ
    space = parameter_space([0.15, 0.3, 0.6])
    @test length(space) == 3
    @test [space(p).χ for p in 1:3] ≈ [0.5, 1.0, 2.0]
end

@testset "the domain and the draw" begin
    # the potential: `cells` cells of degree `degree` on `[0, 2π/κ₀]`
    b = VlasovMethods.basis(potential(BUMP_ON_TAIL))
    @test degree(b) == 3
    @test ncells(b) == 16
    @test first(breakpoints(b)) == 0
    @test last(breakpoints(b)) ≈ 2π / 0.3

    # the particles: drawn after the seed of the options, in the domain, with the total weight of
    # the domain's length
    o = merge(BUMP_ON_TAIL, (; particles = 500))
    P = initial_particles(o)
    @test length(P) == 500
    @test all(0 .≤ vec(P.x) .≤ 2π / 0.3)
    @test sum(P.w) ≈ 2π / 0.3
    rand(100)
    @test initial_particles(o).x == P.x
    @test initial_particles(merge(o, (; seed = 1))).x != P.x
end

@testset "growth rate" begin
    # `W` grows at γ = 0.05 up to its fourth maximum, the first near t = 0, and is flat after it
    W(t) = exp(0.1 * min(t, 3π)) * cos(t)^2
    t = 0:0.001:20
    @test isapprox(fit_growth_rate(t, W.(t)), 0.05; rtol = 1e-2)
    # the defaults of `growth_rate` fit the second to fifth maxima, which reach the flat part
    @test growth_rate(t, W.(t)) < 0.045

    # The fit takes the fourth maximum: flat after the third, at 2π, the line through the four
    # points `(kπ, 0.1·min(kπ, 2π))`, k = 0…3, has the slope `0.07`, where three give `0.1`.
    V(t) = exp(0.1 * min(t, 2π)) * cos(t)^2
    @test isapprox(fit_growth_rate(t, V.(t)), 0.035; rtol = 2e-2)

    # three maxima at window 5: the fit takes all three
    t3 = 0:0.001:(2.5π)
    @test isapprox(fit_growth_rate(t3, W.(t3)), 0.05; rtol = 1e-2)

    # A one-sample spike on a rising slope, above the two samples after it and below the third, is
    # a maximum at a window of 2 and not at 5: the count is at the window of the fit.
    w = W.(t3)
    i = findfirst(≥(2.5), t3)
    w[i] += 2.5 * (w[i + 1] - w[i])
    @test length(VlasovMethods._find_maxima(w, 2)) == 4
    @test isapprox(fit_growth_rate(t3, w), 0.05; rtol = 1e-2)

    # two maxima, near 0 and near π, are the fewest that the fit takes
    t2 = 0:0.001:(1.5π)
    @test isapprox(fit_growth_rate(t2, W.(t2)), 0.05; rtol = 1e-2)

    # one maximum has no fit, and the rate is `NaN`
    t1 = 0:0.001:2.0
    @test isnan(fit_growth_rate(t1, W.(t1)))

    # a matrix gives one rate per column
    @test fit_growth_rate(t3, [W.(t3) W.(t3)]) ≈ fill(fit_growth_rate(t3, W.(t3)), 2)
end

@testset "p39 files" begin
    # a two-sample fixture: a cotangent-lift basis of `k = 3` and the reduced data of `nₛ = 4`
    # saves of two samples
    Random.seed!(3)
    np, ns, k = 20, 4, 3
    space = ParameterSpace(Parameter(:χ, [0.5, 1.5]))
    SS = Snapshots(Float64, 1, np, 16, ns, space)
    SS.X .= rand(1, np, ns, 2)
    SS.V .= randn(1, np, ns, 2)
    SS.A .= randn(1, np, ns, 2)
    basis = particle_basis(CotangentLiftEVD(), SS; k = k)
    X, A = reduced_data(basis, SS)
    @test size(X) == size(A) == (k, ns, 2)
    @test X[:, 2, 1] ≈ basis.Φ' * SS.X[1, :, 2, 1]
    @test A[:, 3, 2] ≈ basis.Φ' * SS.A[1, :, 3, 2]

    commit = "0123456789abcdef0123456789abcdef01234567"
    mktempdir() do dir
        projections = joinpath(dir, "fixture_projections.h5")
        vectorfield = joinpath(dir, "fixture_vectorfield.h5")
        write_p39_files(projections, vectorfield, basis, X, A;
            commit = commit, command = "--seed 7", seed = 7)

        for f in (projections, vectorfield)
            h5open(f) do io
                a = attributes(io)
                @test read(a["commit"]) == commit
                @test read(a["command"]) == "--seed 7"
                @test read(a["seed"]) == 7
            end
        end
        @test ReducedBasis(projections).Φ == basis.Φ
        h5open(io -> @test(!haskey(io, "field")), projections)
        h5open(vectorfield) do io
            @test read(io["X"]) == X
            @test read(io["A"]) == A
        end

        # the field basis, when one is given, is in the group `field` of the projections file
        fields = field_basis(PODEVD(), SS; k = 2)
        write_p39_files(projections, vectorfield, basis, X, A;
            commit = commit, command = "", seed = 7, field_basis = fields)
        @test ReducedBasis(projections).Φ == basis.Φ
        @test ReducedBasis(projections, "field").Φ == fields.Φ
    end

    # the run id comes from the parameters
    o = merge(BUMP_ON_TAIL, (; seed = 1234))
    @test run_id(o) != run_id(merge(o, (; particles = 1000)))
    @test run_id(o) != run_id(merge(o, (; seed = 1)))
    @test occursin(string(o.particles), run_id(o))
end

@testset "docs page" begin
    make = read(joinpath(REPOSITORY, "docs", "make.jl"), String)
    @test occursin("\"paper.md\"", make)

    page = read(joinpath(REPOSITORY, "docs", "src", "paper.md"), String)
    for f in ("bump_on_tail_training.jl", "bump_on_tail_projections.jl",
        "bump_on_tail_testing.jl", "analytic_quadratic.jl")
        @test occursin(f, page)
    end
    @test !occursin("@example", page)
    @test !occursin("@eval", page)
    @test !occursin("![", page)
end
