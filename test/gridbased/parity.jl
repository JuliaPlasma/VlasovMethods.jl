# The parity split and merge, the grading of the grid model's bracket term, and the
# construction, the runs and the output of `ParitySplitting`. The long runs, energy over 10⁴
# steps and the order, are in `test/verification/parity.jl`.
#
# `H = hx hv Σ f v²/2 + (hx/2) ρᵀ ϕ`, with `ϕ` the backend's field for `ρ = hv Σⱼ f`, is computed
# here and not in `src`, and `P(f, h) = poisson_apply(Arakawa(nx, nv, hx, hv), f, h) / (hx hv)` is
# the grid model's bracket term, written out from GeometricBrackets.

using GeometricBrackets: Arakawa, poisson_apply
using HDF5: h5open, h5read
using LinearAlgebra: dot, norm
using PoissonSolvers: FFTWBasis, Potential, solve
using Random: Xoshiro
using Test
using VlasovMethods

# standard-normal draws in Float64, rounded to T
draw(rng, T, dims...) = T.(randn(rng, dims...))

# the Maxwellian e^{-v²/2}/√(2π) at the nodes `v`
maxwellian(T, v) = exp.(-v .^ 2 ./ 2) ./ sqrt(T(2π))

# The grid model on x ∈ [0, 4π) with 16 nodes and v ∈ `vd` with `nv` nodes, at the values
# `shape(x, v)`, by default the perturbed Maxwellian (1 + 0.2 cos(x/2)) e^{-v²/2}/√(2π).
function grid_model(T, nv = 32; vd = (-T(7.5), T(7.5)),
        shape = (x, v) -> (1 + T(0.2) * cos(x / 2)) * exp(-v^2 / 2) / sqrt(T(2π)))
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, 16, nv), xd, T.(vd))
    dist.values .= shape.(dist.x, dist.v')
    return VlasovPoisson(dist, Potential(FFTWBasis(xd, 16)))
end

# the grid model's bracket term on the grid of `dist`
function bracket(dist)
    nx, nv = size(dist.values)
    hx, hv = step(dist.x), step(dist.v)
    a = Arakawa(nx, nv, hx, hv)
    return (f, h) -> reshape(poisson_apply(a, vec(f), vec(h)), nx, nv) ./ (hx * hv)
end

function energy(model, F)
    dist = model.distribution
    hx, hv = step(dist.x), step(dist.v)
    ρ = hv .* vec(sum(F; dims = 2))
    return hx * hv * sum(F .* (dist.v' .^ 2 ./ 2)) +
           (hx / 2) * dot(ρ, solve(model.potential.solver, ρ))
end

# `z`, `t` and the returned distribution of `run!(method, file)` on a new file
function run_file(method)
    mktempdir() do dir
        file = joinpath(dir, "run.h5")
        dist = run!(method, file)
        return h5read(file, "z"), h5read(file, "t"), dist
    end
end

# the exception that `f()` throws, or `nothing`
function thrown(f)
    try
        f()
    catch e
        return e
    end
    return nothing
end

@testset "parity: the parts have their parities, $T, nv = $nv" for T in (Float32, Float64),
    nv in (32, 33)

    rng = Xoshiro(1)
    f = draw(rng, T, 16, nv)
    g = draw(rng, T, nv)
    f0 = (g .+ reverse(g)) ./ 2
    fp, fm = parity_split(f, f0)
    @test eltype(fp) == eltype(fm) == T
    @test size(fp) == size(fm) == size(f)
    @test fm[:, end:-1:1] == -fm
    @test fp[:, end:-1:1] == fp
    @test norm(fm) > 0 && norm(fp) > 0
    isodd(nv) && @test all(iszero, fm[:, (nv + 1) ÷ 2])
end

@testset "round-trip: the merge inverts the split, $T, nv = $nv" for T in (Float32, Float64),
    nv in (32, 33)

    rng = Xoshiro(2)
    for _ in 1:5
        f = draw(rng, T, 16, nv)
        g = draw(rng, T, nv)
        f0 = (g .+ reverse(g)) ./ 2
        merged = parity_merge(f0, parity_split(f, f0)...)
        @test eltype(merged) == T
        @test maximum(abs, merged .- f) ≤ 4 * eps(T) * maximum(abs, f)
    end
end

@testset "parity_split: a background of the wrong length" for T in (Float32, Float64)
    f = zeros(T, 16, 32)
    guard = r"^DimensionMismatch: parity_split takes a background of size\(f, 2\) = 32 values"
    @test_throws guard parity_split(f, zeros(T, 31))
    @test_throws guard parity_split(f, zeros(T, 33))
    # a single value would broadcast over every column
    @test_throws guard parity_split(f, zeros(T, 1))
end

@testset "parity_merge: parts and a background of the wrong size" for T in (Float32, Float64)
    fp, fm, f0 = zeros(T, 16, 32), zeros(T, 16, 32), zeros(T, 32)
    guard = r"^DimensionMismatch: parity_merge takes two parts of one size and a background"
    # each of these but the last would broadcast without the size check
    @test_throws guard parity_merge(f0, fp, fm[:, 1:1])
    @test_throws guard parity_merge(f0, fp[1:1, :], fm)
    @test_throws guard parity_merge(f0[1:1], fp, fm)
    @test_throws guard parity_merge(f0[1:31], fp, fm)
end

@testset "parity_split: values at floatmax(T), $T" for T in (Float32, Float64)
    # (a ± b)/2 overflows for |a ± b| > floatmax(T); a/2 ± b/2 does not
    rng = Xoshiro(4)
    nv = 32
    f = draw(rng, T, 16, nv)
    g = draw(rng, T, nv)
    f0 = (g .+ reverse(g)) ./ 2
    # opposite signs at mirrored nodes reach f⁻, equal signs reach f⁺
    f[1, 3], f[1, nv - 2] = floatmax(T), -floatmax(T)
    f[2, 5], f[2, nv - 4] = floatmax(T), floatmax(T)
    fp, fm = parity_split(f, f0)
    @test eltype(fp) == eltype(fm) == T
    @test all(isfinite, fp) && all(isfinite, fm)
    merged = parity_merge(f0, fp, fm)
    @test eltype(merged) == T
    @test all(isfinite, merged)
    @test maximum(abs, merged .- f) ≤ 4 * eps(T) * maximum(abs, f)
    # the rows without a value at floatmax(T) to the bound of the round-trip testset
    rows = 3:16
    @test maximum(abs, merged[rows, :] .- f[rows, :]) ≤
          4 * eps(T) * maximum(abs, f[rows, :])
end

@testset "parity_split: a background of another element type" begin
    rng = Xoshiro(5)
    f = draw(rng, Float32, 16, 32)
    g = randn(rng, 32)
    f0 = (g .+ reverse(g)) ./ 2
    fp, fm = parity_split(f, f0)
    @test eltype(fp) == eltype(fm) == Float32
    merged = parity_merge(Float32.(f0), fp, fm)
    @test eltype(merged) == Float32
    @test isapprox(merged, f; rtol = 8 * eps(Float32))
end

@testset "grading: J(even, even) is odd, $T, nv = $nv" for T in (Float32, Float64),
    nv in (32, 33)

    dist = grid_model(T, nv).distribution
    P = bracket(dist)
    rng = Xoshiro(3)
    z = zeros(T, nv)
    a = parity_split(draw(rng, T, 16, nv), z)[1]
    b = parity_split(draw(rng, T, 16, nv), z)[1]
    c = parity_split(draw(rng, T, 16, nv), z)[2]
    Pab, Pcb = P(a, b), P(c, b)
    @test eltype(Pab) == eltype(Pcb) == T
    @test norm(Pab) > 0 && norm(Pcb) > 0
    # P(even, even) is odd: its even part vanishes; P(odd, even) is even: its odd part vanishes
    @test norm(parity_split(Pab, z)[1]) ≤ 64 * eps(T) * norm(Pab)
    @test norm(parity_split(Pcb, z)[2]) ≤ 64 * eps(T) * norm(Pcb)
end

@testset "construction: the checks of ParitySplitting, $T" for T in (Float32, Float64)
    model = grid_model(T)
    v = collect(model.distribution.v)
    M = maxwellian(T, v)

    # a v-grid that is not symmetric
    @test thrown(() -> ParitySplitting(grid_model(T; vd = (-5, 6)), (0, 1), 0.1)) isa
          ArgumentError
    # an end node 100 eps off its mirror is outside the bound of 10 eps
    vd = (-T(7.5), T(7.5) * (1 + 100 * eps(T)))
    @test thrown(() -> ParitySplitting(grid_model(T; vd), (0, 1), 0.1)) isa ArgumentError

    # a background that is not even, one that depends on x, one of the wrong length
    @test thrown(() -> ParitySplitting(model, (0, 1), 0.1;
        background = M .* (1 .+ T(0.1) .* v))) isa ArgumentError
    # the peak value 100 eps off its mirror is outside the bound of 10 eps of the largest value
    j = argmax(M)
    b = copy(M)
    b[j] *= 1 + 100 * eps(T)
    @test thrown(() -> ParitySplitting(model, (0, 1), 0.1; background = b)) isa
          ArgumentError
    # 2 eps off is inside the bound, and the stored background is symmetrised
    b = copy(M)
    b[j] *= 1 + 2 * eps(T)
    @test b != reverse(b)
    f0 = ParitySplitting(model, (0, 1), 0.1; background = b).background
    @test f0 == reverse(f0)
    @test eltype(f0) == T
    x = collect(model.distribution.x)
    e = thrown(() -> ParitySplitting(model, (0, 1), 0.1;
        background = (1 .+ T(0.2) .* cos.(x ./ 2)) .* M'))
    @test e isa ArgumentError && occursin("function of v", e.msg)
    @test thrown(() -> ParitySplitting(model, (0, 1), 0.1; background = M[1:(end - 1)])) isa
          DimensionMismatch
    @test thrown(() -> ParitySplitting(model, (0, 1), 0.1; background = 1)) isa
          ArgumentError

    # non-finite input, each with a message that names it
    for bad in (T(NaN), T(Inf))
        m = grid_model(T)
        m.distribution.values[3, 5] = bad
        e = thrown(() -> ParitySplitting(m, (0, 1), 0.1))
        @test e isa ArgumentError && occursin("finite", e.msg)
    end
    b = copy(M)
    b[4] = T(NaN)
    e = thrown(() -> ParitySplitting(model, (0, 1), 0.1; background = b))
    @test e isa ArgumentError && occursin("finite", e.msg)

    # the step, the time span, the tolerance and the iterations
    # 1e39 is finite in Float64 and overflows to Inf in Float32
    for τ in (T == Float32 ? (-0.1, 0, NaN, Inf, -Inf, 1e39) : (-0.1, 0, NaN, Inf, -Inf))
        e = thrown(() -> ParitySplitting(model, (0, 1), τ))
        @test e isa ArgumentError && occursin("tstep", e.msg)
    end
    for tspan in ((0, NaN), (Inf, 1), (0, 1, 2))
        e = thrown(() -> ParitySplitting(model, tspan, 0.1))
        @test e isa ArgumentError && occursin("tspan is not two finite times", e.msg)
    end
    # a span that runs backward is refused by a check of the constructor, which does not rely on
    # the time grid of `run!`; a span of equal endpoints is accepted
    e = thrown(() -> ParitySplitting(model, (1, 0), 0.1))
    @test e isa ArgumentError && occursin("tspan does not point forward", e.msg)
    @test ParitySplitting(model, (1, 1), 0.1).tspan === (T(1), T(1))
    # a zero background is even, at the edge of the evenness bound
    @test ParitySplitting(model, (0, 1), 0.1; background = zeros(T, 32)).background ==
          zeros(T, 32)
    # a step count that does not fit an `Int` throws at construction, in the time grid of `run!`
    @test thrown(() -> ParitySplitting(model, (0, 1), floatmin(T))) isa InexactError
    for tol in (-eps(T), NaN)
        e = thrown(() -> ParitySplitting(model, (0, 1), 0.1; tol))
        @test e isa ArgumentError && occursin("tol", e.msg)
    end
    # the solver stores the iterations as an Int32
    for max_iterations in (-1, Int(typemax(Int32)) + 1)
        e = thrown(() -> ParitySplitting(model, (0, 1), 0.1; max_iterations))
        @test e isa ArgumentError && occursin("max_iterations", e.msg)
    end

    # the edges that are accepted: a zero tolerance, no iteration, an odd nv, an even background
    @test ParitySplitting(model, (0, 1), 0.1; tol = 0) isa VlasovMethods.VlasovMethod
    @test ParitySplitting(model, (0, 1), 0.1; max_iterations = 0) isa ParitySplitting
    @test ParitySplitting(grid_model(T, 33), (0, 1), 0.1) isa ParitySplitting
    @test ParitySplitting(model, (0, 1), 0.1; background = M) isa ParitySplitting
    @test ParitySplitting(model, (0, 1), 0.1; background = Float64.(M)) isa ParitySplitting
end

@testset "construction: a background and values at floatmax(T), $T" for T in (Float32, Float64)
    # the symmetrisation and the default background halve before they add, so neither overflows
    model = grid_model(T)
    nx, nv = size(model.distribution.values)
    b = fill(floatmax(T), nv)
    f0 = ParitySplitting(model, (0, 1), 0.1; background = b).background
    @test eltype(f0) == T
    @test all(isfinite, f0)
    @test isapprox(f0, b; rtol = 8 * eps(T))

    # a column pair at floatmax(T) in one row, with the default background
    model.distribution.values[1, 3] = floatmax(T)
    model.distribution.values[1, nv - 2] = floatmax(T)
    f0 = ParitySplitting(model, (0, 1), 0.1).background
    @test eltype(f0) == T
    @test all(isfinite, f0)
    @test f0 == reverse(f0)
    @test isapprox(f0[3], floatmax(T) / nx; rtol = 8 * eps(T))
end

@testset "construction: the default background of a state with an odd part, $T" for T in (Float32,
    Float64)
    shape = (x, v) -> (1 + T(0.2) * cos(x / 2)) * exp(-v^2 / 2) / sqrt(T(2π)) *
                      (1 + T(0.3) * v)
    model = grid_model(T; shape)
    f = copy(model.distribution.values)
    nx = size(f, 1)
    z, = run_file(ParitySplitting(model, (0, 0.1), 0.1))
    f0 = vec(sum((f .+ f[:, end:-1:1]) ./ 2; dims = 1)) ./ nx
    model.distribution.values .= f
    zf0, = run_file(ParitySplitting(model, (0, 0.1), 0.1; background = f0))
    @test size(z) == (16 * 32, 2)
    @test eltype(z) == T
    @test maximum(abs, z .- zf0) ≤ 8 * eps(T) * maximum(abs, z)
    # the step changes the state, so the comparison is not the one of two initial states
    @test maximum(abs, z[:, 2] .- z[:, 1]) > 1000 * eps(T) * maximum(abs, z)
end

@testset "f⁺ solve: the JVP is the derivative of the residual, $T" for T in (Float32, Float64)
    # The residual is affine in y, so F(y + w) - F(y) is its derivative applied to w, to round-off.
    model = grid_model(T; shape = (x, v) -> (1 + T(0.2) * cos(x / 2)) *
                                            exp(-(v - 1)^2 / 2) / sqrt(T(2π)))
    nx, nv = size(model.distribution.values)
    τ = T(0.1)
    m = ParitySplitting(model, (0, 1), τ)
    kinetic = repeat(collect(model.distribution.v) .^ 2 ./ 2; inner = nx)
    F0 = repeat(m.background; inner = nx)
    fp, fm = vec.(parity_split(m.initial, m.background))
    fm = VlasovMethods._parity_fminus_step(m, kinetic, F0, fp, fm, τ)
    rng = Xoshiro(6)
    y, w = draw(rng, T, nx * nv), draw(rng, T, nx * nv)
    F(y) = VlasovMethods._parity_fplus_residual!(similar(y), y, m, kinetic, F0, fp, fm, τ)
    Jw = VlasovMethods._parity_fplus_jvp!(similar(w), w, m, fm, τ)
    @test eltype(Jw) == eltype(F(y)) == T
    # the P term of the derivative is not negligible against the bound
    @test norm(Jw .- w) > 1000 * 16 * eps(T) * (norm(F(y .+ w)) + norm(F(y)))
    @test norm(Jw .- (F(y .+ w) .- F(y))) ≤ 16 * eps(T) * (norm(F(y .+ w)) + norm(F(y)))
end

@testset "stationary: a homogeneous state, $T" for T in (Float32, Float64)
    model = grid_model(T; shape = (x, v) -> exp(-v^2 / 2) / sqrt(T(2π)))
    f = copy(model.distribution.values)
    z, = run_file(ParitySplitting(model, (0, 1), 0.1))
    @test size(z) == (16 * 32, 11)
    @test eltype(z) == T
    @test maximum(abs, z .- vec(f)) ≤ 64 * eps(T) * maximum(f)
end

@testset "odd nv: energy over 40 steps, nv = 33, $T" for T in (Float32, Float64)
    model = grid_model(T, 33)
    nx, nv = size(model.distribution.values)
    z, = run_file(ParitySplitting(model, (0, 2), 0.05))
    @test size(z) == (nx * nv, 41)
    @test eltype(z) == T
    H = [energy(model, reshape(z[:, n], nx, nv)) for n in axes(z, 2)]
    @test eltype(H) == T
    @test maximum(abs, H .- H[1]) / abs(H[1]) ≤ (T == Float64 ? 1e-12 : 1e-4)
end

@testset "failure: a solve that does not converge, $T" for T in (Float32, Float64)
    model = grid_model(T)
    f = copy(model.distribution.values)
    method = ParitySplitting(model, (0, 1), 0.05; max_iterations = 0)
    mktempdir() do dir
        file = joinpath(dir, "run.h5")
        e = thrown(() -> run!(method, file))
        @test e isa ErrorException
        @test occursin("MAXITERS", e.msg) && occursin("‖F‖", e.msg) &&
              occursin("iterations", e.msg)
        # the file is closed: it opens again, for reading and for writing
        @test h5open(h5 -> haskey(h5, "z"), file, "r")
        @test h5open(h5 -> true, file, "w")
    end
    @test model.distribution.values == f
end

@testset "failure: the step test of the solve is off, $T" for T in (Float32, Float64)
    # With tol = 0 the residual test is not met, so the solve takes every one of its
    # max_iterations steps; a step test would end it earlier, as STALLED.
    model = grid_model(T)
    method = ParitySplitting(model, (0, 0.05), 0.05; tol = 0, max_iterations = 10)
    e = mktempdir(dir -> thrown(() -> run!(method, joinpath(dir, "run.h5"))))
    @test e isa ErrorException
    @test occursin("MAXITERS", e.msg) && occursin("after 10 iterations", e.msg)
end

@testset "run: the time grid, the layout and the precision, $T" for T in (Float32, Float64)
    model = grid_model(T)
    f = copy(model.distribution.values)
    method = ParitySplitting(model, (0.0, 1.0), 0.3)
    z, t, dist = run_file(method)
    @test size(z) == (16 * 32, 4)
    @test t == collect(T(0):T(0.3):T(1))
    @test eltype(z) == eltype(t) == eltype(dist.values) == T
    @test dist === model.distribution
    @test maximum(abs, z[:, 1] .- vec(f)) ≤ 4 * eps(T) * maximum(abs, f)
    @test vec(dist.values) == z[:, end]
    @test maximum(abs, z[:, end] .- vec(f)) > 1000 * eps(T) * maximum(abs, f)
end

@testset "run: two runs give one result, $T" for T in (Float32, Float64)
    model = grid_model(T)
    method = ParitySplitting(model, (0, 0.5), 0.1)
    z1, = run_file(method)
    z2, = run_file(method)
    @test size(z1) == (16 * 32, 6)
    @test maximum(abs, z1 .- z2) ≤ 8 * eps(T) * maximum(abs, z1)
end

@testset "run: the residual bound scales with the state, $T" for T in (Float32, Float64)
    # On a state scaled by 10⁻⁴, an absolute bound of 4 eps(T) on the residual of the f⁺ solve is
    # loose, and the energy can drift past its bound; the bound 4 eps(T) ‖F0 + f⁺‖₂ scales with
    # the state.
    s = T(1e-4)
    model = grid_model(T; shape = (x, v) -> s * (1 + T(0.2) * cos(x / 2)) * exp(-v^2 / 2))
    nx, nv = size(model.distribution.values)
    z, = run_file(ParitySplitting(model, (0, 2), 0.05))
    @test size(z) == (nx * nv, 41)
    @test eltype(z) == T
    H = [energy(model, reshape(z[:, n], nx, nv)) for n in axes(z, 2)]
    @test maximum(abs, z[:, end] .- z[:, 1]) > 1000 * eps(T) * maximum(abs, z)
    @test maximum(abs, H .- H[1]) / abs(H[1]) ≤ (T == Float64 ? 1e-12 : 1e-4)
end

@testset "run: the tolerance reaches the solve, $T" for T in (Float32, Float64)
    # A tolerance of 1 is met by the start of each solve, so f⁺ does not move, and the state
    # differs from the converged run.
    model = grid_model(T)
    f = copy(model.distribution.values)
    z, = run_file(ParitySplitting(model, (0, 0.5), 0.1))
    model.distribution.values .= f
    zloose, = run_file(ParitySplitting(model, (0, 0.5), 0.1; tol = 1))
    @test maximum(abs, z .- zloose) > 1000 * eps(T) * maximum(abs, z)
end
