using LinearAlgebra
using Test
using VlasovMethods
using VlasovMethods: AnalyticHamiltonian

include("../helpers/analytic_vlasov.jl")

# The phrase of the exception a call raises, or a fixed word when it raises nothing. A phrase
# assertion is what tells two `ArgumentError`s of the same call apart.
function refusal(f)
    try
        f()
        return "no exception"
    catch thrown
        return sprint(showerror, thrown)
    end
end

# `A` with the entry `I` moved by `δ`: the difference `A′ - A` is `δ` at `I` and zero elsewhere.
function moved(A::AbstractArray{T}, I::CartesianIndex, δ) where {T}
    B = copy(A)
    B[I] += T(δ)
    return B
end

@testset "check 3: the trajectory errors, $T" for T in (Float64, Float32)
    nₚ, nₛ, nparam = 7, 5, 3
    δ = T(2)

    # one particle differs, so `x2` is `abs(δ)/√nₚ` and `xmax` is `abs(δ)` there and zero elsewhere
    X = zeros(T, nₚ, nₛ, nparam)
    Xh = moved(X, CartesianIndex(3, 2, 2), δ)
    Vh = moved(X, CartesianIndex(5, 1, 3), δ)
    err = trajectory_errors(Xh, Vh, X, X)

    @test size(err.x2) == (nₛ, nparam)
    @test size(err.xmax) == (nₛ, nparam)
    @test eltype(err.x2) === T && eltype(err.xmax) === T
    @test eltype(err.v2) === T && eltype(err.vmax) === T
    @test err.x2[2, 2] == abs(δ) / sqrt(T(nₚ))
    @test err.xmax[2, 2] == abs(δ)
    @test err.v2[1, 3] == abs(δ) / sqrt(T(nₚ))
    @test err.vmax[1, 3] == abs(δ)
    @test count(!iszero, err.x2) == 1
    @test count(!iszero, err.xmax) == 1
    @test count(!iszero, err.v2) == 1
    @test count(!iszero, err.vmax) == 1

    # `xmax` is a maximum and not a 2-norm: two particles that differ give a norm above the
    # largest difference, which the one-particle fixture cannot tell apart
    Y = zeros(T, nₚ, 1, 1)
    Y[1, 1, 1] = one(T)
    Y[2, 1, 1] = one(T)
    both = trajectory_errors(Y, Y, zeros(T, nₚ, 1, 1), zeros(T, nₚ, 1, 1))
    @test both.xmax[1, 1] == one(T)
    @test both.x2[1, 1] == sqrt(T(2)) / sqrt(T(nₚ))

    # a `(nₚ, nₛ)` matrix gives a vector of length `nₛ`
    M = zeros(T, nₚ, nₛ)
    flat = trajectory_errors(moved(M, CartesianIndex(2, 4), δ), M, M, M)
    @test flat.x2 isa Vector{T}
    @test length(flat.x2) == nₛ
    @test flat.x2[4] == abs(δ) / sqrt(T(nₚ))
    @test flat.xmax[4] == abs(δ)

    # the two sizes are named when the HFM and the ROM arrays disagree
    @test_throws DimensionMismatch trajectory_errors(
        Xh, Vh, zeros(T, nₚ + 1, nₛ, nparam), X)
    @test_throws DimensionMismatch trajectory_errors(
        Xh, Vh, X, zeros(T, nₚ, nₛ, nparam + 1))
    positions = refusal(() -> trajectory_errors(
        Xh, Vh, zeros(T, nₚ + 1, nₛ, nparam), X))
    @test occursin("(7, 5, 3)", positions) && occursin("(8, 5, 3)", positions)
    velocities = refusal(() -> trajectory_errors(
        Xh, Vh, X, zeros(T, nₚ, nₛ, nparam + 1)))
    @test occursin("(7, 5, 3)", velocities) && occursin("(7, 5, 4)", velocities)
end

@testset "check 3b: x2 and v2 are root mean squares, independent of nₚ, $T" for T in (Float64, Float32)
    # `m` copies of one particle set, each particle with the same error: the root mean square is
    # unchanged, where the root sum of squares grows with `√m`
    Xh = T[1, 2, 3, 4]
    Xr = T[1.5, 2, 2, 4.5]
    Vh = T[0, 1, 0, 1]
    Vr = T[1, 1, -1, 1]
    base = trajectory_errors(reshape(Xh, :, 1), reshape(Vh, :, 1),
        reshape(Xr, :, 1), reshape(Vr, :, 1))
    # the squared errors are 1/4, 0, 1, 1/4 and 1, 0, 1, 0, with means 3/8 and 1/2
    @test base.x2[1] ≈ sqrt(T(3) / 8)
    @test base.v2[1] ≈ sqrt(T(1) / 2)

    for m in (2, 5)
        copies = trajectory_errors(
            reshape(repeat(Xh, m), :, 1), reshape(repeat(Vh, m), :, 1),
            reshape(repeat(Xr, m), :, 1), reshape(repeat(Vr, m), :, 1))
        @test copies.x2[1] ≈ base.x2[1]
        @test copies.v2[1] ≈ base.v2[1]
        @test copies.xmax[1] == base.xmax[1]
        @test copies.vmax[1] == base.vmax[1]
    end
end

@testset "check 4: the saved indices" begin
    @test saved_indices(10001) == (10001, 7501, 5001)
    @test saved_indices(4) == (4, 3, 2)
    @test saved_indices(5) == (5, 4, 3)   # `div(12, 4) + 1` and `div(4, 2) + 1`
    @test saved_indices(3) == (3, 2, 2)   # below four saved times two of the states coincide
    @test saved_indices(2) == (2, 1, 1)
    @test_throws ArgumentError saved_indices(1)
    @test_throws ArgumentError saved_indices(0)
end

@testset "check 5: the potential error, $T" for T in (Float64, Float32)
    nh, nₛ, nparam = 6, 4, 2
    Φ = reshape(collect(T, 1:(nh * nₛ * nparam)), nh, nₛ, nparam)

    # the same array gives zeros, of the shape of the trailing axes and the element type
    Z = potential_error(Φ, Φ)
    @test Z == zeros(T, nₛ, nparam)
    @test eltype(Z) === T

    # one coefficient differs: the absolute 2-norm, with no divisor
    δ = T(2)
    Φr = moved(Φ, CartesianIndex(3, 2, 1), -δ)
    E = potential_error(Φ, Φr)
    @test size(E) == (nₛ, nparam)
    @test E[2, 1] == abs(δ)

    # two coefficients that differ give the root sum of squares, which a single-point fixture
    # cannot tell from a maximum
    Φr2 = moved(moved(Φ, CartesianIndex(1, 1, 1), -one(T)), CartesianIndex(4, 1, 1), -one(T))
    @test potential_error(Φ, Φr2)[1, 1] == sqrt(T(2))

    @test_throws DimensionMismatch potential_error(Φ, zeros(T, nh + 1, nₛ, nparam))
    sizes = refusal(() -> potential_error(Φ, zeros(T, nh + 1, nₛ, nparam)))
    @test occursin("(6, 4, 2)", sizes) && occursin("(7, 4, 2)", sizes)
end

@testset "check 6: the field-energy error, $T" for T in (Float64, Float32)
    nₛ, nparam = 4, 2
    Wh = T[1, 2, 4, 8, 1, 2, 4, 8]
    Wr = T[1, 4, 3, 8, 2, 2, 8, 4]
    E = field_energy_error(reshape(Wh, nₛ, nparam), reshape(Wr, nₛ, nparam))

    @test eltype(E) === T
    @test E ==
          abs.(reshape(Wh, nₛ, nparam) .- reshape(Wr, nₛ, nparam)) ./
          reshape(Wh, nₛ, nparam)
    @test E[2, 1] == one(T)   # `abs(2 - 4)/2`

    # a `Wh` of zero propagates: `abs(0 - 0)/0` is a `NaN`
    @test isnan(field_energy_error(reshape(T[0, 1], 2, 1), reshape(T[0, 1], 2, 1))[1, 1])
    @test isinf(field_energy_error(reshape(T[0, 1], 2, 1), reshape(T[1, 1], 2, 1))[1, 1])
    @test_throws DimensionMismatch field_energy_error(zeros(T, 3), zeros(T, 4))
end

@testset "check 7: the energy drift, $T" for T in (Float64, Float32)
    @test energy_drift(T[1, 2, 4]) == T[0, 1, 3]
    @test energy_drift(T[5]) == T[0]
    @test eltype(energy_drift(T[1, 2])) === T

    # the exact analytic solution: `H` is constant, so only `sin` and `cos` rounding remains
    β = T(ANALYTIC_β)
    q0, p0 = analytic_initial_condition(T; n = 200)
    ts = collect(T, T(0):T(10):T(1000))
    Q, P = analytic_solution(T, β, q0, p0, ts)
    h = AnalyticHamiltonian(β)
    H = [h(zero(T), vcat(Q[:, j], P[:, j]), NamedTuple()) for j in eachindex(ts)]
    bound = T === Float64 ? 1e-13 : 1e-5
    @test maximum(abs, energy_drift(H)) ≤ bound

    @test_throws ArgumentError energy_drift(T[])
    @test occursin("the initial energy is zero", refusal(() -> energy_drift(T[0, 1])))
end

@testset "check 8: the maxima of a history, $T" for T in (Float64, Float32)
    # a maximum is a point no lower than the `window` points on each side
    @test VlasovMethods._find_maxima(T[0, 1, 2, 2, 1, 0], 1) == [3]
    # a flat maximum counts once, at its first index: a second equal sample is no maximum of its
    # own, and a fit that counts both weights the plateau twice
    @test VlasovMethods._find_maxima(T[0, 2, 2, 2, 0], 1) == [2]
    # two equal peaks that are apart are two maxima: only a point equal to the sample on its own
    # left is a repeat, so an equal value elsewhere in the window does not drop a peak
    @test VlasovMethods._find_maxima(T[0, 0, 3, 0, 3, 0, 0, 0], 2) == [3, 5]
    # the fixture of the decided edge, whose three maxima the refusals below name
    @test VlasovMethods._find_maxima(T[1, 2, 1, 3, 1, 4, 1], 1) == [2, 4, 6]

    # a rising shoulder is no maximum: the first sample that differs from it on its right is
    # higher, so the history is still climbing. The shoulder of `1` before the peak of `2` is
    # no maximum of a window of one, and neither is the longer shoulder of a window of two, whose
    # equal neighbours alone do not rule it out
    @test VlasovMethods._find_maxima(T[0, 1, 1, 2, 0], 1) == [4]
    @test VlasovMethods._find_maxima(T[0, 0, 1, 1, 1, 2, 0, 0], 2) == [6]
    # a plateau that runs to the end of the history never falls, so it is no maximum either
    @test isempty(VlasovMethods._find_maxima(T[0, 1, 2, 2, 2], 1))

    # a plateau that falls after several equal samples is one maximum at its first index, and the
    # fall may come after more equal samples than the window holds
    @test VlasovMethods._find_maxima(T[0, 2, 2, 2, 2, 1, 0], 1) == [2]
    @test VlasovMethods._find_maxima(T[0, 0, 2, 2, 2, 2, 1, 0, 0], 2) == [3]

    # a shoulder before a peak does not enter the fit: the rate of the history with the shoulder
    # is that of the same history without it
    ts = collect(T, 0:7)
    @test growth_rate(ts, T[1, 5, 1, 2, 2, 8, 1, 1]; window = 1, maxima = 1:2) ==
          growth_rate(ts, T[1, 5, 1, 2, 3, 8, 1, 1]; window = 1, maxima = 1:2)

    # one maximum per plateau: a plateau and its first sample alone give the same rate, so the
    # repeated sample of the plateau does not enter the fit twice
    tp = collect(T, 0:5)
    @test growth_rate(tp, T[1, 5, 5, 1, 3, 1]; window = 1, maxima = 1:2) ==
          growth_rate(tp, T[1, 5, 1, 1, 3, 1]; window = 1, maxima = 1:2)

    # two equal peaks that are apart carry the same height, so the fit through them is flat
    @test growth_rate(collect(T, 0:6), T[0, 0, 2, 0, 2, 0, 0];
        window = 2, maxima = 1:2) == zero(T)
end

@testset "check 8: the growth rate, $T" for T in (Float64, Float32)
    γ = T(0.05)
    t = collect(T, T(0):T(1e-3):T(100))
    W = exp.(2γ .* t) .* cos.(t) .^ 2

    g = growth_rate(t, W)
    @test g isa T
    @test isapprox(g, γ; rtol = 1e-4)

    # a matrix gives one rate per column, of the element type of `W`; the two columns hold
    # different rates, so a method that fits the first column for every column fails
    γ2 = T(0.03)
    W2 = exp.(2γ2 .* t) .* cos.(t) .^ 2
    gm = growth_rate(t, hcat(W, W2))
    @test gm isa Vector{T}
    @test length(gm) == 2
    @test isapprox(gm[1], γ; rtol = 1e-4)
    @test isapprox(gm[2], γ2; rtol = 1e-4)

    # `window = 1` finds the maxima of the sampled `cos²`: without it the flat maximum `Float32`
    # reaches at `h = 1e-3` is near enough to the true one for the fit to stay in its bound
    w1 = growth_rate(t, W; window = 1)
    @test isapprox(w1, γ; rtol = 1e-4)

    # `maxima` selects the maxima by their order along `t`, so it drops the first few: the
    # height of a dropped maximum does not enter the fit
    tm = collect(T, 0:11)
    Wm = T[1, 3, 1, 5, 1, 2, 1, 4, 1, 6, 1, 2]
    higher = copy(Wm)
    higher[2] = T(7)
    @test growth_rate(tm, Wm; window = 1, maxima = 2:5) ==
          growth_rate(tm, higher; window = 1, maxima = 2:5)

    # a `W` of zero or below at a selected maximum is refused, where `log` is not defined;
    # `window = 1` puts maxima at 2 and 5, and `W[2]` is the first of them
    for bad in (T(-2), zero(T))
        @test occursin("not positive",
            refusal(() -> growth_rate(
                T[0, 1, 2, 3, 4, 5, 6], T[-3, bad, -3, -1, 5, -1, 5]; window = 1,
                maxima = 1:2)))
    end
    @test occursin("needs two",
        refusal(() -> growth_rate(
            T[0, 1, 2, 3, 4, 5, 6], T[-3, -2, -3, -1, 5, -1, 5]; window = 1, maxima = 1:1)))
    @test occursin("beyond the 2 maxima found",
        refusal(() -> growth_rate(
            T[0, 1, 2, 3, 4, 5, 6], T[-3, -2, -3, -1, 5, -1, 5]; window = 1)))

    # the refusals of the decided edge, on its own fixture: `window = 1` finds three maxima at
    # 2, 4 and 6, an index may not enter the fit twice, and both count refusals name the count
    # found rather than the count selected
    a = collect(T, 0:6)
    Wa = T[1, 2, 1, 3, 1, 4, 1]
    for repeated in ([1, 2, 2], [2, 2])
        @test occursin("repeats the index 2",
            refusal(() -> growth_rate(a, Wa; window = 1, maxima = repeated)))
    end
    @test occursin("3 maxima found",
        refusal(() -> growth_rate(a, Wa; window = 1, maxima = 2:2)))

    # an index below 1, a window below 1 and times of another length are refused by name, not
    # by a `BoundsError` from inside the fit
    @test occursin("counted from 1",
        refusal(() -> growth_rate(a, Wa; window = 1, maxima = 0:2)))
    for window in (0, -1)
        @test occursin("window of 1 or more",
            refusal(() -> growth_rate(a, Wa; window = window, maxima = 1:2)))
    end
    @test_throws DimensionMismatch growth_rate(a[1:6], Wa; window = 1, maxima = 1:2)
    @test occursin("hold 6 values",
        refusal(() -> growth_rate(a[1:6], Wa; window = 1, maxima = 1:2)))
end

@testset "check 8: the best and worst parameters, $T" for T in (Float64, Float32)
    γh = T[1, 2, 3, 4, 5]
    γr = T[1, 5, 3, 2, 5]
    d = abs.(γr .- γh)
    best, worst = best_and_worst(γh, γr)

    @test d == T[0, 3, 0, 2, 0]
    @test best == argmin(d) == 1    # the first index of a tie
    @test worst == argmax(d) == 2
    @test best != worst
end

@testset "check 8: the growth_rate docstring states the rate it returns" begin
    # The rate is half the slope of `log W`, and the docstring says so without naming a source
    # outside the package
    doc = lowercase(string(@doc growth_rate))
    @test occursin("half the least-squares", doc)
    @test !occursin("notebook", doc)
end
