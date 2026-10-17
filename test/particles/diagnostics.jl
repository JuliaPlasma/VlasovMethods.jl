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

    # one particle differs, so `x2` is `abs(δ)/nₚ` and `xmax` is `abs(δ)` there and zero elsewhere
    X = zeros(T, nₚ, nₛ, nparam)
    Xh = moved(X, CartesianIndex(3, 2, 2), δ)
    Vh = moved(X, CartesianIndex(5, 1, 3), δ)
    err = trajectory_errors(Xh, Vh, X, X)

    @test size(err.x2) == (nₛ, nparam)
    @test size(err.xmax) == (nₛ, nparam)
    @test eltype(err.x2) === T && eltype(err.xmax) === T
    @test eltype(err.v2) === T && eltype(err.vmax) === T
    @test err.x2[2, 2] == abs(δ) / nₚ
    @test err.xmax[2, 2] == abs(δ)
    @test err.v2[1, 3] == abs(δ) / nₚ
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
    @test both.x2[1, 1] == sqrt(T(2)) / nₚ

    # a `(nₚ, nₛ)` matrix gives a vector of length `nₛ`
    M = zeros(T, nₚ, nₛ)
    flat = trajectory_errors(moved(M, CartesianIndex(2, 4), δ), M, M, M)
    @test flat.x2 isa Vector{T}
    @test length(flat.x2) == nₛ
    @test flat.x2[4] == abs(δ) / nₚ
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

@testset "check 4: the saved indices" begin
    @test saved_indices(10001) == (10001, 7501, 5001)
    @test saved_indices(4) == (4, 3, 2)
    @test saved_indices(5) == (5, 4, 3)   # `div(12, 4) + 1` and `div(4, 2) + 1`
    @test saved_indices(3) isa NTuple{3, Int}
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

@testset "check 8: the growth rate, $T" for T in (Float64, Float32)
    γ = T(0.05)
    t = collect(T, T(0):T(1e-3):T(100))
    W = exp.(2γ .* t) .* cos.(t) .^ 2

    g = growth_rate(t, W)
    @test g isa T
    @test isapprox(g, γ; rtol = 1e-4)

    # a matrix gives one rate per column, of the element type of `W`
    gm = growth_rate(t, hcat(W, W))
    @test gm isa Vector{T}
    @test gm == [g, g]

    # `window = 1` finds the maxima of the sampled `cos²`: without it the flat maximum `Float32`
    # reaches at `h = 1e-3` is near enough to the true one for the fit to stay in its bound
    w1 = growth_rate(t, W; window = 1)
    @test isapprox(w1, γ; rtol = 1e-4)

    # the range of the maxima
    @test occursin("not positive",
        refusal(() -> growth_rate(
            T[0, 1, 2, 3, 4, 5, 6], T[-3, -2, -3, -1, 5, -1, 5]; window = 1, maxima = 1:2)))
    @test occursin("needs two",
        refusal(() -> growth_rate(
            T[0, 1, 2, 3, 4, 5, 6], T[-3, -2, -3, -1, 5, -1, 5]; window = 1, maxima = 1:1)))
    @test occursin("beyond the 2 maxima found",
        refusal(() -> growth_rate(
            T[0, 1, 2, 3, 4, 5, 6], T[-3, -2, -3, -1, 5, -1, 5]; window = 1)))
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
