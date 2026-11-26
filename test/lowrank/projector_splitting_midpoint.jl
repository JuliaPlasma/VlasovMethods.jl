using GeometricBrackets: Arakawa, poisson_apply
using LinearAlgebra: Diagonal, I, dot, norm, qr, svd
using PoissonSolvers: FFTWBasis, Potential
using Random: Xoshiro, randn
using Test
using VlasovMethods
using VlasovMethods: BoundedVelocityClosure, LowRankFields1D1V

const VM = VlasovMethods

# Projector splitting with implicit-midpoint substeps of the projected-gradient fields. Its
# oracles are the conserved H and L² norm, the plain-gradient control that does not keep H, the
# orders and the symmetry of Strang by successive steps, the agreement of Lie and Strang, the
# dense bracket for the factored fields, the midpoint equation for the substep solves, and the
# central difference for their products. The rank-4 runs carry f to the ends of the 8 × 16 grid
# of v ∈ [−6, 6] by their own dynamics, and a wider velocity grid does not prevent it; energy and
# L² are kept by the antisymmetry of the doubly periodic bracket, which holds for any f, so the
# test needs no f that vanishes at ±vmax.

rel(a, b) = norm(a .- b) / norm(b)

# columns orthonormal in the weighted product h Σ aᵢ bᵢ
orth(A, h) = Matrix(qr(sqrt(h) .* A).Q)[:, 1:size(A, 2)] ./ sqrt(h)

# the grid model on [0, 4π) × [−6, 6] and its fields
function problem(T, nx, nv)
    xd = (zero(T), T(4π))
    vd = (T(-6), T(6))
    dist = GridDistribution(zeros(T, nx, nv), xd, vd)
    model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
    space = GridFactorSpace(dist)
    w = LowRankFields1D1V(space, model.work)
    return (space = space, w = w, G = model.work, x = collect(dist.x), v = collect(dist.v),
        hx = step(dist.x), hv = step(dist.v), nx = nx, nv = nv)
end

const PROBLEMS = Dict{Any, Any}()
cached_problem(T, nx = 8, nv = 16) = get!(() -> problem(T, nx, nv), PROBLEMS, (T, nx, nv))

# f = M(v) [1 + 0.3 cos kx + 0.2 v sin kx + 0.1 (v² − 1) cos 2kx + 0.05 (v³ − 3v) sin 2kx],
# k = 1/2, whose weighted SVD on 8 × 16 has the separated singular values 1.7, 0.18, 0.084, 0.048,
# so S of rank 4 is nonsingular; truncated to rank r in Float64 and rounded to T
function initial_state(P, r)
    T = eltype(P.x)
    x, v = Float64.(P.x), Float64.(P.v)
    hx, hv = Float64(P.hx), Float64(P.hv)
    k = 0.5
    M = exp.(-v .^ 2 ./ 2) ./ sqrt(2π)
    f = (ones(length(x)) .* M') .+ 0.3 .* (cos.(k .* x) .* M') .+
        0.2 .* (sin.(k .* x) .* (v .* M)') .+
        0.1 .* (cos.(2k .* x) .* ((v .^ 2 .- 1) .* M)') .+
        0.05 .* (sin.(2k .* x) .* ((v .^ 3 .- 3v) .* M)')
    F = svd(sqrt(hx * hv) .* f)
    X = F.U[:, 1:r] ./ sqrt(hx)
    V = F.V[:, 1:r] ./ sqrt(hv)
    S = Matrix(Diagonal(F.S[1:r])) ./ sqrt(hx * hv)
    return LowRankDistribution(P.space, T.(X), T.(S), T.(V))
end

dense(Y) = Y.X * Y.S * Y.V'

# the Hamiltonian H = hx hv Σ (v²/2) f + (hx/2) ρᵀGρ, ρ = hv Σⱼ f, and the L² norm hx hv Σ f²
density(P, f) = P.hv .* vec(sum(f; dims = 2))
function energy(P, f)
    ρ = density(P, f)
    return P.hx * P.hv * sum((P.v .^ 2 ./ 2)' .* f) + P.hx / 2 * dot(ρ, P.G * ρ)
end
l2(P, f) = P.hx * P.hv * sum(abs2, f)

# the solve tolerance of the tests: 1e-14 in Float64, the default in Float32
testtol(::Type{Float64}) = 1e-14
testtol(::Type{Float32}) = nothing

function midpoint(T, c; kwargs...)
    ProjectorSplitting(c; substeps = ImplicitMidpointSubsteps(; tol = testtol(T), kwargs...))
end

# the drift bound of a run, 8 · (substeps per step) · (steps) · eps(T): 8 eps(T) per substep,
# summed over the substeps of the run
substeps(c) = c === :lie ? 3 : 5
bound(T, c, n) = 8 * substeps(c) * n * eps(T)

# the largest relative drift of H and of the L² norm over n steps of τ
function drifts(P, Y, method, n, τ)
    f = dense(Y)
    H₀, C₀ = energy(P, f), l2(P, f)
    dH = dC = zero(eltype(f))
    for s in 1:n
        lowrank_step!(Y, method, P.w, (s - 1) * τ, τ)
        f = dense(Y)
        dH = max(dH, abs(energy(P, f) - H₀) / abs(H₀))
        dC = max(dC, abs(l2(P, f) - C₀) / C₀)
    end
    return dH, dC
end

function eltypes_are(Y, T)
    eltype(Y.X) === eltype(Y.S) === eltype(Y.V) === T
end

factors(Y) = (copy(Y.X), copy(Y.S), copy(Y.V))
unchanged(Y, F) = Y.X == F[1] && Y.S == F[2] && Y.V == F[3]

const PRECISIONS = (Float64, Float32)

@testset "the substeps keyword of ProjectorSplitting" begin
    sub = ImplicitMidpointSubsteps()
    @test sub.projected === true
    @test sub.tol === nothing
    @test sub.max_iterations == 50
    for c in (:lie, :strang)
        @test ProjectorSplitting(c) isa ProjectorSplitting{c, Nothing}
        @test ProjectorSplitting(c).substeps === nothing
        m = ProjectorSplitting(c; substeps = sub)
        @test m isa ProjectorSplitting{c, typeof(sub)}
        @test m.substeps === sub
    end
    # the default substeps keep the RK4 methods, and the implicit-midpoint substeps take their own
    P = cached_problem(Float64)
    Y = initial_state(P, 4)
    for c in (:lie, :strang)
        args = (Y, ProjectorSplitting(c), P.w, 0.0, 0.05)
        @test which(lowrank_step!, typeof.(args)) ==
              which(lowrank_step!, Tuple{typeof(Y), ProjectorSplitting{c}, typeof(P.w),
            Float64, Float64})
        margs = (Y, ProjectorSplitting(c; substeps = sub), P.w, 0.0, 0.05)
        @test which(lowrank_step!, typeof.(margs)) != which(lowrank_step!, typeof.(args))
    end
    @test_throws r"^ArgumentError: the substeps of a projector splitting are nothing or ImplicitMidpointSubsteps; got 3$" ProjectorSplitting(
        :lie; substeps = 3)
end

@testset "the factored field of a separable gradient" begin
    for T in PRECISIONS
        P = cached_problem(T)
        (; nx, nv, w) = P
        hx, hv = Float64(P.hx), Float64(P.hv)
        # the dense bracket in Float64 on the grid of T
        a64 = Arakawa(nx, nv, hx, hv)
        J(f, h) = reshape(poisson_apply(a64, vec(f), vec(h)), nx, nv) ./ (hx * hv)
        tol = T == Float64 ? 1e-12 : 4096 * eps(Float32)
        for seed in 1:5
            rng = Xoshiro(seed)
            r = 4
            K, X, L = (T.(randn(rng, nx, r)) for _ in 1:3)
            V, Lv = (T.(randn(rng, nv, r)) for _ in 1:2)
            S = T.(randn(rng, r, r))
            a₁, a₂ = T.(randn(rng, nx)), T.(randn(rng, nx))
            b₁, b₂ = T.(randn(rng, nv)), T.(randn(rng, nv))
            h = Float64.(a₁) * Float64.(b₁)' .+ Float64.(a₂) * Float64.(b₂)'
            K64, X64, V64, S64, L64 = Float64.(K), Float64.(X), Float64.(V), Float64.(S),
            Float64.(Lv)

            K̇ = VM._separable_kfield!(similar(K), K, V, w, a₁, b₁, a₂, b₂)
            @test eltype(K̇) === T
            @test rel(K̇, J(K64 * V64', h) * (hv .* V64)) ≤ tol

            Ṡ = VM._separable_sfield!(similar(S), X, S, V, w, a₁, b₁, a₂, b₂)
            @test eltype(Ṡ) === T
            @test rel(Ṡ, -hx * hv .* X64' * J(X64 * S64 * V64', h) * V64) ≤ tol

            L̇ = VM._separable_lfield!(similar(Lv), X, Lv, w, a₁, b₁, a₂, b₂)
            @test eltype(L̇) === T
            @test rel(L̇, Matrix((hx .* X64' * J(X64 * L64', h))')) ≤ tol
        end
        # at a₁ = ϕ, b₁ = 1, a₂ = 1, b₂ = v²/2, the K-field `_kfield!` of the Hamiltonian
        rng = Xoshiro(6)
        K, V = T.(randn(rng, nx, 4)), T.(randn(rng, nv, 4))
        ϕ = T.(randn(rng, nx))
        K̇ = VM._separable_kfield!(similar(K), K, V, w, ϕ, ones(T, nv), ones(T, nx),
            P.v .^ 2 ./ 2)
        @test isapprox(K̇, VM._kfield!(similar(K), K, V, w, ϕ); rtol = 8 * eps(T))
    end
end

# the field F and the product dF(y) w of the K-, S- and L-substep, at fixed X and V
function substep_maps(P, X, V, projected)
    w = P.w
    return (
        K = ((out, K) -> VM._midpoint_kfield!(out, K, V, w, projected),
            (out, K, W) -> VM._midpoint_kproduct!(out, K, W, V, w, projected)),
        S = ((out, S) -> VM._midpoint_sfield!(out, X, S, V, w, projected),
            (out, S, W) -> VM._midpoint_sproduct!(out, X, S, W, V, w, projected)),
        L = ((out, L) -> VM._midpoint_lfield!(out, X, L, w, projected),
            (out, L, W) -> VM._midpoint_lproduct!(out, X, L, W, w, projected)))
end

# The substep fields against the dense bracket at the dense gradient h = Gρ ⊗ 1 + 1 ⊗ v²/2 of the
# state of the substep, with the weighted projections P_V g = hv g V Vᵀ and P_X g = hx X Xᵀ g
@testset "the substep fields equal their dense form" begin
    for T in PRECISIONS, projected in (true, false)

        P = cached_problem(T)
        (; nx, nv) = P
        hx, hv = Float64(P.hx), Float64(P.hv)
        a64 = Arakawa(nx, nv, hx, hv)
        J(f, h) = reshape(poisson_apply(a64, vec(f), vec(h)), nx, nv) ./ (hx * hv)
        G, v = Float64.(P.G), Float64.(P.v)
        gradient(f) = (G * (hv .* vec(sum(f; dims = 2)))) .+ (v .^ 2 ./ 2)'
        tol = T == Float64 ? 1e-12 : 4096 * eps(Float32)
        rng = Xoshiro(10)
        r = 4
        X, V = orth(T.(randn(rng, nx, r)), P.hx), orth(T.(randn(rng, nv, r)), P.hv)
        K, S, L = T.(randn(rng, nx, r)), T.(randn(rng, r, r)), T.(randn(rng, nv, r))
        X64, V64, K64, S64, L64 = Float64.(X), Float64.(V), Float64.(K), Float64.(S),
        Float64.(L)
        PV(g) = projected ? g * (hv .* V64 * V64') : g
        PX(g) = projected ? (hx .* X64 * X64') * g : g

        fK = K64 * V64'
        K̇ = VM._midpoint_kfield!(similar(K), K, V, P.w, projected)
        @test eltype(K̇) === T
        @test rel(K̇, J(fK, PV(gradient(fK))) * (hv .* V64)) ≤ tol

        fS = X64 * S64 * V64'
        Ṡ = VM._midpoint_sfield!(similar(S), X, S, V, P.w, projected)
        @test eltype(Ṡ) === T
        @test rel(Ṡ, -hx * hv .* X64' * J(fS, PX(PV(gradient(fS)))) * V64) ≤ tol

        fL = X64 * L64'
        L̇ = VM._midpoint_lfield!(similar(L), X, L, P.w, projected)
        @test eltype(L̇) === T
        @test rel(L̇, Matrix((hx .* X64' * J(fL, PX(gradient(fL))))')) ≤ tol
    end
end

@testset "the product equals the central difference" begin
    for T in PRECISIONS, projected in (true, false)

        P = cached_problem(T)
        rng = Xoshiro(7)
        r = 4
        X, V = orth(T.(randn(rng, P.nx, r)), P.hx), orth(T.(randn(rng, P.nv, r)), P.hv)
        sizes = (K = (P.nx, r), S = (r, r), L = (P.nv, r))
        for (name, (F!, dF!)) in pairs(substep_maps(P, X, V, projected))
            y, z = T.(randn(rng, sizes[name]...)), T.(randn(rng, sizes[name]...))
            F(u) = F!(similar(u), u)
            cd = (F(y .+ z) .- F(y .- z)) ./ 2
            dFz = dF!(similar(y), y, z)
            @test eltype(dFz) === T
            @test rel(dFz, cd) ≤ 64 * eps(T)
        end
    end
end

# The solve returns y₁ with y₁ = y₀ + τ F((y₀ + y₁)/2) to its residual test
# ‖R‖₂ ≤ max(tol ‖y₀‖₂, tol ‖R(y_pred)‖₂), from the predictor y_pred = y₀ + τ F(y₀); the
# residual is evaluated here again, so the bound takes twice the tolerance. A step of −τ, a time
# that runs backward, satisfies the equation of −τ and not this one.
@testset "the substep solve" begin
    for T in PRECISIONS, projected in (true, false), τ in (0.05, -0.05)
        P = cached_problem(T)
        Y = initial_state(P, 4)
        tol = T(something(testtol(T), 32 * eps(T)))
        τT = T(τ)
        y₀s = (K = Y.X * Y.S, S = Y.S, L = Y.V * Y.S')
        for (name, (F!, dF!)) in pairs(substep_maps(P, Y.X, Y.V, projected))
            y₀ = y₀s[name]
            F(u) = F!(similar(u), u)
            R(y) = y .- y₀ .- τT .* F((y₀ .+ y) ./ 2)
            y₁ = VM._midpoint_solve(F!, dF!, y₀, τT, tol, 50, string(name))
            @test eltype(y₁) === T
            @test y₁ != y₀
            limit = 2 * tol * max(norm(y₀), norm(R(y₀ .+ τT .* F(y₀))))
            @test norm(R(y₁)) ≤ limit
            # the residual and the product that the solve takes: R is quadratic, so its central
            # difference with unit step is its product, up to round-off
            residual!, jvp! = VM._midpoint_residual(F!, dF!, y₀, τT)
            x, v = vec(y₁), T.(randn(Xoshiro(9), length(y₀)))
            @test residual!(similar(x), x, nothing) ≈ vec(R(y₁)) rtol = 8 * eps(T)
            cd = (residual!(similar(x), x .+ v, nothing) .-
                  residual!(similar(x), x .- v, nothing)) ./ 2
            @test rel(jvp!(similar(x), x, v, nothing), cd) ≤ 64 * eps(T)
        end
    end
    # the default tolerance is 32 eps(T)
    for T in PRECISIONS
        P = cached_problem(T)
        _, tol = VM._check_midpoint(initial_state(P, 4), ImplicitMidpointSubsteps(), P.w, 0.05)
        @test tol === 32 * eps(T)
    end
end

@testset "energy and L² over 10⁴ steps" begin
    for T in PRECISIONS, c in (:lie, :strang)

        P = cached_problem(T)
        Y = initial_state(P, 4)
        n = 10_000
        dH, dC = drifts(P, Y, midpoint(T, c), n, 0.05)
        @test dH ≤ bound(T, c, n)
        @test dC ≤ bound(T, c, n)
        @test eltypes_are(Y, T)
        @test size(Y.S) == (4, 4)
    end
end

@testset "rank one" begin
    for T in PRECISIONS
        P = cached_problem(T)
        Y = initial_state(P, 1)
        dH, dC = drifts(P, Y, midpoint(T, :lie), 100, 0.05)
        @test dH ≤ bound(T, :lie, 100)
        @test dC ≤ bound(T, :lie, 100)
        @test eltypes_are(Y, T)
    end
end

@testset "the plain-gradient control drifts" begin
    for T in PRECISIONS
        P = cached_problem(T)
        Y = initial_state(P, 4)
        n = 10_000
        dH, _ = drifts(P, Y, midpoint(T, :lie; projected = false), n, 0.05)
        @test dH ≥ 10 * bound(T, :lie, n)
        @test eltypes_are(Y, T)
    end
end

const ORDER_STEPS = (0.1, 0.05, 0.025, 0.0125)

# f at t = 0.5 for each step τ of ORDER_STEPS
function order_states(T, c)
    P = cached_problem(T)
    return map(ORDER_STEPS) do τ
        Y = initial_state(P, 4)
        for s in 1:round(Int, 0.5 / τ)
            lowrank_step!(Y, midpoint(T, c), P.w, (s - 1) * τ, τ)
        end
        @test eltypes_are(Y, T)
        dense(Y)
    end
end

norm_l2(T, g) = sqrt(Float64(l2(cached_problem(T), g)))

# the successive differences ‖f_τ − f_{τ/2}‖, and their fitted slope against log τ
differences(T, fs) = [norm_l2(T, fs[i] - fs[i + 1]) for i in 1:(length(fs) - 1)]
function fitted_order(T, fs)
    x = log.(collect(ORDER_STEPS[1:(end - 1)]))
    return ([ones(length(x)) x] \ log.(differences(T, fs)))[2]
end

@testset "the orders" begin
    states = Dict((T, c) => order_states(T, c) for T in PRECISIONS, c in (:lie, :strang))
    @test abs(fitted_order(Float64, states[(Float64, :lie)]) - 1) ≤ 0.2
    @test abs(fitted_order(Float32, states[(Float32, :lie)]) - 1) ≤ 0.2
    # Float32 Strang reaches the floor of the solve tolerance at τ ≤ 0.05, so it is not fitted
    @test abs(fitted_order(Float64, states[(Float64, :strang)]) - 2) ≤ 0.2
    # Lie and Strang converge to the same flow: at the smallest τ, the error of the first-order
    # Lie is about its difference δ(2τ) = ‖f_2τ − f_τ‖ (measured ratio 0.99 in Float64, 1.1 in
    # Float32), and the error of Strang is smaller; a Strang step whose substeps do not add up to
    # τ converges to another limit (ratio 48 with the L-substep over τ/2)
    for T in PRECISIONS
        lie, strang = states[(T, :lie)], states[(T, :strang)]
        @test norm_l2(T, lie[end] - strang[end]) ≤ 2 * differences(T, lie)[end]
    end
end

@testset "Strang is symmetric" begin
    for T in PRECISIONS
        P = cached_problem(T)
        lim = T == Float64 ? 1e-12 : 1e-4
        err = map((:strang, :lie)) do c
            Y = initial_state(P, 4)
            f₀ = dense(Y)
            lowrank_step!(Y, midpoint(T, c), P.w, 0.0, 0.1)
            lowrank_step!(Y, midpoint(T, c), P.w, 0.1, -0.1)
            @test eltypes_are(Y, T)
            sqrt(l2(P, dense(Y) - f₀) / l2(P, f₀))
        end
        @test err[1] ≤ lim
        @test err[2] > 10 * lim
    end
end

# the exception of a call, or nothing
caught(f) =
    try
        f()
        nothing
    catch e
        e
    end

# the message of an exception, and "" for no exception, so that a missing throw fails a test
message(e) = e isa Exception && hasfield(typeof(e), :msg) ? e.msg : ""

@testset "a solve that does not converge" begin
    # one Newton iteration ends the first solve of Lie, the K-substep over τ, with MAXITERS; the
    # first K-substep of Strang, over τ/2, can converge in it, and a later substep ends
    first_failing = (lie = "K", strang = "[KSL]")
    for T in PRECISIONS, c in (:lie, :strang)

        P = cached_problem(T)
        Y = initial_state(P, 4)
        F = factors(Y)
        e = caught(() -> lowrank_step!(
            Y, midpoint(T, c; max_iterations = 1), P.w, 0.0, 0.05))
        @test e isa ErrorException
        pattern = Regex("^the $(first_failing[c])-substep solve of ImplicitMidpointSubsteps " *
                        "returned MAXITERS, with ‖R‖₂/‖y₀‖₂ = ")
        @test occursin(pattern, message(e))
        @test unchanged(Y, F)
    end
    # a state far outside the step size of the solve: in Float32 the first solve ends NONFINITE
    P = cached_problem(Float32)
    Y = initial_state(P, 4)
    Y.S = Y.S .* 3.0f8
    F = factors(Y)
    e = caught(() -> lowrank_step!(Y, midpoint(Float32, :lie), P.w, 0.0, 0.05))
    @test e isa ErrorException
    @test occursin(
        r"^the [KSL]-substep solve of ImplicitMidpointSubsteps returned NONFINITE",
        message(e))
    @test unchanged(Y, F)
end

@testset "a small state" begin
    for T in PRECISIONS
        P = cached_problem(T)
        Y = initial_state(P, 4)
        Y.S = Y.S .* T(3e-8)
        dH, dC = drifts(P, Y, midpoint(T, :lie), 100, 0.05)
        @test dH ≤ bound(T, :lie, 100)
        @test dC ≤ bound(T, :lie, 100)
        @test eltypes_are(Y, T)
    end
end

@testset "inputs refused before a solve" begin
    @test_throws r"^ArgumentError: ImplicitMidpointSubsteps takes max_iterations ≥ 1; got 0$" ImplicitMidpointSubsteps(;
        max_iterations = 0)
    for T in PRECISIONS, c in (:lie, :strang)

        P = cached_problem(T)
        Y = initial_state(P, 4)
        F = factors(Y)
        method = midpoint(T, c)
        for τ in (NaN, Inf, -Inf)
            e = caught(() -> lowrank_step!(Y, method, P.w, 0.0, τ))
            @test e isa ArgumentError
            @test startswith(message(e), "the step τ is not finite in the element type $T")
            @test unchanged(Y, F)
        end
        e = caught(() -> lowrank_step!(Y, method, P.w, 0.0, 0.0))
        @test e isa ArgumentError
        @test startswith(message(e), "ImplicitMidpointSubsteps takes a nonzero step τ; got τ = 0.0")
        @test unchanged(Y, F)
        for name in (:X, :S, :V)
            Z = initial_state(P, 4)
            A = copy(getfield(Z, name))
            A[2, 1] = T(NaN)
            setfield!(Z, name, A)
            G = factors(Z)
            e = caught(() -> lowrank_step!(Z, method, P.w, 0.0, 0.05))
            @test e isa ArgumentError
            @test startswith(message(e),
                "ImplicitMidpointSubsteps takes a state with finite entries; the factor $name")
            @test isequal(Z.X, G[1]) && isequal(Z.S, G[2]) && isequal(Z.V, G[3])
        end
        for tol in (1.0, 2.0, Inf, NaN, eps(T) / 2)
            m = ProjectorSplitting(c; substeps = ImplicitMidpointSubsteps(; tol))
            e = caught(() -> lowrank_step!(Y, m, P.w, 0.0, 0.05))
            @test e isa ArgumentError
            @test startswith(message(e), "ImplicitMidpointSubsteps takes a tol in [eps($T), 1)")
            @test unchanged(Y, F)
        end
        # the edge eps(T) itself is taken
        m = ProjectorSplitting(c; substeps = ImplicitMidpointSubsteps(; tol = eps(T)))
        @test caught(() -> lowrank_step!(initial_state(P, 4), m, P.w, 0.0, 0.05)) ===
              nothing
        # a field of the bounded velocity closure is refused
        wb = LowRankFields1D1V(P.space, P.G; closure = BoundedVelocityClosure())
        e = caught(() -> lowrank_step!(Y, method, wb, 0.0, 0.05))
        @test e isa ArgumentError
        @test startswith(message(e),
            "ImplicitMidpointSubsteps takes a field of the PeriodicVelocityClosure")
        @test unchanged(Y, F)
    end
    # tol = 1e-14 on a Float32 state, below eps(Float32)
    P = cached_problem(Float32)
    Y = initial_state(P, 4)
    F = factors(Y)
    m = ProjectorSplitting(:lie; substeps = ImplicitMidpointSubsteps(; tol = 1e-14))
    e = caught(() -> lowrank_step!(Y, m, P.w, 0.0, 0.05))
    @test e isa ArgumentError
    @test startswith(message(e), "ImplicitMidpointSubsteps takes a tol in [eps(Float32), 1)")
    @test unchanged(Y, F)
end

@testset "the zero state" begin
    for T in PRECISIONS, c in (:lie, :strang)

        P = cached_problem(T)
        Y = initial_state(P, 4)
        Y.S = zero(Y.S)
        lowrank_step!(Y, midpoint(T, c), P.w, 0.0, 0.05)
        @test iszero(Y.S)
        @test all(isfinite, Y.X) && all(isfinite, Y.V)
        @test maximum(abs, P.hx .* (Y.X' * Y.X) - I) ≤ 64 * eps(T)
        @test maximum(abs, P.hv .* (Y.V' * Y.V) - I) ≤ 64 * eps(T)
        @test eltypes_are(Y, T)
    end
end

# the least of five readings of one call, after one warm-up call
function least_allocated(f!::F, args::Vararg{Any, N}) where {F, N}
    f!(args...)
    return minimum(@allocated(f!(args...)) for _ in 1:5)
end

@testset "no f in the substep fields" begin
    T = Float64
    P = problem(T, 256, 512)
    rng = Xoshiro(8)
    r = 4
    X, V = orth(randn(rng, P.nx, r), P.hx), orth(randn(rng, P.nv, r), P.hv)
    K, Kw = randn(rng, P.nx, r), randn(rng, P.nx, r)
    S, Sw = randn(rng, r, r), randn(rng, r, r)
    L, Lw = randn(rng, P.nv, r), randn(rng, P.nv, r)
    w = P.w
    limit = 256 * 512 * 8
    # the control: the helper sees an array of f
    @test least_allocated((K, V) -> K * V', K, V) ≥ limit
    @test least_allocated(VM._midpoint_kfield!, similar(K), K, V, w, true) < limit
    @test least_allocated(VM._midpoint_kproduct!, similar(K), K, Kw, V, w, true) < limit
    @test least_allocated(VM._midpoint_sfield!, similar(S), X, S, V, w, true) < limit
    @test least_allocated(VM._midpoint_sproduct!, similar(S), X, S, Sw, V, w, true) < limit
    @test least_allocated(VM._midpoint_lfield!, similar(L), X, L, w, true) < limit
    @test least_allocated(VM._midpoint_lproduct!, similar(L), X, L, Lw, w, true) < limit
end
