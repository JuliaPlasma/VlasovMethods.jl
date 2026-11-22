using GeometricBrackets: poisson_apply
using GeometricIntegrators: GeometricIntegrators
using LinearAlgebra: I, mul!, norm, qr
using PoissonSolvers: FFTWBasis, Potential
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V

const VM = VlasovMethods

# The oracles of the step: the full-order RK4 step of `_grid_1d1v_rhs!` (Lie at full rank), the
# full-order step of the frozen-field composition of Einkemmer–Lubich's Algorithm 1 (Strang at full
# rank), the exactness theorem of Lubich–Oseledets on an analytic path of rank r (Lie), and the
# orders of the two compositions by successive differences.

rel(a, b) = norm(a .- b) / norm(b)

# the unit Maxwellian, in the precision of `v`
maxwellian(v) = exp(-v^2 / 2) / sqrt(oftype(v, 2π))

# columns orthonormal in the weighted product h Σ aᵢ bᵢ
orth(A, h) = Matrix(qr(sqrt(h) .* A).Q)[:, 1:size(A, 2)] ./ sqrt(h)

# the grid model on [0, 4π) × [−10, 10], the parameters of its integrator and the fields
function problem(T, nx, nv)
    xd = (zero(T), T(4π))
    vd = (T(-10), T(10))
    dist = GridDistribution(zeros(T, nx, nv), xd, vd)
    model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
    int = GeometricIntegrator(model, (zero(T), one(T)), T(1) / 10, GeometricIntegrators.RK4())
    space = GridFactorSpace(dist)
    w = LowRankFields1D1V(space, model.work)
    return (dist = dist, space = space, p = int.equation.parameters, w = w,
        x = collect(dist.x), v = collect(dist.v), hx = step(dist.x), hv = step(dist.v))
end

const PROBLEMS = Dict{Any, Any}()
cached_problem(T, nx, nv) = get!(() -> problem(T, nx, nv), PROBLEMS, (T, nx, nv))

# The fixture factors: the weighted-orthonormalised 1, cos(x/2), sin(x/2), cos(x), … of x and
# vʲ M(v), j = 0, …, r − 1, of v, with S the weighted projection hx hv Xᵀ f₀ V of f₀.
function fourier(x, r)
    T = eltype(x)
    cols = [ones(T, length(x))]
    m = 1
    while length(cols) < r
        push!(cols, cos.(m .* x ./ 2))
        length(cols) < r && push!(cols, sin.(m .* x ./ 2))
        m += 1
    end
    return reduce(hcat, cols)
end
hermite(v, r) = reduce(hcat, [v .^ j .* maxwellian.(v) for j in 0:(r - 1)])

function state(P, f0, r; nfixed = 0)
    X = orth(fourier(P.x, r), P.hx)
    V = orth(hermite(P.v, r), P.hv)
    S = (P.hx * P.hv) .* (X' * f0 * V)
    return LowRankDistribution(P.space, X, S, V; nfixed)
end

landau_f0(P, α) = (1 .+ eltype(P.x)(α) .* cos.(P.x ./ 2)) .* maxwellian.(P.v)'

dense(Y) = Y.X * Y.S * Y.V'

# F(f) of the model's `_grid_1d1v_rhs!`
function rhs(f, p)
    ḟ = similar(vec(f))
    VM._grid_1d1v_rhs!(ḟ, zero(eltype(f)), copy(vec(f)), p)
    return reshape(ḟ, size(f))
end

# F_ϕ(f) = [ϕ ⊗ 1 + 1 ⊗ v²/2, f] / (hx hv) at a given ϕ, by the bracket call of the base's test
function rhs_at(f, ϕ, p)
    a = p.arakawa
    h = vec(ϕ .+ p.v' .^ 2 ./ 2)
    return reshape(poisson_apply(a, vec(f), h) ./ (a.hx * a.hv), size(f))
end

potential(f, p) = p.G * (p.arakawa.hv .* vec(sum(f; dims = 2)))

function rk4(F, y, τ)
    k1 = F(y)
    k2 = F(y .+ (τ / 2) .* k1)
    k3 = F(y .+ (τ / 2) .* k2)
    k4 = F(y .+ τ .* k3)
    return y .+ (τ / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

# the full-order step of the Strang composition: an RK4 half step at the frozen start potential
# predicts the midpoint potential, and an RK4 full step runs at it
function strang_full(f, τ, p)
    ϕ0 = potential(f, p)
    ϕh = potential(rk4(g -> rhs_at(g, ϕ0, p), f, τ / 2), p)
    return rk4(g -> rhs_at(g, ϕh, p), f, τ)
end

fullorder(::Val{:lie}, f, τ, p) = rk4(g -> rhs(g, p), f, τ)
fullorder(::Val{:strang}, f, τ, p) = strang_full(f, τ, p)

fullrank_bound(::Type{Float64}) = 1e-10
fullrank_bound(::Type{Float32}) = 1e-5

function orthonormality(Y, P)
    max(maximum(abs, P.hx .* (Y.X' * Y.X) - I), maximum(abs, P.hv .* (Y.V' * Y.V) - I))
end

function check_factors(Y, P, r)
    T = eltype(P.x)
    nx, nv = size(P.space)
    @test orthonormality(Y, P) ≤ 8 * r * eps(T)
    @test size(Y.X) == (nx, r)
    @test size(Y.S) == (r, r)
    @test size(Y.V) == (nv, r)
    @test eltype(Y.X) === eltype(Y.S) === eltype(Y.V) === T
    @test Y.nfixed == 0
end

const COMPOSITIONS = (:lie, :strang)

# the exception that `f()` throws, or `nothing`
thrown(f) =
    try
        f()
        nothing
    catch err
        err
    end

@testset "The method object" begin
    lie, strang = ProjectorSplitting(:lie), ProjectorSplitting(:strang)
    @test typeof(lie) != typeof(strang)
    for bad in (:midpoint, "lie")
        err = thrown(() -> ProjectorSplitting(bad))
        @test err isa ArgumentError
        @test err.msg ==
              "the composition of a projector splitting is :lie or :strang; got " *
              repr(bad)
    end
end

@testset "Full rank reproduces the full-order model ($T, τ = $τ)" for T in (Float64, Float32),
    τ in (0.01, -0.01)

    nx, nv = 16, 32
    P = cached_problem(T, nx, nv)
    f0 = landau_f0(P, 0.5)
    for c in COMPOSITIONS
        Y = state(P, f0, nx)
        method = ProjectorSplitting(c)
        f = copy(f0)
        e = zero(T)
        t = zero(T)
        for _ in 1:20
            @test lowrank_step!(Y, method, P.w, t, T(τ)) === Y
            f = fullorder(Val(c), f, T(τ), P.p)
            e = max(e, rel(dense(Y), f))
            t += T(τ)
        end
        @test e ≤ fullrank_bound(T)
        @test eltype(Y.X) === eltype(Y.S) === eltype(Y.V) === T
    end
end

# An analytic field F(Y, t) = Ȧ(t) of the path A(t) = (X₀ + t X₁) S₀ (V₀ + t V₁)ᵀ, on a factor
# space of nx × nv nodes with the weights hx and hv and no grid.
struct WeightedSpace{T} <: VM.AbstractFactorSpace{T}
    nx::Int
    nv::Int
    hx::T
    hv::T
end

Base.size(s::WeightedSpace) = (s.nx, s.nv)
function weighted_qr(A, h)
    Q = orth(A, h)
    return Q, h .* (Q' * A)
end
VM.orthonormalize_x(s::WeightedSpace, A::AbstractMatrix) = weighted_qr(A, s.hx)
VM.orthonormalize_v(s::WeightedSpace, A::AbstractMatrix) = weighted_qr(A, s.hv)

struct PathField{T}
    X0::Matrix{T}
    X1::Matrix{T}
    S0::Matrix{T}
    V0::Matrix{T}
    V1::Matrix{T}
    hx::T
    hv::T
end

path(F::PathField, t) = (F.X0 .+ t .* F.X1) * F.S0 * (F.V0 .+ t .* F.V1)'
function rate(F::PathField, t)
    F.X1 * F.S0 * (F.V0 .+ t .* F.V1)' .+ (F.X0 .+ t .* F.X1) * F.S0 * F.V1'
end

function VM._lowrank_kfield!(K̇, K, V, F::PathField, t)
    K̇ .= F.hv .* (rate(F, t) * V)
end
function VM._lowrank_sfield!(Ṡ, X, S, V, F::PathField, t)
    Ṡ .= (F.hx * F.hv) .* (X' * rate(F, t) * V)
end
function VM._lowrank_lfield!(L̇, X, L, F::PathField, t)
    L̇ .= F.hx .* (rate(F, t)' * X)
end

# The field depends on no potential, so the Strang composition reads none.
VM._lowrank_potential(X, S, V, F::PathField, t) = nothing
VM._lowrank_potential(X, L, F::PathField, t) = nothing
VM._lowrank_lfield!(L̇, X, L, F::PathField, t, ϕ) = VM._lowrank_lfield!(L̇, X, L, F, t)

function path_problem(T, r)
    nx, nv = 20, 30
    hx, hv = T(0.3), T(0.2)
    X0, X1 = T.(cos.((1:nx) .* (1:r)' ./ 7)), T.(sin.((1:nx) .* (1:r)' ./ 5) ./ 3)
    V0, V1 = T.(cos.((1:nv) .* (1:r)' ./ 9)), T.(sin.((1:nv) .* (1:r)' ./ 4) ./ 3)
    S0 = T[2 0.3 0; 0.1 1 0.2; 0 0.4 0.5][1:r, 1:r]
    F = PathField(X0, X1, S0, V0, V1, hx, hv)
    space = WeightedSpace(nx, nv, hx, hv)
    X, Ŝ = weighted_qr(X0, hx)
    V, R = weighted_qr(V0, hv)
    Y = LowRankDistribution(space, X, Ŝ * S0 * R', V)
    return Y, F
end

@testset "Exactness of Lie on a rank-$r path ($T)" for T in (Float64, Float32), r in (3, 1)

    Y, F = path_problem(T, r)
    @test rel(dense(Y), path(F, zero(T))) ≤ 100 * eps(T)
    method = ProjectorSplitting(:lie)
    t = zero(T)
    for _ in 1:5
        lowrank_step!(Y, method, F, t, T(0.1))
        t += T(0.1)
        @test rel(dense(Y), path(F, t)) ≤ 100 * eps(T)
    end
end

# The time intervals of the Strang substeps. With the potential ignored, the Strang step is the
# K–S half step over [t, t + τ/2], the L step over [t, t + τ] with X fixed, which equals two L half
# steps because the field does not depend on L, and the S–K half step over [t + τ/2, t + τ]: the
# Lie step over [t, t + τ/2] and its adjoint over [t + τ/2, t + τ]. Each is exact on the path, as
# the Lie step is (Lubich–Oseledets, Theorem 4.1, with RK4 exact for the affine Ȧ), the adjoint by
# the same theorem for Aᵀ, so a substep on another interval moves the state off A(t).
@testset "The time intervals of the Strang substeps ($T, r = $r)" for T in (Float64, Float32),
    r in (3, 1)

    Y, F = path_problem(T, r)
    method = ProjectorSplitting(:strang)
    t = zero(T)
    for _ in 1:5
        lowrank_step!(Y, method, F, t, T(0.1))
        t += T(0.1)
        @test rel(dense(Y), path(F, t)) ≤ 100 * eps(T)
    end
end

# A field that logs the kind and the time of each call, and forwards it to the path field.
struct TimeLog{F}
    field::F
    calls::Vector{Tuple{Symbol, Any}}
end
TimeLog(field) = TimeLog(field, Tuple{Symbol, Any}[])

function VM._lowrank_kfield!(K̇, K, V, F::TimeLog, t)
    push!(F.calls, (:K, t))
    VM._lowrank_kfield!(K̇, K, V, F.field, t)
end
function VM._lowrank_sfield!(Ṡ, X, S, V, F::TimeLog, t)
    push!(F.calls, (:S, t))
    VM._lowrank_sfield!(Ṡ, X, S, V, F.field, t)
end
function VM._lowrank_lfield!(L̇, X, L, F::TimeLog, t)
    push!(F.calls, (:L, t))
    VM._lowrank_lfield!(L̇, X, L, F.field, t)
end
function VM._lowrank_lfield!(L̇, X, L, F::TimeLog, t, ϕ)
    push!(F.calls, (:Lϕ, t))
    VM._lowrank_lfield!(L̇, X, L, F.field, t, ϕ)
end
function VM._lowrank_potential(X, S, V, F::TimeLog, t)
    push!(F.calls, (:ϕ, t))
    VM._lowrank_potential(X, S, V, F.field, t)
end
function VM._lowrank_potential(X, L, F::TimeLog, t)
    push!(F.calls, (:ϕ, t))
    VM._lowrank_potential(X, L, F.field, t)
end

# the times of the four stages of one RK4 step over [t, t + h]
stages(kind, t, h) = [(kind, t), (kind, t + h / 2), (kind, t + h / 2), (kind, t + h)]

# The step reads the field at the times of the substeps' intervals, in the state's element type,
# for a `t` given in another float type.
@testset "The times at which the step reads the field ($T)" for T in (Float64, Float32)
    t, τ = T(0.3), T(0.1)
    for (c, expected) in (
        (:lie, [stages(:K, t, τ); stages(:S, t, τ); stages(:L, t, τ)]),
        (:strang,
        [stages(:K, t, τ / 2); stages(:S, t, τ / 2); (:ϕ, t); stages(:Lϕ, t, τ / 2);
         (:ϕ, t + τ / 2); stages(:Lϕ, t, τ); stages(:S, t + τ / 2, τ / 2);
         stages(:K, t + τ / 2, τ / 2)]))
        Y, F = path_problem(T, 3)
        logged = TimeLog(F)
        lowrank_step!(Y, ProjectorSplitting(c), logged, 0.3, 0.1)
        @test logged.calls == expected
        @test all(call -> call[2] isa T, logged.calls)
    end
end

# the rank-4 state of the order runs on 16 × 64 nodes of [0, 4π) × [−10, 10]
function order_f0(P, α)
    T = eltype(P.x)
    x, v = P.x, P.v
    M = maxwellian.(v)
    a = T(α)
    return ones(T, length(x)) .* (M .+ (a / 8) .* v .^ 3 .* M)' .+
           a .* cos.(x ./ 2) .* M' .+ (a / 2) .* sin.(x ./ 2) .* (v .* M)' .+
           (a / 4) .* cos.(x) .* (v .^ 2 .* M)'
end

function slope(lτ, le)
    mτ, me = sum(lτ) / length(lτ), sum(le) / length(le)
    return sum((lτ .- mτ) .* (le .- me)) / sum(abs2, lτ .- mτ)
end

@testset "The orders ($T)" for T in (Float64, Float32)
    P = cached_problem(T, 16, 64)
    f0 = order_f0(P, 0.3)
    τs = (0.1, 0.05, 0.025)
    for (c, p) in ((:lie, 1), (:strang, 2))
        method = ProjectorSplitting(c)
        Ys = map(τs) do τ
            Y = state(P, f0, 4)
            t = zero(T)
            for _ in 1:round(Int, 0.5 / τ)
                lowrank_step!(Y, method, P.w, t, T(τ))
                t += T(τ)
            end
            if τ == first(τs)
                check_factors(Y, P, 4)
            end
            dense(Y)
        end
        e = [norm(Ys[i] - Ys[i + 1]) for i in 1:(length(Ys) - 1)]
        s = slope(log.(collect(τs[1:(end - 1)])), log.(Float64.(e)))
        @test abs(s - p) ≤ 0.2
    end
end

@testset "The factors after a step ($T)" for T in (Float64, Float32)
    nx, nv = 16, 32
    P = cached_problem(T, nx, nv)
    f0 = landau_f0(P, 0.5)
    for c in COMPOSITIONS
        Y = state(P, f0, nx; nfixed = 1)
        @test Y.nfixed == 1
        method = ProjectorSplitting(c)
        t = zero(T)
        for _ in 1:20
            lowrank_step!(Y, method, P.w, t, T(0.01))
            t += T(0.01)
        end
        check_factors(Y, P, nx)
    end
end

@testset "Two states, one method and one field ($T)" for T in (Float64, Float32)
    nx, nv = 16, 32
    P = cached_problem(T, nx, nv)
    for c in COMPOSITIONS
        method = ProjectorSplitting(c)
        Y1 = state(P, landau_f0(P, 0.5), nx)
        Y2 = state(P, landau_f0(P, 0.2), nx)
        lowrank_step!(Y1, method, P.w, zero(T), T(0.01))
        X, S, V = copy(Y1.X), copy(Y1.S), copy(Y1.V)
        lowrank_step!(Y2, method, P.w, zero(T), T(0.01))
        @test Y1.X == X
        @test Y1.S == S
        @test Y1.V == V
    end
end

@testset "The decided edges of the step ($T)" for T in (Float64, Float32)
    nx, nv = 16, 32
    P = cached_problem(T, nx, nv)
    f0 = landau_f0(P, 0.5)
    pattern = r"^the step τ of a projector splitting is finite; got (NaN|Inf|-Inf)(32)?$"
    for c in COMPOSITIONS
        method = ProjectorSplitting(c)

        # a non-finite step
        for τ in (T(NaN), T(Inf), T(-Inf), NaN, Inf, -Inf)
            Y = state(P, f0, nx)
            X, S, V = copy(Y.X), copy(Y.S), copy(Y.V)
            err = thrown(() -> lowrank_step!(Y, method, P.w, zero(T), τ))
            @test err isa ArgumentError
            @test err isa ArgumentError && occursin(pattern, err.msg)
            @test Y.X == X && Y.S == S && Y.V == V
        end

        # a zero step
        Y = state(P, f0, nx)
        Y0 = dense(Y)
        lowrank_step!(Y, method, P.w, zero(T), zero(T))
        @test rel(dense(Y), Y0) ≤ 8 * eps(T)

        # a step of another float type
        if T == Float32
            Ya, Yb = state(P, f0, nx), state(P, f0, nx)
            lowrank_step!(Ya, method, P.w, 0, 0.05)
            lowrank_step!(Yb, method, P.w, 0, 0.05f0)
            @test Ya.X == Yb.X && Ya.S == Yb.S && Ya.V == Yb.V
            @test eltype(Ya.X) === eltype(Ya.S) === eltype(Ya.V) === Float32
        end

        # the zero state
        Y = state(P, f0, nx)
        Y.S = zero(Y.S)
        lowrank_step!(Y, method, P.w, zero(T), T(0.05))
        @test iszero(Y.S)
        @test orthonormality(Y, P) ≤ 8 * nx * eps(T)

        # a non-finite state
        Y = state(P, f0, nx)
        Y.S[1, 1] = T(NaN)
        lowrank_step!(Y, method, P.w, zero(T), T(0.01))
        @test !all(isfinite, Y.S)
    end
end
