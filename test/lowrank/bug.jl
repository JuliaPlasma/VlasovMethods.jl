using GeometricIntegrators: GeometricIntegrators
using HDF5: attrs, h5open
using LinearAlgebra: I, diag, isdiag, norm, svdvals
using PoissonSolvers: FFTWBasis, Potential
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V

include("../helpers/dense_bug.jl")
include("../helpers/rate_fits.jl")

# Every fixture lies on x ∈ [0, 4π), v ∈ [-10, 10], with `G` the Green's matrix of the model's FFT
# backend. On 64 velocity nodes f at v = ±10 stays at round-off of its largest value over each run;
# on the 16 and 24 nodes of the edge, rank-guard and driver fixtures it does not
# (KNOWN_ISSUES.md K110). The expected values come from the dense transcription of each paper's
# algorithm in `test/helpers/dense_bug.jl`, from the package's full-order RK4 integrator, or from
# the published orders of the three methods.

rel(a, b) = norm(a .- b) / norm(b)
reconstruction(Y) = Y.X * Y.S * Y.V'
factors(Y) = (copy(Y.X), copy(Y.S), copy(Y.V))
copy_state(Y) = LowRankDistribution(Y.space, factors(Y)...; nfixed = Y.nfixed)

maxwellian(v) = exp(-v^2 / 2) / sqrt(2 * oftype(v, π))

# the profile f₀(x, v), by default (1 + α cos(x/2)) M(v), on nx × nv nodes, its model, its
# factorisation `LowRankDistribution(dist, r)`, the workspace of the fields and the parameters of
# the dense transcription
function fixture(
        T, nx, nv; r, α = 0.3, f₀ = (x, v) -> (1 + T(α) * cos(x / 2)) * maxwellian(v),
        field = true, fixed = nothing)
    xd, vd = (zero(T), T(4π)), (T(-10), T(10))
    dist = GridDistribution(zeros(T, nx, nv), xd, vd)
    dist.values .= [f₀(x, v) for x in dist.x, v in dist.v]
    model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
    G = field ? model.work : zeros(T, nx, nx)
    Y = fixed === nothing ? LowRankDistribution(dist, r) :
        LowRankDistribution(dist, r; fixed = fixed)
    w = LowRankFields1D1V(Y, G)
    d = (p = (arakawa = w.arakawa, v = dist.v, G = G), hx = step(dist.x), hv = step(dist.v))
    return (; dist, model, Y, w, d)
end

# the largest |f| on the two end rows of the velocity grid, relative to the largest |f|
edge_ratio(f) = maximum(abs, f[:, [1, end]]) / maximum(abs, f)

function orthonormality(Y)
    max(norm(step(Y.space.x) .* (Y.X' * Y.X) - I),
        norm(step(Y.space.v) .* (Y.V' * Y.V) - I))
end

@testset "the method objects" begin
    @test BUG() isa BUG
    @test MidpointBUG() isa MidpointBUG
    @test AugmentedBUG() isa AugmentedBUG
    m = AugmentedBUG(1e-6; rmin = 2, rmax = 5)
    @test m isa AugmentedBUG
    @test (m.tol, m.rmin, m.rmax) == (1e-6, 2, 5)
    @test AugmentedBUG(0; rmax = 1).rmin == 1

    tolerance_phrase = r"^the tolerance of AugmentedBUG is a finite real ≥ 0; got "
    for tol in (-1e-3, NaN, Inf, -Inf, 1im)
        err = try
            AugmentedBUG(tol; rmax = 4)
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(tolerance_phrase, err.msg)
    end
    bounds_phrase = r"^AugmentedBUG takes integers 1 ≤ rmin ≤ rmax; got rmin = "
    for (rmin, rmax) in ((0, 4), (5, 4), (-1, 0), (2.5, 4), (1, 4.0))
        err = try
            AugmentedBUG(1e-3; rmin, rmax)
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(bounds_phrase, err.msg)
    end
end

for T in (Float64, Float32)
    @testset "one step against the dense transcription, $T" begin
        # the shifted profile (1 + 0.3 cos(x/2)) M(v - sin(x/2)) on 16 × 64 nodes, r = 2, one step
        # of τ = 0.1: the smallest weighted singular values of [K(τ) X₀] and [L(τ) V₀] are 9.1e-3
        # and 6.8e-3, and those of the midpoint's augmentation 1.8e-5 and 6.2e-8, so no augmented
        # direction comes from round-off (at r = 1 and, for the midpoint, at r = 3 one does, and
        # one ulp of the input then moves the dense transcription past the Float64 bound); f at
        # ±10 is at most 2.1e-15 (Float64) and 2e-7 (Float32) of its largest value after the step
        shifted = (x, v) -> (1 + T(0.3) * cos(x / 2)) * maxwellian(v - sin(x / 2))
        fx = fixture(T, 16, 64; r = 2, f₀ = shifted)
        bound = T === Float64 ? 1e-12 : 1e-4
        τ = T(0.1)
        X, S, V = factors(fx.Y)
        d = fx.d
        cases = (
            (BUG(), false, () -> dense_bug(d, X, S, V, τ)),
            (AugmentedBUG(), true, () -> dense_augmented_bug(d, X, S, V, τ)),
            (AugmentedBUG(1e-6; rmin = 1, rmax = 6), true,
                () -> dense_augmented_bug(d, X, S, V, τ; tol = 1e-6, rmin = 1, rmax = 6)),
            (MidpointBUG(), true, () -> dense_midpoint_bug(d, X, S, V, τ)))
        for (method, truncates, dense) in cases
            Y = copy_state(fx.Y)
            @test lowrank_step!(Y, method, fx.w, zero(T), τ) === Y
            Xd, Sd, Vd = dense()
            @test eltype(Y.X) == eltype(Y.S) == eltype(Y.V) == T
            @test rel(reconstruction(Y), Xd * Sd * Vd') ≤ bound
            @test size(Y.S, 1) == size(Sd, 1)
            @test reconstruction(Y) != X * S * V'
            if truncates
                @test isdiag(Y.S)
                @test issorted(diag(Y.S); rev = true)
                @test all(≥(0), diag(Y.S))
            end
        end
        @test size(dense_bug(d, X, S, V, τ)[2], 1) == 2
        @test size(dense_augmented_bug(d, X, S, V, τ)[2], 1) == 2
        @test size(dense_midpoint_bug(d, X, S, V, τ)[2], 1) == 2
        adaptive = dense_augmented_bug(d, X, S, V, τ; tol = 1e-6, rmin = 1, rmax = 6)
        @test size(adaptive[2], 1) == 4
    end
end

for T in (Float64, Float32)
    @testset "the truncation, $T" begin
        # rank 4 on 16 × 64 nodes, τ = 0.05: the tolerance √(σ₃ σ₄) of the untruncated step lies
        # between two singular values, and the discarded tail lies far above round-off
        fx = fixture(T, 16, 64; r = 4)
        τ = T(0.05)
        Yu = copy_state(fx.Y)
        lowrank_step!(Yu, AugmentedBUG(0; rmin = 1, rmax = 8), fx.w, zero(T), τ)
        σ = diag(Yu.S)
        @test isdiag(Yu.S) && size(Yu.S, 1) == 8

        tol = sqrt(σ[3] * σ[4])
        Yt = copy_state(fx.Y)
        lowrank_step!(Yt, AugmentedBUG(tol; rmin = 1, rmax = 8), fx.w, zero(T), τ)
        r₁ = size(Yt.S, 1)
        discarded = sqrt(fx.d.hx * fx.d.hv) * norm(reconstruction(Yt) - reconstruction(Yu))
        @test r₁ == 3
        @test discarded ≤ tol
        @test abs(discarded - norm(σ[(r₁ + 1):end])) ≤ 64 * eps(T) * norm(σ)
        @test norm(σ[r₁:end]) > tol
    end
end

for T in (Float64, Float32)
    @testset "the rank bounds, $T" begin
        # 20 steps of τ = 0.05 from rank 4 on 16 × 64 nodes, rmin = 2, rmax = 6: tol = 1e-12
        # reaches rmax from the second step on, tol = 0.1 stays at rmin
        fx = fixture(T, 16, 64; r = 4)
        for (tol, expected) in ((1e-12, 6), (0.1, 2))
            Y = copy_state(fx.Y)
            ranks = Int[]
            for _ in 1:20
                lowrank_step!(
                    Y, AugmentedBUG(T(tol); rmin = 2, rmax = 6), fx.w, zero(T), T(0.05))
                push!(ranks, size(Y.S, 1))
            end
            @test all(r -> 2 ≤ r ≤ 6, ranks)
            settled = tol == 0.1 ? ranks : ranks[2:end]
            @test all(==(expected), settled)
        end
    end

    @testset "a rank outside the bounds, $T" begin
        fx = fixture(T, 16, 24; r = 4)
        for (rmin, rmax) in ((5, 8), (1, 3))
            Y = copy_state(fx.Y)
            before = factors(Y)
            err = try
                lowrank_step!(Y, AugmentedBUG(1e-8; rmin, rmax), fx.w, zero(T), T(0.05))
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(r"^the rank 4 of the state is not in rmin:rmax = ", err.msg)
            @test factors(Y) == before
        end
    end
end

for T in (Float64, Float32)
    @testset "full rank, $T" begin
        # r = nx ≤ nv ≤ 2nx for the augmented and the midpoint BUG, r = nx = nv for BUG: every
        # projection of a substep is the identity, and each step is the full-order RK4 step
        bound = T === Float64 ? 1e-10 : 1e-4
        τ = T(0.05)
        for (nx, nv) in ((64, 64), (32, 64))
            fx = fixture(T, nx, nv; r = nx)
            int = VlasovMethods.GeometricIntegrator(fx.model, (zero(T), 20τ), τ,
                GeometricIntegrators.RK4())
            sol = GeometricIntegrators.integrate(int.integrator)
            reference = [reshape(sol.q[s], nx, nv) for s in 1:20]
            T === Float64 && @test maximum(edge_ratio, reference) ≤ 1e-14
            methods = (AugmentedBUG(), AugmentedBUG(0; rmin = nx, rmax = nx), MidpointBUG())
            for method in (nx == nv ? (BUG(), methods...) : methods)
                Y = copy_state(fx.Y)
                err = 0.0
                for s in 1:20
                    lowrank_step!(Y, method, fx.w, (s - 1) * τ, τ)
                    err = max(err, rel(reconstruction(Y), reference[s]))
                end
                @test eltype(Y.S) == eltype(Y.X) == eltype(Y.V) == T
                @test size(Y.S) == (nx, nx)
                @test err ≤ bound
            end
        end
    end
end

@testset "orders on free streaming" begin
    # G = 0: f(0) = (1 + 0.3 cos(x/2)) M(v) keeps a rank of at most 3, so the error is the time
    # error. The reference is the full-order field with G = 0, stepped by RK4 at τ = 1e-3.
    T = Float64
    fx = fixture(T, 32, 64; r = 3, field = false)
    reference = fx.dist.values
    for _ in 1:1000
        reference = dense_rk4(f -> dense_rhs(fx.d, f), reference, 1e-3)
    end
    τs = (0.05, 0.025, 0.0125, 0.00625)
    function errors(method)
        map(τs) do τ
            Y = copy_state(fx.Y)
            for s in 1:round(Int, 1 / τ)
                lowrank_step!(Y, method, fx.w, (s - 1) * τ, τ)
            end
            rel(reconstruction(Y), reference)
        end
    end
    eb, ea, em = errors(BUG()), errors(AugmentedBUG()), errors(MidpointBUG())
    order(e) = slope(log.(collect(τs)), log.(collect(e)))
    @test abs(order(eb) - 1) ≤ 0.2
    @test order(ea) ≥ 0.8
    @test order(em) ≥ 1.8
    @test all(ea ./ em .≥ 5)
end

for T in (Float64, Float32)
    @testset "the edges of the step, $T" begin
        fx = fixture(T, 16, 24; r = 3)
        methods = (BUG(), AugmentedBUG(), MidpointBUG())

        # a non-finite step is refused before any write
        phrase = r"^the step τ is not finite in the element type (Float64|Float32) of the state; got τ = (NaN|Inf|-Inf)(32)?$"
        for method in methods, τ in (NaN, Inf, -Inf, T(NaN), T(Inf), T(-Inf))

            Y = copy_state(fx.Y)
            before = factors(Y)
            err = try
                lowrank_step!(Y, method, fx.w, zero(T), τ)
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(phrase, err.msg)
            @test factors(Y) == before
        end
        if T === Float32
            # a finite Float64 step that overflows in Float32 is refused, and the message names
            # the step as it was given
            for method in methods, (τ, given) in ((1e300, "1.0e300"), (1.0e39, "1.0e39"))

                Y = copy_state(fx.Y)
                before = factors(Y)
                err = try
                    lowrank_step!(Y, method, fx.w, zero(T), τ)
                catch e
                    e
                end
                @test err isa ArgumentError
                @test err.msg ==
                      "the step τ is not finite in the element type Float32 of " *
                      "the state; got τ = $given"
                @test factors(Y) == before
            end
        end

        # τ of another float type is converted to T before its first use
        for method in methods
            Y₁, Y₂ = copy_state(fx.Y), copy_state(fx.Y)
            lowrank_step!(Y₁, method, fx.w, 0, 0.05)
            lowrank_step!(Y₂, method, fx.w, 0, T(0.05))
            @test factors(Y₁) == factors(Y₂)
            @test eltype(Y₁.S) == eltype(Y₁.X) == eltype(Y₁.V) == T
        end

        # a zero step keeps X S Vᵀ
        for method in methods
            Y = copy_state(fx.Y)
            lowrank_step!(Y, method, fx.w, zero(T), zero(T))
            @test rel(reconstruction(Y), reconstruction(fx.Y)) ≤ 8 * eps(T)
        end

        # the truncation of AugmentedBUG(tol) applies at τ = 0: on the core diag(1, 1/2, 1/4),
        # tol = 0.6 admits the tail (1/2, 1/4), so the step keeps rank 1 and changes X S Vᵀ by
        # the norm √(1/4 + 1/16) of the tail in the weighted Frobenius norm
        Y = LowRankDistribution(fx.Y.space, copy(fx.Y.X), T[1 0 0; 0 0.5 0; 0 0 0.25],
            copy(fx.Y.V))
        Y₀ = reconstruction(Y)
        lowrank_step!(Y, AugmentedBUG(0.6; rmin = 1, rmax = 3), fx.w, zero(T), zero(T))
        @test size(Y.S) == (1, 1)
        change = sqrt(step(fx.dist.x) * step(fx.dist.v)) * norm(reconstruction(Y) - Y₀)
        @test change ≈ sqrt(T(5) / 16) rtol=50 * eps(T)

        # the zero state stays at zero, and a rank-deficient core needs no inverse of S: the
        # profile has rank 2, its factorisation at r = 6 a core of rank 2
        rank2(x, v) = (1 + T(0.3) * cos(x / 2)) * maxwellian(v) +
                      T(0.1) * sin(x / 2) * v * maxwellian(v)
        deficient = fixture(T, 16, 24; r = 6, f₀ = rank2)
        σ = svdvals(deficient.Y.S)
        @test count(>(sqrt(eps(T)) * σ[1]), σ) == 2
        for method in methods
            Y = copy_state(fx.Y)
            Y.S .= 0
            lowrank_step!(Y, method, fx.w, zero(T), T(0.05))
            @test Y.S == zeros(T, 3, 3)
            @test orthonormality(Y) ≤ 50 * eps(T)

            Y = copy_state(deficient.Y)
            lowrank_step!(Y, method, deficient.w, zero(T), T(0.05))
            @test all(isfinite, Y.S)
            @test orthonormality(Y) ≤ 50 * eps(T)
        end

        # a non-finite state is not checked: BUG() keeps its rank-r factors, non-finite as they
        # come, and the augmented methods install the first r columns of the augmented bases and
        # an r × r core of NaN; the status all(isfinite, Y.S) is false, and the factors make a
        # state. Rank 3 on 16 × 24 nodes, and rank 16 on 32 × 16 and 16 × 32 nodes, where the
        # augmented bases have more columns in one direction than in the other.
        augmented = (AugmentedBUG(), AugmentedBUG(1e-6; rmax = 16), MidpointBUG())
        for start in (fx, fixture(T, 32, 16; r = 16), fixture(T, 16, 32; r = 16))
            r, (nx, nv) = size(start.Y.S, 1), size(start.Y)
            for method in (BUG(), augmented...), value in (T(NaN), T(Inf))

                Y = copy_state(start.Y)
                Y.S[1, 1] = value
                @test lowrank_step!(Y, method, start.w, zero(T), T(0.05)) === Y
                @test size(Y.S) == (r, r)
                @test size(Y.X) == (nx, r) && size(Y.V) == (nv, r)
                @test eltype(Y.X) == eltype(Y.S) == eltype(Y.V) == T
                @test !all(isfinite, Y.S)
                method isa BUG || @test all(isnan, Y.S)
                @test LowRankDistribution(Y.space, Y.X, Y.S, Y.V) isa LowRankDistribution{T}
            end
        end

        # the truncation of a non-finite core on finite bases, which the steps above do not
        # reach (their bases are NaN as well): the first r columns of the bases, a core of NaN
        X̂, V̂ = T.(cos.((1:16) .* (1:6)' ./ 7)), T.(sin.((1:24) .* (1:6)' ./ 5))
        Ŝ = Matrix{T}(I, 6, 6)
        Ŝ[2, 3] = T(Inf)
        for method in augmented
            X₁, S₁, V₁ = VlasovMethods._truncate(method, X̂, Ŝ, V̂, 3)
            @test X₁ == X̂[:, 1:3] && V₁ == V̂[:, 1:3]
            @test size(S₁) == (3, 3) && eltype(S₁) == T && all(isnan, S₁)
        end

        # nx > nv: at r = nv the augmented x basis has more columns than the v basis, and one
        # step of the augmented and the midpoint BUG is the full-order RK4 step
        wide = fixture(T, 32, 16; r = 16)
        full_order = dense_rk4(f -> dense_rhs(wide.d, f), wide.dist.values, T(0.05))
        for method in (AugmentedBUG(), MidpointBUG())
            Y = copy_state(wide.Y)
            lowrank_step!(Y, method, wide.w, zero(T), T(0.05))
            @test rel(reconstruction(Y), full_order) ≤ (T === Float64 ? 1e-10 : 1e-4)
        end

        # a step counts the fixed modes in the rank and sets nfixed = 0
        withfixed = fixture(T, 16, 24; r = 6,
            fixed = hcat(ones(T, 24), collect(fx.dist.v),
                collect(fx.dist.v) .^ 2))
        @test withfixed.Y.nfixed == 3
        for method in methods
            Y = copy_state(withfixed.Y)
            lowrank_step!(Y, method, withfixed.w, zero(T), T(0.05))
            @test size(Y.S) == (6, 6)
            @test Y.nfixed == 0
            @test length(Y.X) + length(Y.S) + length(Y.V) == 6 * (16 + 24) + 36
        end

        # two states share one method and one workspace: no factor is bound to a buffer
        other = fixture(T, 16, 24; r = 3, α = 0.2)
        for method in methods
            Y₁, Y₂ = copy_state(fx.Y), copy_state(other.Y)
            lowrank_step!(Y₁, method, fx.w, zero(T), T(0.05))
            kept = factors(Y₁)
            lowrank_step!(Y₂, method, fx.w, zero(T), T(0.05))
            @test factors(Y₁) == kept
        end
    end
end

# a method of the test file whose step doubles S, and writes a NaN at its third call
struct DoublingMethod
    calls::Base.RefValue{Int}
end

function VlasovMethods.lowrank_step!(Y, method::DoublingMethod, field, t, τ)
    method.calls[] += 1
    Y.S .*= 2
    method.calls[] == 3 && (Y.S[1, 1] = NaN)
    return Y
end

# a method of the test file whose step grows the rank by one
struct GrowingMethod end

function VlasovMethods.lowrank_step!(Y, ::GrowingMethod, field, t, τ)
    T, r = eltype(Y.S), size(Y.S, 1)
    Y.X = hcat(Y.X, zeros(T, size(Y.X, 1)))
    Y.V = hcat(Y.V, zeros(T, size(Y.V, 1)))
    Y.S = [Y.S zeros(T, r); zeros(T, 1, r) one(T)]
    return Y
end

# a method of the test file whose step doubles S, and throws `exception` at its second call
struct ThrowingMethod{E <: Exception}
    calls::Base.RefValue{Int}
    exception::E
end

function VlasovMethods.lowrank_step!(Y, method::ThrowingMethod, field, t, τ)
    method.calls[] += 1
    method.calls[] == 2 && throw(method.exception)
    Y.S .*= 2
    return Y
end

# an exception type of the test file, whose `typeof` prints with the name of its module
struct ProbeError <: Exception end

for T in (Float64, Float32)
    @testset "the run driver, $T" begin
        dir = mktempdir()
        fx = fixture(T, 16, 24; r = 3)
        nx, nv = size(fx.Y)

        # the datasets of a run, against the same steps by hand
        method = AugmentedBUG()
        Y = copy_state(fx.Y)
        int = LowRankIntegrator(Y, method, fx.w, (0, 0.25), 0.05)
        file = joinpath(dir, "fixed.h5")
        @test run!(int, file) === Y
        Yh = copy_state(fx.Y)
        times = 0:0.05:0.25
        hand = [factors(Yh)]
        for s in 1:5
            lowrank_step!(Yh, method, fx.w, times[s], 0.05)
            push!(hand, factors(Yh))
        end
        @test factors(Y) == hand[end]
        h5open(file, "r") do h5
            @test read(h5["t"]) == times
            @test read(h5["rank"]) == fill(3, 6)
            X, S, V = read(h5["X"]), read(h5["S"]), read(h5["V"])
            @test eltype(X) == eltype(S) == eltype(V) == T
            @test size(X) == (nx, 3, 6) && size(S) == (3, 3, 6) && size(V) == (nv, 3, 6)
            @test all(n -> (X[:, :, n], S[:, :, n], V[:, :, n]) == hand[n], 1:6)
            @test attrs(h5)["status"] == "ok"
        end

        # the adaptive method: R = min(rmax, nx, nv) columns, zero beyond each rank
        file = joinpath(dir, "adaptive.h5")
        Y = copy_state(fx.Y)
        run!(
            LowRankIntegrator(
                Y, AugmentedBUG(1e-6; rmin = 2, rmax = 8), fx.w, (0, 0.25), 0.05),
            file)
        h5open(file, "r") do h5
            X, S, V, rank = read(h5["X"]), read(h5["S"]), read(h5["V"]), read(h5["rank"])
            @test eltype(X) == eltype(S) == eltype(V) == T
            @test size(X) == (nx, 8, 6) && size(S) == (8, 8, 6) && size(V) == (nv, 8, 6)
            @test rank[1] == 3 && rank[end] == size(Y.S, 1)
            @test all(2 .≤ rank .≤ 8)
            @test S[1:rank[end], 1:rank[end], end] == Y.S
            for n in 1:6
                r = rank[n]
                @test iszero(X[:, (r + 1):end, n]) && iszero(V[:, (r + 1):end, n])
                @test iszero(S[(r + 1):end, :, n]) && iszero(S[:, (r + 1):end, n])
            end
            @test attrs(h5)["status"] == "ok"
        end

        # an rmax above min(nx, nv): R = min(nx, nv)
        file = joinpath(dir, "wide.h5")
        Y = copy_state(fx.Y)
        run!(
            LowRankIntegrator(Y, AugmentedBUG(1e-6; rmin = 2, rmax = 100), fx.w, (0, 0.1),
                0.05), file)
        h5open(file, "r") do h5
            R = min(nx, nv)
            @test size(read(h5["X"])) == (nx, R, 3) && size(read(h5["S"])) == (R, R, 3)
            @test size(read(h5["V"])) == (nv, R, 3)
            @test attrs(h5)["status"] == "ok"
        end

        # a state whose rank is outside rmin:rmax is refused before the file is written
        file = joinpath(dir, "refused.h5")
        Y = copy_state(fx.Y)
        before = factors(Y)
        @test_throws r"^ArgumentError: the rank 3 of the state is not in rmin:rmax = 1:2" run!(
            LowRankIntegrator(
                Y, AugmentedBUG(1e-6; rmin = 1, rmax = 2), fx.w, (0, 0.1), 0.05),
            file)
        @test factors(Y) == before
        @test !isfile(file)

        # a step that is zero or not finite in T is refused before the file is written: 1e-50
        # and 1e300 are finite and positive in Float64, and 0 and Inf in Float32, for every
        # method
        if T === Float32
            for method in (BUG(), ProjectorSplitting(:lie)),
                (tspan, tstep) in (((0, 5e-50), 1e-50), ((0, 1e301), 1e300))

                file = joinpath(dir, "step.h5")
                Y = copy_state(fx.Y)
                before = factors(Y)
                @test_throws r"^ArgumentError: the step [^ ]+ is zero or not finite in the element type Float32 of the state$" run!(
                    LowRankIntegrator(Y, method, fx.w, tspan, tstep), file)
                @test factors(Y) == before
                @test !isfile(file)
            end
        end

        # an element type of the times or of the state that HDF5 does not store is refused
        # before the file is opened, so an existing file keeps its bytes
        file = joinpath(dir, "kept.h5")
        h5open(h5 -> (h5["kept"] = collect(1:4)), file, "w")
        bytes = read(file)
        xd, vd = (Float16(0), Float16(4π)), (Float16(-10), Float16(10))
        half = GridFactorSpace(GridDistribution(zeros(Float16, nx, nv), xd, vd))
        Yhalf = LowRankDistribution(half, Float16.(fx.Y.X), Float16.(fx.Y.S), Float16.(fx.Y.V))
        for (Y, method, tspan, tstep, name) in (
            (copy_state(fx.Y), BUG(), (Int128(0), Int128(2)), Int128(1), "Int128"),
            (Yhalf, DoublingMethod(Ref(0)), (0, 0.2), 0.1, "Float16"))
            before = factors(Y)
            err = try
                run!(LowRankIntegrator(Y, method, fx.w, tspan, tstep), file)
            catch e
                e
            end
            @test err isa ArgumentError
            @test err.msg ==
                  "run! writes the times and the state to HDF5, which has no datatype for $name"
            @test read(file) == bytes
            @test factors(Y) == before
        end

        # a method of the test file that keeps the rank, and a run that fails at its third step
        file = joinpath(dir, "failing.h5")
        Y = copy_state(fx.Y)
        S₀ = copy(Y.S)
        @test run!(LowRankIntegrator(Y, DoublingMethod(Ref(0)), fx.w, (0, 0.5), 0.1), file) ===
              Y
        h5open(file, "r") do h5
            @test attrs(h5)["status"] == "failed at step 3"
            S, rank = read(h5["S"]), read(h5["rank"])
            @test S[:, :, 1] == S₀
            @test S[:, :, 2] == 2S₀ && S[:, :, 3] == 4S₀
            @test rank == [3, 3, 3, 0, 0, 0]
            @test iszero(S[:, :, 4:6])
        end

        # a step that throws sets the status with the name of the exception type, and the
        # exception leaves run!; the states before it are written
        for (exception, name) in ((ErrorException("probe"), "ErrorException"),
            (ProbeError(), "ProbeError"))
            file = joinpath(dir, "throwing.h5")
            Y = copy_state(fx.Y)
            S₀ = copy(Y.S)
            @test_throws exception run!(
                LowRankIntegrator(
                    Y, ThrowingMethod(Ref(0), exception), fx.w, (0, 0.5), 0.1),
                file)
            h5open(file, "r") do h5
                @test attrs(h5)["status"] == "failed at step 2: $name"
                S, rank = read(h5["S"]), read(h5["rank"])
                @test eltype(S) == T
                @test S[:, :, 1] == S₀ && S[:, :, 2] == 2S₀
                @test rank == [3, 3, 0, 0, 0, 0]
                @test iszero(S[:, :, 3:6])
            end
        end

        # a step that grows the rank beyond the stored columns throws, and sets the status
        file = joinpath(dir, "growing.h5")
        Y = copy_state(fx.Y)
        S₀ = copy(Y.S)
        @test_throws r"^ArgumentError: step 1 of the run changes the rank to 4, beyond the 3 columns that run! stores" run!(
            LowRankIntegrator(Y, GrowingMethod(), fx.w, (0, 0.2), 0.1), file)
        h5open(file, "r") do h5
            @test attrs(h5)["status"] ==
                  "failed at step 1: rank 4 beyond the 3 stored columns"
            @test read(h5["rank"]) == [3, 0, 0]
            @test read(h5["S"])[:, :, 1] == S₀
            @test iszero(read(h5["S"])[:, :, 2:3])
        end

        # a span of one time holds the initial state alone
        file = joinpath(dir, "single.h5")
        Y = copy_state(fx.Y)
        before = factors(Y)
        run!(LowRankIntegrator(Y, BUG(), fx.w, (0.5, 0.5), 0.1), file)
        @test factors(Y) == before
        h5open(file, "r") do h5
            @test read(h5["t"]) == [0.5]
            @test read(h5["S"])[:, :, 1] == before[2]
            @test attrs(h5)["status"] == "ok"
        end

        # the guards of the constructor; the last span gives the grid t₀:tstep:t₁ that repeats
        # its first time, 16777216f0 twice for Float32, which `allunique` of the range misses
        tbig = T(2)^precision(T)
        @test collect(tbig:one(T):(tbig + 2)) == [tbig, tbig, tbig + 2]
        for (tspan, tstep) in (((0, 1), 0), ((0, 1), -0.1), ((0, 1), NaN), ((0, 1), Inf),
            ((1, 0), 0.1), ((0, Inf), 0.1), ((-Inf, 0), 0.1), ((tbig, tbig + 2), one(T)))
            err = try
                LowRankIntegrator(copy_state(fx.Y), BUG(), fx.w, tspan, tstep)
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(r"^LowRankIntegrator takes ", err.msg)
        end
    end
end
