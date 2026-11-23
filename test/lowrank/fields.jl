using GeometricBrackets: Arakawa, poisson_apply
using GeometricIntegrators: GeometricIntegrators
using LinearAlgebra: mul!, norm, qr
using PoissonSolvers: FFTWBasis, Potential
using Random: Xoshiro
using Test
using VlasovMethods
using VlasovMethods: LowRankFields1D1V, _kfield!, _lfield!, _potential!, _sfield!

const VM = VlasovMethods

# Every expected value comes from the full-order right-hand side `_grid_1d1v_rhs!` of the model,
# or from GeometricBrackets' `poisson_apply`, applied to the grid function `f` that the test
# forms from the factors, and projected onto the other factor. The tolerance `1e-12` in `Float64`
# is 4504 eps(Float64); the same budget in units of eps(T), rounded down to a power of 2, gives
# `4096 eps(Float32)` in `Float32`.

tolerance(::Type{Float64}) = 1e-12
tolerance(::Type{Float32}) = 4096 * eps(Float32)

rel(a, b) = norm(a .- b) / norm(b)

# the grid on [0, 4π) × [-6, 6], its model, the parameters of the model's integrator (with the
# buffers `ρ`, `ϕ` and `h`) and the workspace of the fields, both of the velocity closure `closure`
function problem(T, nx, nv, closure)
    xd = (zero(T), T(4π))
    vd = (T(-6), T(6))
    dist = GridDistribution(zeros(T, nx, nv), xd, vd)
    model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
    int = GeometricIntegrator(
        model, (zero(T), one(T)), T(1) / 10, GeometricIntegrators.RK4();
        closure)
    w = LowRankFields1D1V(Arakawa(nx, nv, step(dist.x), step(dist.v)), model.work, dist.v;
        closure)
    return (dist = dist, model = model, p = int.equation.parameters, w = w)
end

const PROBLEMS = Dict{Any, Any}()
function cached_problem(T, nx, nv, closure = PeriodicVelocityClosure())
    get!(() -> problem(T, nx, nv, closure), PROBLEMS, (T, nx, nv, closure))
end

# F(f) of the model's buffered `_grid_1d1v_rhs!`, which also leaves `G (hv Σⱼ f)` in `p.ϕ`
function dense(f, p)
    a = p.arakawa
    ḟ = similar(vec(f))
    VM._grid_1d1v_rhs!(ḟ, zero(eltype(f)), copy(vec(f)), p)
    return reshape(ḟ, a.nx, a.nv)
end

# F_ϕ(f) = [ϕ ⊗ 1 + 1 ⊗ v²/2, f] at a given ϕ, by GeometricBrackets' `poisson_apply`
function dense_at(f, ϕ, p)
    a = p.arakawa
    h = vec(ϕ .+ p.kinetic')
    return reshape(poisson_apply(a, vec(f), h) ./ (a.hx * a.hv), a.nx, a.nv)
end

# The bounded velocity closure: the unknowns are the nodes j = 2 … nv-1, and the field is
# F_b(f) = χ ∘ [Ph, χ ∘ f] with Ph = χ ∘ (ϕ ⊗ 1 + 1 ⊗ v²/2). `cut` sets columns 1 and nv to zero
# by assignment.
function cut(a::AbstractMatrix)
    b = copy(a)
    b[:, [1, end]] .= 0
    return b
end

# F_b(f) of the model's bounded vector field
function dense_bounded(f, p)
    a = p.arakawa
    ḟ = similar(vec(f))
    VM._grid_1d1v_bounded_rhs!(ḟ, zero(eltype(f)), copy(vec(f)), p)
    return reshape(ḟ, a.nx, a.nv)
end

# F_b(f) at a given ϕ, by GeometricBrackets' `poisson_apply` on the cut data
function dense_bounded_at(f, ϕ, p)
    a = p.arakawa
    Ph = cut(ϕ .+ p.kinetic')
    return cut(reshape(poisson_apply(a, vec(cut(f)), vec(Ph)) ./ (a.hx * a.hv), a.nx, a.nv))
end

# the potential G (hv Σ_{j=2}^{nv-1} f) of the bounded closure
bounded_potential(f, p) = p.G * (p.arakawa.hv .* vec(sum(cut(f); dims = 2)))

# standard-normal draws in Float64, rounded to T
draw(rng, T, dims...) = T.(randn(rng, dims...))

function factors(T, nx, nv, r, seed)
    rng = Xoshiro(seed)
    X, V = draw(rng, T, nx, r), draw(rng, T, nv, r)
    S, K, L = draw(rng, T, r, r), draw(rng, T, nx, r), draw(rng, T, nv, r)
    return X, S, V, K, L
end

# `f!(out, args...)` returns its output `out`, and the test reads the output buffer itself
function written(f!, out, args...)
    @test f!(out, args...) === out
    return out
end

kfield(K, V, w) = written(_kfield!, similar(K), K, V, w)
lfield(X, L, w) = written(_lfield!, similar(L), X, L, w)
sfield(X, S, V, w) = written(_sfield!, similar(S), X, S, V, w)
kfield(K, V, w, ϕ) = written(_kfield!, similar(K), K, V, w, ϕ)
lfield(X, L, w, ϕ) = written(_lfield!, similar(L), X, L, w, ϕ)
sfield(X, S, V, w, ϕ) = written(_sfield!, similar(S), X, S, V, w, ϕ)
potential(K, V, w) = written(_potential!, zeros(eltype(K), size(K, 1)), K, V, w)
potential(X, S, V, w) = written(_potential!, zeros(eltype(X), size(X, 1)), X, S, V, w)

# columns orthonormal in the weighted product h Σ aᵢ bᵢ
orth(A, h) = Matrix(qr(sqrt(h) .* A).Q)[:, 1:size(A, 2)] ./ sqrt(h)

# the DLR field P(f)F(f) = K̇Vᵀ + X L̇ᵀ − X Ṡ Vᵀ from the factored fields at f = X S Vᵀ
function dlr(X, S, V, w)
    K̇ = kfield(X * S, V, w)
    L̇ = lfield(X, V * S', w)
    Ṡ = sfield(X, S, V, w)
    return K̇ * V' .+ X * L̇' .- X * Ṡ * V'
end

const SIZES = ((3, 3, 2), (16, 24, 4), (32, 64, 6), (17, 40, 1))
const SEEDS = 1:5

@testset "the factored fields equal the dense evaluation, $T" for T in (Float64, Float32)
    tol = tolerance(T)
    for (nx, nv, r) in SIZES, seed in SEEDS

        (; p, w) = cached_problem(T, nx, nv)
        hx, hv = p.arakawa.hx, p.arakawa.hv
        X, S, V, K, L = factors(T, nx, nv, r, seed)

        K̇ = kfield(K, V, w)
        L̇ = lfield(X, L, w)
        Ṡ = sfield(X, S, V, w)
        @test eltype(K̇) === eltype(L̇) === eltype(Ṡ) === T
        @test rel(K̇, dense(K * V', p) * (hv .* V)) ≤ tol
        @test rel(L̇, dense(X * L', p)' * (hx .* X)) ≤ tol
        # at r = 1 the S-field vanishes; its own testset checks it against a scale
        r > 1 && @test rel(Ṡ, (hx * hv) .* (X' * dense(X * S * V', p) * V)) ≤ tol
    end
end

@testset "the scale extremes, $T" for T in (Float64, Float32)
    (; p, w) = cached_problem(T, 16, 24)
    hv = p.arakawa.hv
    for s in (1e-8, 1e8), seed in SEEDS

        rng = Xoshiro(seed)
        K, V = T.(s .* randn(rng, 16, 4)), draw(rng, T, 24, 4)
        K̇ = kfield(K, V, w)
        @test all(isfinite, K̇)
        @test rel(K̇, dense(K * V', p) * (hv .* V)) ≤ tolerance(T)
    end
end

@testset "the S-field vanishes at rank one, $T" for T in (Float64, Float32)
    (; p, w) = cached_problem(T, 17, 40)
    hx, hv = p.arakawa.hx, p.arakawa.hv
    for seed in SEEDS
        X, S, V, K, L = factors(T, 17, 40, 1, seed)
        F = dense(X * S * V', p)
        Ṡ = sfield(X, S, V, w)
        @test eltype(Ṡ) === T
        @test norm(Ṡ) ≤ 4096 * eps(T) * hx * hv * norm(X) * norm(F) * norm(V)
        @test rel(kfield(K, V, w), dense(K * V', p) * (hv .* V)) ≤ tolerance(T)
        @test rel(lfield(X, L, w), dense(X * L', p)' * (hx .* X)) ≤ tolerance(T)
    end
end

@testset "the potential of factored data, $T" for T in (Float64, Float32)
    tol = tolerance(T)
    for (nx, nv, r) in SIZES, seed in SEEDS

        (; p, w) = cached_problem(T, nx, nv)
        X, S, V, K, L = factors(T, nx, nv, r, seed)

        dense(K * V', p)
        ϕ = potential(K, V, w)
        @test eltype(ϕ) === T
        @test rel(ϕ, p.ϕ) ≤ tol

        dense(X * L', p)
        @test rel(potential(X, L, w), p.ϕ) ≤ tol

        dense(X * S * V', p)
        ϕ = potential(X, S, V, w)
        @test eltype(ϕ) === T
        @test rel(ϕ, p.ϕ) ≤ tol
    end
end

@testset "the fields at a given potential, $T" for T in (Float64, Float32)
    tol = tolerance(T)
    for (nx, nv, r) in SIZES, seed in SEEDS

        (; p, w) = cached_problem(T, nx, nv)
        hx, hv = p.arakawa.hx, p.arakawa.hv
        X, S, V, K, L = factors(T, nx, nv, r, seed)
        ϕ = draw(Xoshiro(100 + seed), T, nx)

        K̇ = kfield(K, V, w, ϕ)
        L̇ = lfield(X, L, w, ϕ)
        Ṡ = sfield(X, S, V, w, ϕ)
        @test eltype(K̇) === eltype(L̇) === eltype(Ṡ) === T
        @test rel(K̇, dense_at(K * V', ϕ, p) * (hv .* V)) ≤ tol
        @test rel(L̇, dense_at(X * L', ϕ, p)' * (hx .* X)) ≤ tol
        F = dense_at(X * S * V', ϕ, p)
        Ṡref = (hx * hv) .* (X' * F * V)
        if r == 1
            @test norm(Ṡ .- Ṡref) ≤ tol * hx * hv * norm(X) * norm(F) * norm(V)
        else
            @test rel(Ṡ, Ṡref) ≤ tol
        end
    end
end

@testset "the field without ϕ is the field at its own potential, $T" for T in (Float64, Float32)
    (; w) = cached_problem(T, 16, 24)
    for seed in SEEDS
        X, S, V, K, L = factors(T, 16, 24, 4, seed)
        @test rel(kfield(K, V, w), kfield(K, V, w, potential(K, V, w))) ≤ 8 * eps(T)
        @test rel(lfield(X, L, w), lfield(X, L, w, potential(X, L, w))) ≤ 8 * eps(T)
        @test rel(sfield(X, S, V, w), sfield(X, S, V, w, potential(X, S, V, w))) ≤
              8 * eps(T)
    end
end

@testset "the full-rank identity, $T" for T in (Float64, Float32)
    for (nx, nv) in ((8, 12), (8, 8), (3, 5)), seed in SEEDS

        (; p, w) = cached_problem(T, nx, nv)
        rng = Xoshiro(seed)
        X = orth(draw(rng, T, nx, nx), p.arakawa.hx)
        V = orth(draw(rng, T, nv, nx), p.arakawa.hv)
        S = draw(rng, T, nx, nx)
        @test rel(dlr(X, S, V, w), dense(X * S * V', p)) ≤ tolerance(T)
    end
end

@testset "the energy rate of the DLR field, $T" for T in (Float64, Float32)
    (; p, w) = cached_problem(T, 16, 24)
    hx, hv = p.arakawa.hx, p.arakawa.hv
    ip(a, b) = hx * hv * sum(a .* b)
    for seed in SEEDS, modes in (true, false)

        rng = Xoshiro(seed)
        X = orth(draw(rng, T, 16, 4), hx)
        Vr = draw(rng, T, 24, 4)
        if modes
            Vr[:, 1] .= one(T)
            Vr[:, 2] .= collect(p.v) .^ 2
        end
        V = orth(Vr, hv)
        S = draw(rng, T, 4, 4)
        f = X * S * V'
        h = p.G * (hv .* vec(sum(f; dims = 2))) .+ p.kinetic'
        PF = dlr(X, S, V, w)
        rate = abs(ip(h, PF)) / sqrt(ip(h, h) * ip(PF, PF))
        if modes
            @test rate ≤ tolerance(T)
        else
            @test rate ≥ 1e-3
        end
    end
end

# nv = 4 is the smallest grid with a velocity stencil among the unknowns
const BOUNDED_SIZES = ((3, 4, 2), (16, 24, 4), (32, 64, 6), (17, 40, 1))

@testset "the bounded fields equal the dense bounded field, $T" for T in (Float64, Float32)
    # the end rows of V and L are drawn too, so a field that read them would differ
    tol = tolerance(T)
    for (nx, nv, r) in BOUNDED_SIZES, seed in SEEDS, frozen in (false, true)
        (; p, w) = cached_problem(T, nx, nv, BoundedVelocityClosure())
        hx, hv = p.arakawa.hx, p.arakawa.hv
        X, S, V, K, L = factors(T, nx, nv, r, seed)
        ϕ = draw(Xoshiro(100 + seed), T, nx)
        Fb(f) = frozen ? dense_bounded_at(f, ϕ, p) : dense_bounded(f, p)

        K̇ = frozen ? kfield(K, V, w, ϕ) : kfield(K, V, w)
        L̇ = frozen ? lfield(X, L, w, ϕ) : lfield(X, L, w)
        Ṡ = frozen ? sfield(X, S, V, w, ϕ) : sfield(X, S, V, w)
        @test eltype(K̇) === eltype(L̇) === eltype(Ṡ) === T
        @test rel(K̇, Fb(K * V') * (hv .* V)) ≤ tol
        @test rel(L̇, Fb(X * L')' * (hx .* X)) ≤ tol
        @test all(==(0), L̇[[1, nv], :])
        F = Fb(X * S * V')
        Ṡref = (hx * hv) .* (X' * F * V)
        if r == 1
            @test norm(Ṡ .- Ṡref) ≤ tol * hx * hv * norm(X) * norm(F) * norm(V)
        else
            @test rel(Ṡ, Ṡref) ≤ tol
        end

        if !frozen
            ϕK = potential(K, V, w)
            @test eltype(ϕK) === T
            @test rel(ϕK, bounded_potential(K * V', p)) ≤ tol
            @test rel(potential(X, S, V, w), bounded_potential(X * S * V', p)) ≤ tol
        end
    end
end

@testset "the bounded DLR field: full rank and energy, $T" for T in (Float64, Float32)
    # V has zero end rows and is orthonormal in hv Σ, so that P(f) is the orthogonal projector of
    # the tangent space in the inner product over the unknowns
    function cut_orth(rng, p, nv, r; modes = false)
        Vr = draw(rng, T, nv, r)
        if modes
            Vr[:, 1] .= one(T)
            Vr[:, 2] .= collect(p.v) .^ 2
        end
        Vr[[1, nv], :] .= 0
        return orth(Vr, p.arakawa.hv)
    end

    for (nx, nv) in ((8, 12), (8, 10), (3, 5)), seed in SEEDS

        (; p, w) = cached_problem(T, nx, nv, BoundedVelocityClosure())
        rng = Xoshiro(seed)
        X = orth(draw(rng, T, nx, nx), p.arakawa.hx)
        V = cut_orth(rng, p, nv, nx)
        S = draw(rng, T, nx, nx)
        @test rel(dlr(X, S, V, w), dense_bounded(X * S * V', p)) ≤ tolerance(T)
    end

    (; p, w) = cached_problem(T, 16, 24, BoundedVelocityClosure())
    hx, hv = p.arakawa.hx, p.arakawa.hv
    ip(a, b) = hx * hv * sum((a .* b)[:, 2:(end - 1)])
    for seed in SEEDS, modes in (true, false)

        rng = Xoshiro(seed)
        X = orth(draw(rng, T, 16, 4), hx)
        V = cut_orth(rng, p, 24, 4; modes)
        S = draw(rng, T, 4, 4)
        f = X * S * V'
        Ph = cut(bounded_potential(f, p) .+ p.kinetic')
        PF = dlr(X, S, V, w)
        rate = abs(ip(Ph, PF)) / sqrt(ip(Ph, Ph) * ip(PF, PF))
        if modes
            @test rate ≤ tolerance(T)
        else
            @test rate ≥ 1e-3
        end
    end
end

# The bytes of one call of `f(args...)` after a warm call, the minimum over three calls, through
# a barrier with concrete argument types.
function allocations(f::F, args::A) where {F, A}
    f(args...)
    return minimum(@allocated(f(args...)) for _ in 1:3)
end

# the K-field formed through the full-order right-hand side, as a control that forms f
function dense_kfield!(K̇, K, V, p)
    f = K * V'
    ḟ = similar(f)
    VM._grid_1d1v_rhs!(vec(ḟ), zero(eltype(f)), vec(f), p)
    mul!(K̇, ḟ, V)
    K̇ .*= p.arakawa.hv
    return K̇
end

@testset "allocations of the fields, $T" for T in (Float64, Float32)
    # `--check-bounds=yes` changes the compiled code that the bound measures, so the measurement
    # does not run rather than being weakened; it is skipped where a run forces
    # `--check-bounds=yes`, and under coverage instrumentation, as the allocation tests of the
    # other test files are
    if Base.JLOptions().check_bounds != 1 && Base.JLOptions().code_coverage == 0
        for (nx, nv) in ((64, 128), (512, 1024)), r in (4, 8)

            xd = (zero(T), T(4π))
            dist = GridDistribution(zeros(T, nx, nv), xd, (T(-6), T(6)))
            model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
            a = Arakawa(nx, nv, step(dist.x), step(dist.v))
            w = LowRankFields1D1V(a, model.work, dist.v)
            X, S, V, K, L = factors(T, nx, nv, r, 1)
            ϕ = draw(Xoshiro(2), T, nx)
            # a call allocates at most three vectors of length nx (ϕ, ρ and ϕ̃), one of length nv
            # (k̃), two of length r and five r × r matrices; 2048 bytes cover the array headers.
            # One more nx × r array breaks the bound at (512, 1024), one more nv × r at (64, 128).
            bound = (3nx + nv + 2r + 5r^2) * sizeof(T) + 2048

            @test allocations(_kfield!, (similar(K), K, V, w)) ≤ bound
            @test allocations(_lfield!, (similar(L), X, L, w)) ≤ bound
            @test allocations(_sfield!, (similar(S), X, S, V, w)) ≤ bound
            @test allocations(_kfield!, (similar(K), K, V, w, ϕ)) ≤ bound
            @test allocations(_lfield!, (similar(L), X, L, w, ϕ)) ≤ bound
            @test allocations(_sfield!, (similar(S), X, S, V, w, ϕ)) ≤ bound
            @test allocations(_potential!, (similar(ϕ), K, V, w)) ≤ bound
            @test allocations(_potential!, (similar(ϕ), X, S, V, w)) ≤ bound

            if (nx, nv, r) == (512, 1024, 4)
                @test bound < nx * nv * sizeof(T)
                p = (arakawa = a, v = dist.v, G = model.work)
                @test allocations(dense_kfield!, (similar(K), K, V, p)) > bound
            end
        end
    end
end

@testset "the edges of the fields, $T, $(nameof(typeof(closure)))" for T in (Float64, Float32),
    closure in (PeriodicVelocityClosure(), BoundedVelocityClosure())

    nx, nv, r = 16, 24, 4
    (; w) = cached_problem(T, nx, nv, closure)
    X, S, V, K, L = factors(T, nx, nv, r, 1)
    ϕ = draw(Xoshiro(7), T, nx)
    sentinel = T(12345)
    filled(dims...) = fill(sentinel, dims...)

    # sizes that do not fit: a DimensionMismatch before any write, with the whole message
    function misfit(pattern, f!, out, args...)
        @test_throws pattern f!(out, args...)
        @test all(==(sentinel), out)
    end
    kmsg = r"^DimensionMismatch: _kfield! on a 16 × 24 grid takes K̇ and K of 16 rows and V of 24 rows, with \d+ columns each(, and ϕ of length 16)?; got K̇ of size \(\d+, \d+\), K of size \(\d+, \d+\) and V of size \(\d+, \d+\)(, and ϕ of length \d+)?$"
    lmsg = r"^DimensionMismatch: _lfield! on a 16 × 24 grid takes L̇ and L of 24 rows and X of 16 rows, with \d+ columns each(, and ϕ of length 16)?; got L̇ of size \(\d+, \d+\), X of size \(\d+, \d+\) and L of size \(\d+, \d+\)(, and ϕ of length \d+)?$"
    smsg = r"^DimensionMismatch: _sfield! on a 16 × 24 grid takes X of 16 rows, V of 24 rows, and Ṡ and S of size r × r, with r the columns of X and V(, and ϕ of length 16)?; got Ṡ of size \(\d+, \d+\), X of size \(\d+, \d+\), S of size \(\d+, \d+\) and V of size \(\d+, \d+\)(, and ϕ of length \d+)?$"
    pmsg = r"^DimensionMismatch: _potential! on a 16 × 24 grid takes ϕ of length 16, (K of 16 rows and V of 24 rows, with as many columns each|X of 16 rows, V of 24 rows and S of size r × r, with r the columns of X and V); got ϕ of length \d+, (K of size \(\d+, \d+\) and V of size \(\d+, \d+\)|X of size \(\d+, \d+\), S of size \(\d+, \d+\) and V of size \(\d+, \d+\))$"

    misfit(kmsg, _kfield!, filled(nx, r), draw(Xoshiro(1), T, nx + 1, r), V, w)
    misfit(kmsg, _kfield!, filled(nx, r), draw(Xoshiro(1), T, nx + 1, r), V, w, ϕ)
    misfit(kmsg, _kfield!, filled(nx, r), K, draw(Xoshiro(1), T, nv, r + 1), w)
    misfit(kmsg, _kfield!, filled(nx + 1, r), K, V, w)
    misfit(kmsg, _kfield!, filled(nx, r + 1), K, V, w, ϕ)
    misfit(kmsg, _kfield!, filled(nx, r), K, V, w, zeros(T, nx + 1))
    misfit(lmsg, _lfield!, filled(nv, r), X, draw(Xoshiro(1), T, nv - 1, r), w)
    misfit(lmsg, _lfield!, filled(nv, r), X, L, w, zeros(T, nx + 1))
    misfit(lmsg, _lfield!, filled(nv, r), draw(Xoshiro(1), T, nx - 1, r), L, w)
    misfit(lmsg, _lfield!, filled(nv, r), X, draw(Xoshiro(1), T, nv, r + 1), w, ϕ)
    misfit(lmsg, _lfield!, filled(nv - 1, r), X, L, w)
    misfit(lmsg, _lfield!, filled(nv, r - 1), X, L, w, ϕ)
    misfit(smsg, _sfield!, filled(r, r), X, draw(Xoshiro(1), T, r, r + 1), V, w)
    misfit(smsg, _sfield!, filled(r, r), X, draw(Xoshiro(1), T, r + 1, r), V, w, ϕ)
    misfit(smsg, _sfield!, filled(r, r), draw(Xoshiro(1), T, nx + 1, r), S, V, w)
    misfit(smsg, _sfield!, filled(r, r), X, S, draw(Xoshiro(1), T, nv, r - 1), w)
    misfit(smsg, _sfield!, filled(r, r), X, S, V, w, zeros(T, nx - 1))
    misfit(smsg, _sfield!, filled(r, r + 1), X, S, V, w)
    misfit(smsg, _sfield!, filled(r + 1, r), X, S, V, w, ϕ)
    misfit(pmsg, _potential!, filled(nx - 1), K, V, w)
    misfit(pmsg, _potential!, filled(nx), draw(Xoshiro(1), T, nx + 1, r), V, w)
    misfit(pmsg, _potential!, filled(nx), K, draw(Xoshiro(1), T, nv + 1, r), w)
    misfit(pmsg, _potential!, filled(nx), K, draw(Xoshiro(1), T, nv, r - 1), w)
    misfit(pmsg, _potential!, filled(nx - 1), X, S, V, w)
    misfit(pmsg, _potential!, filled(nx), draw(Xoshiro(1), T, nx - 1, r), S, V, w)
    misfit(pmsg, _potential!, filled(nx), X, draw(Xoshiro(1), T, r, r + 1), V, w)
    misfit(pmsg, _potential!, filled(nx), X, S, draw(Xoshiro(1), T, nv + 1, r), w)
    misfit(pmsg, _potential!, filled(nx), X, S, draw(Xoshiro(1), T, nv, r + 1), w)

    # an output that shares memory with an input: an ArgumentError before any write
    # the whole message of the guard, which names the output and the input that it shares
    # memory with
    function aliased(f!, out, input, args...)
        pattern = Regex("^ArgumentError: $(nameof(f!)): the output $out shares memory " *
                        "with the input $input\$")
        @test_throws pattern f!(args...)
    end
    Kc, Lc, Sc = copy(K), copy(L), copy(S)
    aliased(_kfield!, "K̇", "K", Kc, Kc, V, w)
    aliased(_kfield!, "K̇", "K", Kc, Kc, V, w, ϕ)
    aliased(_lfield!, "L̇", "L", Lc, X, Lc, w)
    aliased(_lfield!, "L̇", "L", Lc, X, Lc, w, ϕ)
    aliased(_sfield!, "Ṡ", "S", Sc, X, Sc, V, w)
    aliased(_sfield!, "Ṡ", "S", Sc, X, Sc, V, w, ϕ)
    @test (Kc, Lc, Sc) == (K, L, S)
    # outputs that share memory with the factors of the other direction, at a square grid
    # where the sizes fit
    w3 = cached_problem(T, 3, 3, closure).w
    A3 = draw(Xoshiro(3), T, 3, 3)
    B3 = draw(Xoshiro(4), T, 3, 3)
    C3 = draw(Xoshiro(5), T, 3, 3)
    aliased(_kfield!, "K̇", "V", A3, B3, A3, w3)
    aliased(_lfield!, "L̇", "X", A3, A3, B3, w3)
    aliased(_sfield!, "Ṡ", "X", A3, A3, B3, C3, w3)
    aliased(_sfield!, "Ṡ", "V", A3, B3, C3, A3, w3)
    # an output that shares memory with the given ϕ, or the potential with a factor, through a
    # `Matrix` reshaped from the `Vector` where the lengths fit
    ϕ1 = draw(Xoshiro(8), T, nx)
    V1 = draw(Xoshiro(9), T, nv, 1)
    aliased(_kfield!, "K̇", "ϕ", reshape(ϕ1, nx, 1), draw(Xoshiro(10), T, nx, 1), V1, w, ϕ1)
    aliased(_sfield!, "Ṡ", "ϕ", reshape(ϕ1, r, r), X, S, V, w, ϕ1)
    aliased(_potential!, "ϕ", "K", ϕ1, reshape(ϕ1, nx, 1), V1, w)
    aliased(_potential!, "ϕ", "X", ϕ1, reshape(ϕ1, nx, 1), ones(T, 1, 1), V1, w)
    aliased(_potential!, "ϕ", "S", ϕ1, draw(Xoshiro(10), T, nx, r), reshape(ϕ1, r, r), V, w)
    ϕ3 = draw(Xoshiro(13), T, 3)
    aliased(_lfield!, "L̇", "ϕ", reshape(ϕ3, 3, 1), draw(Xoshiro(14), T, 3, 1),
        draw(Xoshiro(15), T, 3, 1), w3, ϕ3)
    @test ϕ1 == draw(Xoshiro(8), T, nx)
    @test ϕ3 == draw(Xoshiro(13), T, 3)

    # rank zero: the empty output, and a zero potential
    X0, S0, V0, K0, L0 = factors(T, nx, nv, 0, 1)
    @test size(written(_kfield!, similar(K0), K0, V0, w)) == (nx, 0)
    @test size(written(_lfield!, similar(L0), X0, L0, w)) == (nv, 0)
    @test size(written(_sfield!, similar(S0), X0, S0, V0, w)) == (0, 0)
    @test size(written(_kfield!, similar(K0), K0, V0, w, ϕ)) == (nx, 0)
    @test size(written(_lfield!, similar(L0), X0, L0, w, ϕ)) == (nv, 0)
    @test size(written(_sfield!, similar(S0), X0, S0, V0, w, ϕ)) == (0, 0)
    @test written(_potential!, filled(nx), K0, V0, w) == zeros(T, nx)
    @test written(_potential!, filled(nx), X0, S0, V0, w) == zeros(T, nx)

    # zero factors: an exactly zero field and potential, also at a given ϕ
    Z = zero
    @test all(==(0), written(_kfield!, filled(nx, r), Z(K), Z(V), w))
    @test all(==(0), written(_lfield!, filled(nv, r), Z(X), Z(L), w))
    @test all(==(0), written(_sfield!, filled(r, r), Z(X), Z(S), Z(V), w))
    @test all(==(0), written(_kfield!, filled(nx, r), Z(K), Z(V), w, ϕ))
    @test all(==(0), written(_lfield!, filled(nv, r), Z(X), Z(L), w, ϕ))
    @test all(==(0), written(_sfield!, filled(r, r), Z(X), Z(S), Z(V), w, ϕ))
    @test all(==(0), written(_potential!, filled(nx), Z(K), Z(V), w))
    @test all(==(0), written(_potential!, filled(nx), Z(X), Z(S), Z(V), w))

    # a non-finite entry propagates, without an error
    for bad in (T(NaN), T(Inf))
        Kb, Lb, Sb = copy(K), copy(L), copy(S)
        Kb[3, 2] = Lb[5, 1] = Sb[2, 3] = bad
        @test any(!isfinite, written(_kfield!, similar(K), Kb, V, w))
        @test any(!isfinite, written(_kfield!, similar(K), Kb, V, w, ϕ))
        @test any(!isfinite, written(_lfield!, similar(L), X, Lb, w))
        @test any(!isfinite, written(_lfield!, similar(L), X, Lb, w, ϕ))
        @test any(!isfinite, written(_sfield!, similar(S), X, Sb, V, w))
        @test any(!isfinite, written(_sfield!, similar(S), X, Sb, V, w, ϕ))
        @test any(!isfinite, written(_potential!, similar(ϕ), Kb, V, w))
        @test any(!isfinite, written(_potential!, similar(ϕ), X, Sb, V, w))
    end

    # another element type has no method
    U = T === Float32 ? Float64 : Float32
    @test_throws r"^MethodError: no method matching _kfield!\(" _kfield!(
        U.(K), U.(K), U.(V), w)
end

@testset "the workspace" begin
    for T in (Float64, Float32)
        nx, nv, r = 16, 24, 4
        xd = (zero(T), T(4π))
        dist = GridDistribution(zeros(T, nx, nv), xd, (T(-6), T(6)))
        dist.values .= (1 .+ T(0.3) .* cos.(dist.x ./ 2)) .* exp.(-dist.v' .^ 2 ./ 2) .+
                       T(0.2) .* sin.(dist.x) .* dist.v' .* exp.(-(dist.v' .- 1) .^ 2)
        G = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx))).work
        d = LowRankDistribution(dist, r)

        w1 = LowRankFields1D1V(Arakawa(nx, nv, step(dist.x), step(dist.v)), G, dist.v)
        w2 = LowRankFields1D1V(GridFactorSpace(dist), G)
        w3 = LowRankFields1D1V(d, G)
        @test w1.kinetic == collect(dist.v) .^ 2 ./ 2
        @test eltype(w1.kinetic) === T
        X, S, V = d.X, d.S, d.V
        fields(w) = (kfield(X * S, V, w), lfield(X, V * S', w), sfield(X, S, V, w))
        @test fields(w1) == fields(w2) == fields(w3)
        @test all(norm.(fields(w1)) .> 0)

        a = Arakawa(nx, nv, step(dist.x), step(dist.v))
        gmsg(sz) = Regex("^DimensionMismatch: LowRankFields1D1V on a 16 × 24 grid takes G of " *
                         "size \\(16, 16\\); got \\($(sz[1]), $(sz[2])\\)\$")
        @test_throws gmsg((17, 17)) LowRankFields1D1V(a, zeros(T, nx + 1, nx + 1), dist.v)
        @test_throws gmsg((16, 17)) LowRankFields1D1V(a, zeros(T, nx, nx + 1), dist.v)
        @test_throws r"^DimensionMismatch: LowRankFields1D1V on a 16 × 24 grid takes 24 velocity nodes; got 25$" LowRankFields1D1V(
            a, G, range(T(-6), T(6); length = nv + 1))
    end
end

@testset "the bounded workspace, $T" for T in (Float64, Float32)
    nx, nv, r = 16, 24, 4
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, nx, nv), xd, (T(-6), T(6)))
    dist.values .= (1 .+ T(0.3) .* cos.(dist.x ./ 2)) .* exp.(-dist.v' .^ 2 ./ 2) .+
                   T(0.2) .* sin.(dist.x) .* dist.v' .* exp.(-(dist.v' .- 1) .^ 2)
    G = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx))).work
    d = LowRankDistribution(dist, r)
    a = Arakawa(nx, nv, step(dist.x), step(dist.v))
    space = GridFactorSpace(dist)
    bounded = BoundedVelocityClosure()
    periodic = PeriodicVelocityClosure()

    # factors whose end rows of V are not zero, where the two closures differ
    X, S, V = d.X, d.S, d.V
    @test all(!=(0), V[[1, nv], 1])
    fields(w) = (kfield(X * S, V, w), lfield(X, V * S', w), sfield(X, S, V, w),
        potential(X, S, V, w))

    w1 = LowRankFields1D1V(a, G, dist.v; closure = bounded)
    w2 = LowRankFields1D1V(space, G; closure = bounded)
    w3 = LowRankFields1D1V(d, G; closure = bounded)
    @test w1.closure === w2.closure === w3.closure === bounded
    @test fields(w1) == fields(w2) == fields(w3)

    p1 = LowRankFields1D1V(a, G, dist.v)
    p2 = LowRankFields1D1V(space, G)
    p3 = LowRankFields1D1V(d, G)
    @test p1.closure === p2.closure === p3.closure === periodic
    @test fields(p1) == fields(p2) == fields(p3) ==
          fields(LowRankFields1D1V(a, G, dist.v; closure = periodic))
    @test all(fields(p1) .!= fields(w1))

    # the guard's own message, which names the constructor
    tmsg = r"^TypeError: in LowRankFields1D1V, in keyword closure, expected VlasovMethods\.VelocityClosure, got a value of type Symbol$"
    @test_throws tmsg LowRankFields1D1V(a, G, dist.v; closure = :bounded)
    @test_throws tmsg LowRankFields1D1V(space, G; closure = :bounded)
end

@testset "allocations of the bounded fields, $T" for T in (Float64, Float32)
    # skipped where a run forces `--check-bounds=yes`, and under coverage instrumentation, as
    # the allocation testset of the periodic fields is
    if Base.JLOptions().check_bounds != 1 && Base.JLOptions().code_coverage == 0
        for (nx, nv) in ((64, 128), (512, 1024)), r in (4, 8)

            xd = (zero(T), T(4π))
            dist = GridDistribution(zeros(T, nx, nv), xd, (T(-6), T(6)))
            model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
            a = Arakawa(nx, nv, step(dist.x), step(dist.v))
            w = LowRankFields1D1V(a, model.work, dist.v; closure = BoundedVelocityClosure())
            X, S, V, K, L = factors(T, nx, nv, r, 1)
            ϕ = draw(Xoshiro(2), T, nx)
            # an nx × nv array breaks the bound at (512, 1024)
            bound = 32 * (nx + nv) * r * sizeof(T) + 16384

            @test allocations(_kfield!, (similar(K), K, V, w)) ≤ bound
            @test allocations(_lfield!, (similar(L), X, L, w)) ≤ bound
            @test allocations(_sfield!, (similar(S), X, S, V, w)) ≤ bound
            @test allocations(_kfield!, (similar(K), K, V, w, ϕ)) ≤ bound
            @test allocations(_lfield!, (similar(L), X, L, w, ϕ)) ≤ bound
            @test allocations(_sfield!, (similar(S), X, S, V, w, ϕ)) ≤ bound
            @test allocations(_potential!, (similar(ϕ), K, V, w)) ≤ bound
            @test allocations(_potential!, (similar(ϕ), X, S, V, w)) ≤ bound
            (nx, nv, r) == (512, 1024, 4) && @test bound < nx * nv * sizeof(T)
        end
    end
end
