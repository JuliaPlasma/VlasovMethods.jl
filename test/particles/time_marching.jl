using HDF5
using LinearAlgebra
using ParticleMethods: ParticleList
using Random
using ReducedBasisMethods: get_DEIM_indices
using ReducedComplexityModeling: ParameterSpace, Parameter
using Test
using VlasovMethods
using VlasovMethods: PoissonField, coefficients, efield!, energy, full_leapfrog_step!,
                     h5save, leapfrog_step!, save_full_solution, save_solution, update!

include("../helpers/allocation_barriers.jl")
include("../helpers/poisson_field.jl")

# A seeded random orthonormal `n × k` basis, from the QR of a Gaussian matrix.
function orthogonal(::Type{T}, n, k; seed = 42) where {T}
    Random.seed!(seed)
    return Matrix(qr(randn(T, n, k)).Q)
end

"The one-sample `ParameterSpace` of the march, which writes into its sample 1."
space(::Type{T}) where {T} = ParameterSpace(Parameter(:β, T(6), T(6), 1))

"""
The parts of one march, so that one of them can be made wrong: the particle list `plt`, the basis
`Ψ`, the reduced field `efield`, the snapshots `SS`, the parameters `ip`, the cache `ic` and the
caller's `params`, which holds `χ`. `ns` defaults to `nₜ + 1`, which saves every step, and `nssnap`
defaults to `ns`, so that the saved times of the snapshots and the count the parameters ask for can
be made to disagree.
"""
function march_fixture(
        ::Type{T}; n = 40, k = 4, nt = 1000, ns = nt + 1, nd = 1, dt = T(1e-4),
        χ = 1, seed = 5, nssnap = ns) where {T}
    field, x, w = poisson_field(T; n = n)
    Ψ = orthogonal(T, n, k; seed = seed)
    X = Ψ' * x
    efield = ReducedElectricField(field, Ψ, X)
    plt = ParticleList(reshape(x, 1, n), reshape(T(0.5) .* x, 1, n), reshape(w, 1, n))
    # `nh` is the degree count of the field's potential, which `save_solution` copies the solved
    # coefficients into
    SS = Snapshots(T, nd, n, length(coefficients(efield)), nssnap, space(T))
    ip = IntegratorParameters(dt, nt, ns, length(coefficients(efield)), n, 1)
    ic = ReducedIntegratorCache(ip, k)
    return (; plt, Ψ, efield, SS, ip, ic, params = (χ = χ,))
end

"The snapshots of one march of the fixture, with the timestep `dt` and the scale `χ`."
function march(::Type{T}, plt, Ψ, efield, k, nt; dt = T(1e-4), χ = 1) where {T}
    n = size(plt.x, 2)
    ns, nh = nt + 1, length(coefficients(efield))
    SS = Snapshots(T, 1, n, nh, ns, space(T))
    ip = IntegratorParameters(dt, nt, ns, nh, n, 1)
    ic = ReducedIntegratorCache(ip, k)
    reduced_integrate_vp(plt, Ψ, (χ = χ,), efield, SS, ip, ic, 1)
    return SS
end

@testset "check 2: the χ scale is a timestep scale, $T" for T in (Float64, Float32)
    # The fixture is large enough that a one-ulp difference of `Δt` in `Float32` — `T(dt·χ)` against
    # `dt·T(χ)` at `χ = 0.3` — grows into saved states that differ; at `n = 40` and `k = 4` it does
    # not.
    n, k, nt = 1000, 20, 1000
    field, x, w = poisson_field(T; n = n)
    Ψ = orthogonal(T, n, k; seed = 3)
    X = Ψ' * x
    plt = ParticleList(reshape(x, 1, n), reshape(T(0.5) .* x, 1, n), reshape(w, 1, n))

    # the two wrappers of `src/particles/electric_field.jl`, whose `update!`s ignore `t`
    for efield in (ReducedElectricField(field, Ψ, X),
        DEIMElectricField(field, Ψ, Ψ, get_DEIM_indices(Ψ), X))

        # `(dt, χ = 2)` is `(2dt, χ = 1)` state for state: every product is exact at `χ = 2`, so
        # this is the assertion that `χ` reaches the timestep at all
        @test march(T, plt, Ψ, efield, k, nt; dt = T(1e-4), χ = 2) ==
              march(T, plt, Ψ, efield, k, nt; dt = T(2e-4), χ = 1)

        # `χ` is converted to `T` before the product, so both sides step by the same float. In
        # `Float32` a `Float64` `χ = 0.3` tells the two apart: the product rounded after the
        # conversion is a different float, and the two marches of this fixture diverge.
        @test march(T, plt, Ψ, efield, k, nt; dt = T(1e-4), χ = 0.3) ==
              march(T, plt, Ψ, efield, k, nt; dt = T(1e-4) * T(0.3), χ = 1)

        # The same pair from a zero position, where `zₓ` holds the step's own half-drifts and
        # nothing else, so that a last-ulp difference of the half-drift is not rounded away by the
        # sum into a `~0.5` position.
        zplt = ParticleList(zeros(T, 1, n), reshape(T(0.5) .* x, 1, n), reshape(w, 1, n))
        @test march(T, zplt, Ψ, efield, k, nt; dt = T(1e-4), χ = 0.3) ==
              march(T, zplt, Ψ, efield, k, nt; dt = T(1e-4) * T(0.3), χ = 1)
    end
end

@testset "the guards of the march, $T" for T in (Float64, Float32)
    f = march_fixture(T)
    march!(f, params) = reduced_integrate_vp(
        f.plt, f.Ψ, params, f.efield, f.SS, f.ip, f.ic, 1)

    # the accepted call: this is the reference the guards are read against, and the diagnostics it
    # saves are the ones of the state it wrote: the kinetic energy `½∑ⱼwⱼvⱼ²` and the momentum
    # `∑ⱼwⱼvⱼ` of the lifted step 0, from the weights the march was handed. `K` is summed in another
    # order in the march, so the two agree to the rounding of a sum rather than bitwise.
    @test march!(f, f.params) === nothing
    v0 = vec(f.SS.V[:, :, 1, 1])
    @test isapprox(f.SS.K[1, 1], sum(f.ic.w .* abs2.(v0)) / 2; rtol = sqrt(eps(T)))
    @test f.SS.K[1, 1] > zero(T)
    @test f.SS.M[1, 1] == dot(f.ic.w, v0)
    @test f.SS.M[1, 1] != zero(T)

    # `χ` is a positive timestep scale, and the comparison is written `> 0` so that `NaN` is
    # refused with the same exception
    for χ in (0.0, -1.0, NaN)
        @test_throws ArgumentError march!(f, (χ = χ,))
    end

    # the snapshots hold the `nₛ` saved times the parameters ask for
    @test_throws DimensionMismatch march!(
        march_fixture(T; nt = 4, ns = 5, nssnap = 4), (χ = 1,))

    # the scalar momentum needs one component
    @test_throws ArgumentError march!(march_fixture(T; nd = 2), (χ = 1,))

    # a basis of the wrong width has no reduced state of the length the cache holds: neither a
    # zero-width basis nor one with more columns than rows
    n = 40
    @test_throws ArgumentError reduced_integrate_vp(
        f.plt, zeros(T, n, 0), f.params, f.efield,
        f.SS, f.ip, f.ic, 1)
    @test_throws ArgumentError reduced_integrate_vp(
        f.plt, zeros(T, n, n + 1), f.params, f.efield,
        f.SS, f.ip, f.ic, 1)

    # the accepted side of the last two: `k = 1` with `nₚ = 1` runs
    g = march_fixture(T; n = 1, k = 1)
    @test reduced_integrate_vp(g.plt, g.Ψ, g.params, g.efield, g.SS, g.ip, g.ic, 1) ===
          nothing

    # the save is decided on the step index, so a march of a stride `nsave = 2` writes the same
    # states an every-step march writes at steps `0, 2, 4`: `nₜ = 4` with `nₛ = 3` saves `1:3`
    strided = march_fixture(T; nt = 4, ns = 3)
    @test reduced_integrate_vp(
        strided.plt, strided.Ψ, strided.params, strided.efield, strided.SS,
        strided.ip, strided.ic, 1) === nothing
    every = march_fixture(T; nt = 4, ns = 5)
    reduced_integrate_vp(
        every.plt, every.Ψ, every.params, every.efield, every.SS, every.ip, every.ic, 1)
    @test strided.SS.V[:, :, 1, 1] == every.SS.V[:, :, 1, 1]
    @test strided.SS.V[:, :, 2, 1] == every.SS.V[:, :, 3, 1]
    @test strided.SS.V[:, :, 3, 1] == every.SS.V[:, :, 5, 1]
    @test any(!iszero, strided.SS.V[:, :, 3, 1])
end

@testset "the integrator-parameters round trip, $T" for T in (Float64, Float32)
    IP = IntegratorParameters(T(1e-4), 10^4, 1001, 16, 1000, 6)
    @test length(IP.t) == 1001
    @test IP.t[1] == zero(T)
    @test IP.t[end] == T(10^4) * T(1e-4)

    mktempdir() do dir
        fpath = joinpath(dir, "integrator.h5")

        # the root group, which `h5save` finds and does not create
        h5open(fpath, "w") do h5
            h5save(h5, IP)
        end
        @test IntegratorParameters(fpath) == IP

        # a path that no group holds yet, which `h5save` makes with HDF5's own group creation
        h5open(fpath, "r+") do h5
            h5save(h5, IP; path = "integrator")
        end
        h5open(fpath, "r") do h5
            @test IntegratorParameters(h5, "integrator") == IP
        end

        # an `Int32` count is not an `Int`: every count of the parameters is `Int`, so the
        # conversion is refused rather than silently widened
        @test_throws MethodError IntegratorParameters(
            T(1e-4), Int32(10^4), 1001, 16, 1000, 6)

        # the HDF5 constructor checks the stored counts: `nₜ = 10` is not a multiple of `nₛ - 1 = 4`
        h5open(fpath, "r+") do h5
            g = create_group(h5, "bad")
            for (key, value) in (("dt", T(1e-4)), ("nt", 10), ("ns", 5), ("nh", 16),
                ("np", 1000), ("nparam", 6))
                attributes(g)[key] = value
            end
        end
        h5open(fpath, "r") do h5
            @test_throws ArgumentError IntegratorParameters(h5, "bad")
        end
    end
end

@testset "the counts of the integrator parameters, $T" for T in (Float64, Float32)
    # the smallest accepted sweep: one step, the initial state and the final one
    @test IntegratorParameters(T(1e-4), 1, 2, 16, 1000, 6).t == [zero(T), T(1e-4)]

    # `nₜ < 1` has no step to take; `nₛ = 2` keeps the `nₛ` check out of the way
    @test_throws ArgumentError IntegratorParameters(T(1e-4), 0, 2, 16, 1000, 6)
    @test_throws ArgumentError IntegratorParameters(T(1e-4), -10, 2, 16, 1000, 6)

    # `nₛ < 2` has no save after the initial state. At `dt = 0` the `range` of the times builds for
    # any length, so the count check is what raises
    @test_throws ArgumentError IntegratorParameters(zero(T), 10, 1, 16, 1000, 6)
    @test_throws ArgumentError IntegratorParameters(zero(T), 10, 0, 16, 1000, 6)

    # `nₜ = 10` is not a multiple of `nₛ - 1 = 4`, so the saves would not reach `nₛ`
    @test_throws ArgumentError IntegratorParameters(T(1e-4), 10, 5, 16, 1000, 6)
end

@testset "the step and the save allocate nothing, $T" for T in (Float64, Float32)
    n, k = 200, 20
    f = march_fixture(T; n = n, k = k, nt = 10)
    Δt, t, ts = T(1e-4), zero(T), 1

    @test @inferred(leapfrog_step!(f.ic, f.efield, Δt, t)) === f.ic
    @test @inferred(save_solution(f.SS, f.ic, f.Ψ, f.efield, 1, t, ts)) === nothing

    # coverage instrumentation allocates on every counted line, so under it the assertion does not
    # run rather than being weakened
    if Base.JLOptions().code_coverage == 0
        @test _alloc4(leapfrog_step!, f.ic, f.efield, Δt, t) == 0
        @test _alloc_save(save_solution, f.SS, f.ic, f.Ψ, f.efield, 1, t, ts) == 0
    end

    # the controls: the barrier sees an allocation when there is one. `similar` rather than
    # `copy`: on Julia 1.11 a `copy` whose result is dead is elided and reads zero, while
    # `similar` always allocates.
    @test _alloc4((a, b, c, d) -> similar(a.zₓ), f.ic, f.efield, Δt, t) > 0
    @test _alloc_save((a, b, c, d, e, g, h) -> similar(a.X), f.SS, f.ic, f.Ψ,
        f.efield, 1, t, ts) > 0
end

# The full-order march's own two barriers. `test/helpers/allocation_barriers.jl` carries the
# reduced march's `_alloc4` and `_alloc_save`, whose argument lists are those of
# `leapfrog_step!` and `save_solution`; the full march's buffers are separate arguments, so its
# step and its save take their own fixed-arity barrier.
function _alloc_full_step(f::F, x::X, v::V, a::A, w::W, e::E, Δt::D,
        t::TT) where {F, X, V, A, W, E, D, TT}
    (f(x, v, a, w, e, Δt, t); @allocated f(x, v, a, w, e, Δt, t))
end
function _alloc_full_save(f::F, SS::S, x::X, v::V, w::W, e::E, p::P, t::TT,
        ts::I) where {F, S, X, V, W, E, P, TT, I}
    (f(SS, x, v, w, e, p, t, ts); @allocated f(SS, x, v, w, e, p, t, ts))
end

"""
A seeded bump-on-tail list of `n` particles: the positions uniform on `(0, 2π/κ)` and the
velocities the two-population mix of `src/examples/bumpontail.jl`, with equal weights. The
potential is a periodic spline on the same domain, so the deposit covers it.
"""
function bump_on_tail(::Type{T}; n = 64, κ = T(0.3), α = T(0.1), σ = T(0.5), v₀ = T(4.5),
        seed = 7, ncells = 16) where {T}
    rng = Xoshiro(seed)
    x = (T(2π) / κ) .* rand(rng, T, n)
    v = sqrt(T(2)) .* randn(rng, T, n)
    for i in 1:n
        rand(rng) < α && (v[i] = σ * v[i] + v₀)
    end
    w = fill(T(1) / n, n)
    plt = ParticleList(reshape(x, 1, n), reshape(v, 1, n), reshape(w, 1, n))
    field = PoissonField(Potential(PeriodicBasisSpline((zero(T), T(2π) / κ), 3, ncells)))
    return plt, field
end

"""
The parts of one full-order march, so that one of them can be made wrong: the bump-on-tail list
`plt`, its `field`, the snapshots `SS`, the parameters `ip` and the caller's `params`, which
holds `χ`. `ns` defaults to `nₜ + 1`, which saves every step, and `nssnap` defaults to `ns`, so
that the saved times of the snapshots and the count the parameters ask for can disagree.
"""
function full_fixture(::Type{T}; n = 64, nt = 20, ns = nt + 1, nd = 1, dt = T(1e-3),
        χ = 1, seed = 7, nssnap = ns) where {T}
    plt, field = bump_on_tail(T; n = n, seed = seed)
    nh = length(coefficients(field))
    SS = Snapshots(T, nd, n, nh, nssnap, space(T))
    ip = IntegratorParameters(dt, nt, ns, nh, n, 1)
    return (; plt, field, SS, ip, params = (χ = χ,))
end

"The snapshots of one full-order march of `plt` on `field`, over `nt` steps saved as `ns`."
function full_march(::Type{T}, plt, field, nt, ns; dt = T(1e-3), χ = 1) where {T}
    n, nh = size(plt.x, 2), length(coefficients(field))
    SS = Snapshots(T, 1, n, nh, ns, space(T))
    ip = IntegratorParameters(dt, nt, ns, nh, n, 1)
    integrate_vp(plt, (χ = χ,), field, SS, ip, 1)
    return SS
end

@testset "check 2: the full-order march, $T" for T in (Float64, Float32)
    n, χ = 64, 2
    plt, field = bump_on_tail(T; n = n)
    x0, v0 = copy(plt.x), copy(plt.v)

    # (a) `χ` is the timestep scale: `(dt, χ = 2)` is `(2dt, χ = 1)` state for state
    @test full_march(T, plt, field, 20, 21; dt = T(1e-3), χ = χ) ==
          full_march(T, plt, field, 20, 21; dt = T(2e-3), χ = 1)

    # (b) the full march is the reduced march of the identity basis, bit for bit: a product with
    # the identity is exact in IEEE arithmetic
    Ψ = Matrix{T}(I, n, n)
    nh = length(coefficients(field))
    SSr = Snapshots(T, 1, n, nh, 21, space(T))
    IPr = IntegratorParameters(T(1e-3), 20, 21, nh, n, 1)
    ef = ReducedElectricField(field, Ψ, Ψ' * vec(plt.x))
    reduced_integrate_vp(
        plt, Ψ, (χ = χ,), ef, SSr, IPr, ReducedIntegratorCache(IPr, n), 1)
    full = full_march(T, plt, field, 20, 21; dt = T(1e-3), χ = χ)
    @test full.X == SSr.X
    @test full.V == SSr.V
    @test full.Φ == SSr.Φ
    @test full.W == SSr.W
    @test full.K == SSr.K
    @test full.M == SSr.M
    @test full == SSr

    # (c) slot 1 holds the initial state and slot `nₛ` the state after `IP.nₜ` steps: a march of
    # four steps saved as three holds the initial state, the state after two and that after four
    long = full_march(T, plt, field, 4, 3; dt = T(1e-3), χ = χ)
    @test long.X[:, :, 1, 1] == reshape(vec(x0), 1, n)
    @test long.V[:, :, 1, 1] == reshape(vec(v0), 1, n)
    short = full_march(T, plt, field, 2, 2; dt = T(1e-3), χ = χ)
    @test short.X[:, :, 2, 1] == long.X[:, :, 2, 1]
    @test short.V[:, :, 2, 1] == long.V[:, :, 2, 1]
    @test long.X[:, :, 3, 1] != long.X[:, :, 1, 1]

    # the march leaves the caller's particle list alone
    @test plt.x == x0 && plt.v == v0

    # (d) one step of the inner loop allocates nothing
    f = full_fixture(T)
    x, v, a, w = collect(vec(f.plt.x)), collect(vec(f.plt.v)), zeros(T, n),
    collect(vec(f.plt.w))
    Δt, t, ts = T(1e-3), zero(T), 1
    @test @inferred(full_leapfrog_step!(x, v, a, w, f.field, Δt, t)) === nothing
    @test @inferred(save_full_solution(f.SS, x, v, w, f.field, 1, t, ts)) === nothing
    if Base.JLOptions().code_coverage == 0
        @test _alloc_full_step(full_leapfrog_step!, x, v, a, w, f.field, Δt, t) == 0
        @test _alloc_full_save(save_full_solution, f.SS, x, v, w, f.field, 1, t, ts) == 0
    end
    @test _alloc_full_step((x, v, a, w, e, Δt, t) -> similar(x), x, v, a, w, f.field,
        Δt, t) > 0
    @test _alloc_full_save((SS, x, v, w, e, p, t, ts) -> similar(SS.X), f.SS, x, v, w,
        f.field, 1, t, ts) > 0
end

@testset "the marches take the scale of a ScaledField from params.χ, $T" for T in (Float64, Float32)
    n, nt, dt, χ = 64, 20, T(1e-3), 2
    plt, field = bump_on_tail(T; n = n)
    _, plain = bump_on_tail(T; n = n)
    nh = length(coefficients(field))

    # The reference holds no `ScaledField`: the leapfrog of the plain field at the step `dt·χ`, with
    # the acceleration `E/χ²` and the energy `W/χ²` written out.
    x, v, w = copy(vec(plt.x)), copy(vec(plt.v)), copy(vec(plt.w))
    a = zeros(T, n)
    Δt = dt * T(χ)
    for it in 1:nt
        x .+= (Δt / 2) .* v
        update!(plain, x, w, zero(T))
        efield!(plain, a, x)
        a ./= T(χ)^2
        v .+= Δt .* a
        x .+= (Δt / 2) .* v
    end
    update!(plain, x, w, zero(T))
    W = energy(plain) / T(χ)^2

    # a field built at the parameters' `χ` marches at that `χ`, whether the stored scale is the
    # parameters' integer or its own element type
    for stored in (χ, T(χ))
        scaled = ScaledField(field, stored)
        SS = Snapshots(T, 1, n, nh, 2, space(T))
        ip = IntegratorParameters(dt, nt, 2, nh, n, 1)
        integrate_vp(plt, (χ = χ,), scaled, SS, ip, 1)
        @test SS.X[1, :, 2, 1] == x
        @test SS.V[1, :, 2, 1] == v
        @test SS.W[2, 1] == W
    end

    # the scale reaches the field: the march of the plain field at the same `χ` differs
    SSp = full_march(T, plt, field, nt, 2; dt = dt, χ = χ)
    @test SSp.V[1, :, 2, 1] != v

    # the reduced march of the identity basis agrees with it bit for bit
    Ψ = Matrix{T}(I, n, n)
    ef = ReducedElectricField(ScaledField(field, T(χ)), Ψ, Ψ' * vec(plt.x))
    SSr = Snapshots(T, 1, n, nh, 2, space(T))
    IPr = IntegratorParameters(dt, nt, 2, nh, n, 1)
    reduced_integrate_vp(
        plt, Ψ, (χ = χ,), ef, SSr, IPr, ReducedIntegratorCache(IPr, n), 1)
    @test SSr.X[1, :, 2, 1] == x
    @test SSr.V[1, :, 2, 1] == v
    @test SSr.W[2, 1] == W
end

# The whole message of the exception a call raises, or a word when it raises nothing.
function _message(f)
    try
        f()
        return "no exception"
    catch thrown
        return sprint(showerror, thrown)
    end
end

@testset "a χ that a march would discard is refused, $T" for T in (Float64, Float32)
    # The step is `dt·params.χ` and the field is that of the `ScaledField`; a stored `χ` of its own
    # would put the two scales of one march out of step. The comparison is exact, so a stored `χ`
    # one ulp away is refused as well. Each refusal names both values, and the phrase is asserted
    # because `_check_march` and `_check_full_march` refuse a χ of the wrong sign with the same
    # `ArgumentError`: the type alone does not tell the guards apart.
    n, nt, dt, χ = 64, 20, T(1e-3), T(2)
    plt, field = bump_on_tail(T; n = n)
    nh = length(coefficients(field))

    for stored in (5, T(0.25), nextfloat(χ))
        scaled = ScaledField(field, stored)
        SS = Snapshots(T, 1, n, nh, 2, space(T))
        ip = IntegratorParameters(dt, nt, 2, nh, n, 1)
        @test_throws ArgumentError integrate_vp(plt, (χ = χ,), scaled, SS, ip, 1)
        msg = _message(() -> integrate_vp(plt, (χ = χ,), scaled, SS, ip, 1))
        @test occursin("the field is a ScaledField at χ", msg)
        @test !occursin("wraps", msg)
        @test occursin(string(stored), msg) && occursin(string(χ), msg)

        # the reduced path refuses it too, and with a phrase of its own: `_with_scale` of a
        # `ReducedElectricField` reads the scale of the field it wraps, so a reduced march names
        # the wrapped field's `χ` and not the stored `χ` of a scale it does not carry
        Ψ = Matrix{T}(I, n, n)
        ef = ReducedElectricField(scaled, Ψ, Ψ' * vec(plt.x))
        SSr = Snapshots(T, 1, n, nh, 2, space(T))
        IPr = IntegratorParameters(dt, nt, 2, nh, n, 1)
        ic = ReducedIntegratorCache(IPr, n)
        @test_throws ArgumentError reduced_integrate_vp(
            plt, Ψ, (χ = χ,), ef, SSr, IPr, ic, 1)
        msg_r = _message(() -> reduced_integrate_vp(
            plt, Ψ, (χ = χ,), ef, SSr, IPr, ic, 1))
        @test occursin("the reduced field wraps a ScaledField at χ", msg_r)
        @test occursin(string(stored), msg_r) && occursin(string(χ), msg_r)
    end
end

@testset "the guards of the full-order march, $T" for T in (Float64, Float32)
    f = full_fixture(T)
    march!(f, params) = integrate_vp(f.plt, params, f.field, f.SS, f.ip, 1)

    @test march!(f, f.params) === nothing

    # `χ` is a finite, positive timestep scale, and the comparison is written `> 0` so that `NaN`
    # is refused with the same exception
    for χ in (0, 0.0, -1, -1.0, NaN, Inf)
        @test_throws ArgumentError march!(f, (χ = χ,))
    end

    # the snapshots hold the `nₛ` saved times the parameters ask for
    @test_throws DimensionMismatch march!(
        full_fixture(T; nt = 4, ns = 5, nssnap = 4), (χ = 1,))

    # the scalar momentum needs one component
    @test_throws ArgumentError march!(full_fixture(T; nd = 2), (χ = 1,))

    # the particle count of the list is the particle extent of the snapshots, which the reduced
    # march has no guard for
    g = full_fixture(T)
    bad = Snapshots(T, 1, 63, length(coefficients(g.field)), g.ip.nₛ, space(T))
    @test_throws DimensionMismatch integrate_vp(g.plt, g.params, g.field, bad, g.ip, 1)

    # the accepted sides: the smallest march, and a march of one particle
    small = full_fixture(T; n = 1)
    @test integrate_vp(
        small.plt, small.params, small.field, small.SS, small.ip, 1) === nothing
end
