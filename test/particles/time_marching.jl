using HDF5
using LinearAlgebra
using ParticleMethods: ParticleList
using Random
using ReducedBasisMethods: get_DEIM_indices
using ReducedComplexityModeling: ParameterSpace, Parameter
using Test
# `import` and not `using`: the field interface and the march's internals are extended below.
import VlasovMethods: update!, efield!, energy, coefficients
using VlasovMethods
using VlasovMethods: h5save, leapfrog_step!, save_solution

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
    # The fixture is large enough to resolve the χ mutants at the ulp level: `Δt = T(IP.dt·params.χ)`
    # rounds the product after the conversion and differs from the shipped `IP.dt·T(params.χ)` by one
    # ulp of `Δt` at `χ = 0.3` (`3.0e-5` against `3.0000001e-5`), which is `Float32` round-off about
    # a `Float64` χ and so is invisible in `Float64`, where the two products are the same float. At
    # `n = 40` and `k = 4` the two marches agree in all 80080 saved entries and the mutant survives;
    # at `n = 1000` and `k = 20` — check 3's particle count and basis size — they first differ at
    # save 177 (`ReducedElectricField`) and 55 (`DEIMElectricField`), in 160229 and 99864 of the
    # 2002000 entries, and by up to `1.4e-1` and `3.2e-2` relative.
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
    end
end

@testset "the guards of the march, $T" for T in (Float64, Float32)
    f = march_fixture(T)
    march!(f, params) = reduced_integrate_vp(
        f.plt, f.Ψ, params, f.efield, f.SS, f.ip, f.ic, 1)

    # the accepted call: this is the reference the guards are read against, and the diagnostics it
    # saves are the ones of the state it wrote: the kinetic energy `½∑ⱼwⱼvⱼ²` and the momentum
    # `∑ⱼwⱼvⱼ` of the lifted step 0, from the weights the march was handed. `K` is a `dims`-free
    # reduction in the march, so the two agree to the rounding of a sum rather than bitwise.
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

    # `nₜ = 0` has no step to take. `nₛ = 2` keeps the `nₛ` guard out of the way, so that this one is
    # the guard that raises
    @test_throws ArgumentError march!(march_fixture(T; nt = 0, ns = 2), (χ = 1,))

    # `nₛ = 1` has no save after the initial state. The guard is reached only at `dt = 0`: a nonzero
    # `dt` with `nₛ = 1` throws from `IntegratorParameters` first, whose `range(0, stop = dt·nₜ,
    # length = 1)` has unequal endpoints. At `dt = 0` the parameters build and the march's own guard
    # is the one that raises
    @test_throws ArgumentError march!(march_fixture(T; ns = 1, dt = zero(T)), (χ = 1,))

    # `nₜ = 10` is not a multiple of `nₛ - 1 = 4`, so the saves would not reach `nₛ`
    @test_throws ArgumentError march!(march_fixture(T; nt = 10, ns = 5), (χ = 1,))

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

        # an `Int32` count is not an `Int`: the two counts of the parameters are `Int`, so the
        # conversion is refused rather than silently widened
        @test_throws MethodError IntegratorParameters(
            T(1e-4), Int32(10^4), 1001, 16, 1000, 6)
    end
end

@testset "the step and the save allocate nothing, $T" for T in (Float64, Float32)
    n, k = 200, 20
    f = march_fixture(T; n = n, k = k, nt = 10)
    Δt, t, it, nsave, save = T(1e-4), zero(T), 0, 1, true

    @test @inferred(leapfrog_step!(f.ic, f.efield, Δt, t)) === f.ic
    @test @inferred(save_solution(f.SS, f.ic, f.Ψ, f.efield, f.ic.w, 1, t, it, nsave,
        save)) === nothing

    # coverage instrumentation allocates on every counted line, so under it the assertion does not
    # run rather than being weakened
    if Base.JLOptions().code_coverage == 0
        @test _alloc_step(leapfrog_step!, f.ic, f.efield, Δt, t) == 0
        @test _alloc_save(
            save_solution, f.SS, f.ic, f.Ψ, f.efield, f.ic.w, 1, t, it, nsave,
            save) == 0
    end

    # the controls: the barrier sees an allocation when there is one. `similar` rather than
    # `copy`: on Julia 1.11 a `copy` whose result is dead is elided and reads zero, while
    # `similar` always allocates.
    @test _alloc_step((a, b, c, d) -> similar(a.zₓ), f.ic, f.efield, Δt, t) > 0
    @test _alloc_save((a, b, c, d, e, g, h, i, j, k) -> similar(a.X), f.SS, f.ic, f.Ψ,
        f.efield, f.ic.w, 1, t, it, nsave, save) > 0
end
