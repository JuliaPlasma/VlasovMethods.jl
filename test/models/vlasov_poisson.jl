using Random
using Test
using VlasovMethods: s_advection!, v_advection!

include("../helpers/allocation_barriers.jl")

# The advection half of the particle splitting of `VlasovPoisson`: its vector field and its
# solution are inferred and allocate nothing, at `N = 256` particles. Each reading follows a warm
# call and an `@inferred` call at the same argument types, and is the minimum of three readings
# through a barrier of fixed arity. Coverage instrumentation allocates on every counted line, so
# the assertions do not run under it rather than being weakened.
@testset "the advection is inferred and allocates nothing, $T" for T in (Float32, Float64)
    N = 256
    z̄ = randn(Xoshiro(7), T, 2, N)
    z = similar(z̄)
    ż = similar(z̄)
    params = (w = fill(one(T) / N, 1, N),)
    t̄, t = T(0.3), T(0.4)

    v_advection!(ż, t, z̄, params)
    s_advection!(z, t, z̄, t̄, params)
    @test @inferred(v_advection!(ż, t, z̄, params)) === nothing
    @test @inferred(s_advection!(z, t, z̄, t̄, params)) === nothing

    # the values: `ẋ = v`, `v̇ = 0`, and `x ← x̄ + (t - t̄) v̄`, `v` unchanged
    @test ż[1, :] == z̄[2, :]
    @test all(iszero, ż[2, :])
    @test z[1, :] == z̄[1, :] .+ (t - t̄) .* z̄[2, :]
    @test z[2, :] == z̄[2, :]

    if Base.JLOptions().code_coverage == 0
        @test minimum(_ -> _alloc4(v_advection!, ż, t, z̄, params), 1:3) == 0
        @test minimum(_ -> _alloc5(s_advection!, z, t, z̄, t̄, params), 1:3) == 0
    end
    # the controls: each barrier sees an allocation, so the zeros above are not vacuous
    @test _alloc4((ż, t, z, p) -> copy(z), ż, t, z̄, params) > 0
    @test _alloc5((z, t, z̄, t̄, p) -> copy(z̄), z, t, z̄, t̄, params) > 0
end
