using Random
using Test
using VlasovMethods: _apply_∫dv!, _apply_∫dvᵀ!, _apply_∫dv_μ!, _apply_∫dvᵀ_μ!,
                     _apply_∫vdv!, _apply_∫v²dv!

include("../helpers/allocation_barriers.jl")

# The six velocity moments of a grid distribution are inferred and allocate nothing, on two grids.
# Each reading follows a warm call and an `@inferred` call at the same argument types, and is the
# minimum of three readings through a barrier of fixed arity. Coverage instrumentation allocates
# on every counted line, so the assertions do not run under it rather than being weakened.
@testset "the moments are inferred and allocate nothing, $T, $nx × $nv" for T in (Float32, Float64),
    (nx, nv) in ((16, 64), (64, 256))

    rng = Xoshiro(11)
    ci = CartesianIndices((nx, nv))
    li = LinearIndices((nx, nv))
    h₁, h₂ = one(T) / nx, T(6) / (nv - 1)
    v = collect(range(T(-3), T(3); length = nv))
    μ = rand(rng, T, nv)
    f = rand(rng, T, nx * nv)
    ψ = rand(rng, T, nx)
    y = zeros(T, nx)
    Y = zeros(T, nx * nv)
    F = reshape(f, nx, nv)
    tol = sqrt(eps(T))

    # each moment against its sum over the grid, written out here
    @test @inferred(_apply_∫dv!(y, f, ci, li, h₁, h₂)) === nothing
    @test y ≈ [sum(F[i, j] * h₂ for j in 1:nv) for i in 1:nx] rtol = tol
    @test @inferred(_apply_∫dv_μ!(y, f, μ, ci, li, h₁, h₂)) === nothing
    @test y ≈ [sum(F[i, j] * h₂ * μ[j] for j in 1:nv) for i in 1:nx] rtol = tol
    @test @inferred(_apply_∫vdv!(y, f, v, ci, li, h₁, h₂)) === nothing
    @test y ≈ [sum(F[i, j] * h₂ * v[j] for j in 1:nv) for i in 1:nx] rtol = tol
    @test @inferred(_apply_∫v²dv!(y, f, v, ci, li, h₁, h₂)) === nothing
    @test y ≈ [sum(F[i, j] * h₂ * v[j]^2 for j in 1:nv) for i in 1:nx] rtol = tol
    @test @inferred(_apply_∫dvᵀ!(Y, ψ, ci, li, h₁, h₂)) === nothing
    @test Y == vec([ψ[i] * h₂ for i in 1:nx, j in 1:nv])
    @test @inferred(_apply_∫dvᵀ_μ!(Y, ψ, μ, ci, li, h₁, h₂)) === nothing
    @test Y == vec([ψ[i] * h₂ * μ[j] for i in 1:nx, j in 1:nv])

    if Base.JLOptions().code_coverage == 0
        @test minimum(_ -> _alloc6(_apply_∫dv!, y, f, ci, li, h₁, h₂), 1:3) == 0
        @test minimum(_ -> _alloc6(_apply_∫dvᵀ!, Y, ψ, ci, li, h₁, h₂), 1:3) == 0
        @test minimum(_ -> _alloc_save(_apply_∫dv_μ!, y, f, μ, ci, li, h₁, h₂), 1:3) == 0
        @test minimum(_ -> _alloc_save(_apply_∫dvᵀ_μ!, Y, ψ, μ, ci, li, h₁, h₂), 1:3) == 0
        @test minimum(_ -> _alloc_save(_apply_∫vdv!, y, f, v, ci, li, h₁, h₂), 1:3) == 0
        @test minimum(_ -> _alloc_save(_apply_∫v²dv!, y, f, v, ci, li, h₁, h₂), 1:3) == 0
    end
    # the controls: each barrier sees an allocation, so the zeros above are not vacuous
    @test _alloc6((y, f, ci, li, h₁, h₂) -> copy(f), y, f, ci, li, h₁, h₂) > 0
    @test _alloc_save((y, f, v, ci, li, h₁, h₂) -> copy(f), y, f, v, ci, li, h₁, h₂) > 0
end
