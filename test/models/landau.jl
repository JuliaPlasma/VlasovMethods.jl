using LinearAlgebra
using ParticleMethods: ParticleList
using Random
using SimpleSplines: BSplineBasis, UniformMesh, Free, TensorProductBasis,
                     TensorProductQuadrature
using StaticArrays: SMatrix, SVector
using Test
using VlasovMethods

const VM = VlasovMethods

include("../helpers/allocation_barriers.jl")

# `n1²` particles of uniform weight at the cell centres of an `n1 × n1` grid on the velocity
# domain `[-2, 2]²`, so that the L² projection stays positive (KNOWN_ISSUES.md, K8), and a cubic
# clamped basis of 8 × 8 cells: M = 121 basis functions, Q = 1600 quadrature points.
function landau_setup(::Type{T}, n1; γ = -3) where {T}
    npart = n1^2
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pdist = ParticleDistribution(1, 2, ParticleList(zeros(T, 4, npart); variables = vars))
    g = [T(-2) + (T(i) - T(1) / 2) * T(4) / n1 for i in 1:n1]
    pdist.particles.v[1, :] .= repeat(g; outer = n1)
    pdist.particles.v[2, :] .= repeat(g; inner = n1)
    pdist.particles.w[1, :] .= one(T) / npart
    axis = BSplineBasis(UniformMesh(8, T(-2) .. T(2)), 3, Free())
    B = TensorProductBasis((axis, axis))
    sdist = SplineDistribution(1, 2, B, TensorProductQuadrature(B), zeros(T, size(B)...))
    return Landau(pdist, CollisionEntropy(sdist); γ), pdist
end

# The bytes of one call of three arguments after a warm call, through a barrier whose arguments
# have concrete types; the shared helper has no barrier of this arity.
function bytes3(f::F, a::A, b::B, c::C) where {F, A, B, C}
    f(a, b, c)
    return @allocated f(a, b, c)
end

# The controls: calls of each arity that allocate.
copy_first(a, b, c, d, e) = copy(a)
copy_first(a, b, c, d, e, g) = copy(a)
# A least-squares solve of the vector field, on the four-argument barrier of the vector field.
solve(K, LJ, _, _) = K \ LJ

@testset "Landau buffers, $T, Np = $(n1^2), γ = $γ" for T in (Float32, Float64),
    n1 in (16, 24), γ in (-3, 0)
    landau, pdist = landau_setup(T, n1; γ)
    @test landau.γ == γ
    cache = landau.cache[T]
    sdist = cache.sdist
    @test length(sdist) == 121
    @test size(cache.D[1]) == (121, 1600)
    v = copy(pdist.particles.v)
    v̇ = similar(v)

    vf = _alloc4(VM.collisional_vectorfield!, v̇, v, nothing, landau)
    # the two solves on the matrices that the evaluation leaves in the cache
    ceiling = _alloc4(solve, cache.K1, cache.LJ, nothing, nothing) +
              _alloc4(solve, cache.K2, cache.LJ, nothing, nothing)

    @test @inferred(VM.collisional_vectorfield!(v̇, v, nothing, landau)) === nothing
    @test @inferred(VM.compute_L!(cache.L, sdist, landau)) === cache.L
    @test @inferred(VM.compute_J!(cache.J, sdist, landau)) === cache.J
    @test @inferred(VM.compute_K!(cache.K1, cache.K2, v, sdist, landau)) ==
          (cache.K1, cache.K2)
    @test @inferred(VM._projection!(
        cache.rhs, cache.vals, cache.mwork, v, pdist, sdist)) ===
          sdist.spline

    # the vector field is the minimum-norm solution of the two least-squares problems
    @test isapprox(v̇[1, :], (cache.K1 \ cache.LJ) .* landau.ν; rtol = 8 * eps(T))
    @test isapprox(v̇[2, :], (cache.K2 \ cache.LJ) .* landau.ν; rtol = 8 * eps(T))

    # coverage instrumentation allocates on every counted line, so under it the assertions
    # do not run rather than being weakened
    if Base.JLOptions().code_coverage == 0
        @test bytes3(VM.compute_L!, cache.L, sdist, landau) == 0
        @test bytes3(VM.compute_J!, cache.J, sdist, landau) == 0
        @test _alloc5(VM.compute_K!, cache.K1, cache.K2, v, sdist, landau) == 0
        @test _alloc6(
            VM._projection!, cache.rhs, cache.vals, cache.mwork, v, pdist, sdist) ==
              0
        @test 0 < vf ≤ ceiling
    end
    # the controls: each barrier sees an allocation
    @test bytes3(projection, v, pdist, sdist) > 0
    @test ceiling > 0
    @test _alloc5(copy_first, v, v, v, v, v) > 0
    @test _alloc6(copy_first, v, v, v, v, v, v) > 0
end

@testset "the in-place 2D projection is the allocating one" begin
    landau, pdist = landau_setup(Float64, 16)
    cache = landau.cache[Float64]
    v = copy(pdist.particles.v)
    v[1, 1:7] .+= 0.1
    projection(v, pdist, cache.sdist)
    expected = copy(cache.sdist.coefficients)
    fill!(cache.sdist.coefficients, 0)
    VM._projection!(cache.rhs, cache.vals, cache.mwork, v, pdist, cache.sdist)
    @test cache.sdist.coefficients == expected
    # the load array is cleared on each call, so a second call gives the same coefficients
    VM._projection!(cache.rhs, cache.vals, cache.mwork, v, pdist, cache.sdist)
    @test cache.sdist.coefficients == expected
end

@testset "compute_L! reads the gradient tabulations stored in the cache" begin
    landau, pdist = landau_setup(Float64, 16)
    cache = landau.cache[Float64]
    sdist = cache.sdist
    projection(pdist.particles.v, pdist, sdist)

    # `LandauCache` stores the tabulations of its own quadrature grid
    @test cache.D == map(Matrix, VM.gradient_tabulations(sdist))

    L = VM.compute_L!(zeros(121, 121), sdist, landau)
    # `L` is quadratic in the tabulations, and a factor of 2 is exact in floating point, so
    # doubling the stored tables gives exactly 4L when `compute_L!` reads them
    foreach(D -> D .*= 2, cache.D)
    @test VM.compute_L!(zeros(121, 121), sdist, landau) == 4 .* L
    foreach(D -> D ./= 2, cache.D)
    @test VM.compute_L!(zeros(121, 121), sdist, landau) == L
end

@testset "compute_L! rescales the values of X in place and keeps its sparsity pattern" begin
    landau, pdist = landau_setup(Float64, 16)
    cache = landau.cache[Float64]
    sdist = cache.sdist
    projection(pdist.particles.v, pdist, sdist)
    D = map(copy, cache.D)

    VM.compute_L!(zeros(121, 121), sdist, landau)
    for c in 1:2
        # the same structure as `D`, and the values of `D` scaled by the weighted density
        @test cache.X[c].colptr == D[c].colptr
        @test cache.X[c].rowval == D[c].rowval
        @test cache.D[c] == D[c]
        @test cache.X[c] == D[c] * Diagonal(cache.s)
    end

    # a second call with a different density rescales the values again, on the same pattern
    colptr = map(X -> copy(X.colptr), cache.X)
    rowval = map(X -> copy(X.rowval), cache.X)
    sdist.coefficients .*= 2
    VM.compute_L!(zeros(121, 121), sdist, landau)
    for c in 1:2
        @test cache.X[c].colptr == colptr[c]
        @test cache.X[c].rowval == rowval[c]
        @test cache.X[c] == D[c] * Diagonal(cache.s)
    end
end

@testset "a Landau on a one-dimensional velocity space does not build" begin
    # `LandauCache` builds the buffers of a two-dimensional velocity space, and no method of
    # `compute_L!` takes another, so the constructor raises a `MethodError`
    pdist1 = ParticleDistribution(1, 1, 10)
    sdist1 = SplineDistribution(1, 1, 21, 4, (-6.0, 6.0), 0, Free())
    @test_throws MethodError Landau(pdist1, CollisionEntropy(sdist1))
end

@testset "compute_L! refuses a spline distribution on another quadrature" begin
    pd = ParticleDistribution(1, 2, 100)
    g = v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π) + 0.05
    sd1 = SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free())
    sd2 = SplineDistribution(1, 2, 4, 2, (-2.5, 2.5), 0, Free())
    project_function(g, sd1)
    project_function(g, sd2)
    landau1 = Landau(pd, CollisionEntropy(sd1))
    M = length(sd1)
    @test length(sd2) == M

    # the same sizes on another domain: the stored tabulations are not those of `sd2`'s grid
    err = try
        VM.compute_L!(zeros(M, M), sd2, landau1)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test startswith(err.msg, "compute_L!: `sdist` does not share the quadrature")

    # `similar` keeps the quadrature, so a copy of the cache's distribution is accepted
    sd1b = similar(Float64, sd1)
    sd1b.coefficients .= sd1.coefficients
    @test VM.compute_L!(zeros(M, M), sd1b, landau1) ==
          VM.compute_L!(zeros(M, M), sd1, landau1)

    # the same basis with another quadrature: sharing the basis is not enough
    sd3 = SplineDistribution(1, 2, sd1.basis, TensorProductQuadrature(sd1.basis; nq = 4),
        copy(sd1.coefficients))
    @test_throws ArgumentError VM.compute_L!(zeros(M, M), sd3, landau1)
end

@testset "compute_J! refuses, and compute_K! takes, a spline distribution on another quadrature" begin
    pd = ParticleDistribution(1, 2, 100)
    pd.particles.v .= [2 * sin(i * α) for i in 1:2, α in 1:100]
    pd.particles.w .= 1 / 100
    g = v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π) + 0.05
    sd1 = SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free())
    project_function(g, sd1)
    landau1 = Landau(pd, CollisionEntropy(sd1))
    M1 = length(sd1)
    N = size(pd.particles.v, 2)
    refuses_J(sd) =
        try
            VM.compute_J!(zeros(length(sd)), sd, landau1)
            nothing
        catch e
            e isa ArgumentError &&
                startswith(e.msg, "compute_J!: `sdist` does not share the quadrature")
        end

    # the same sizes on another domain, and another number of cells
    foreign = (SplineDistribution(1, 2, 4, 2, (-2.5, 2.5), 0, Free()),
        SplineDistribution(1, 2, 6, 2, (-3.0, 3.0), 0, Free()))
    for sd in foreign
        project_function(g, sd)
        @test refuses_J(sd) === true
    end

    # `similar` keeps the quadrature, so a copy of the cache's distribution is accepted
    sd1b = similar(Float64, sd1)
    sd1b.coefficients .= sd1.coefficients
    @test VM.compute_J!(zeros(M1), sd1b, landau1) == VM.compute_J!(zeros(M1), sd1, landau1)

    # the same basis with another quadrature: sharing the basis is not enough
    sd3 = SplineDistribution(1, 2, sd1.basis, TensorProductQuadrature(sd1.basis; nq = 4),
        copy(sd1.coefficients))
    @test refuses_J(sd3) === true

    # `compute_K!` reads no quadrature, so it takes any `sdist` of the cache's basis order and
    # gives the `K` of a `Landau` built on that distribution
    for sd in (foreign..., sd3)
        M = length(sd)
        K = VM.compute_K!(zeros(M, N), zeros(M, N), pd.particles.v, sd, landau1)
        landau = Landau(pd, CollisionEntropy(sd))
        @test K == VM.compute_K!(zeros(M, N), zeros(M, N), pd.particles.v, sd, landau)
        @test any(!iszero, K[1]) && any(!iszero, K[2])
    end
end

# A `Landau` of the exponent `γ` on a small basis, for the tests that read only its kernel.
function kernel_model(γ)
    pd = ParticleDistribution(1, 2, 100)
    sd = SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free())
    return Landau(pd, CollisionEntropy(sd); γ)
end

@testset "the kernel, $T" for T in (Float32, Float64)
    # 25 directions, and radii from 1e-2 to 1e2
    us = [SVector{2, T}(cos(2π * k / 25 + 0.1) * 10.0^((k - 12) / 6),
              sin(2π * k / 25 + 0.1) * 10.0^((k - 12) / 6)) for k in 0:24]
    vβ = zero(SVector{2, T})
    γs = (-5.5, -3, -1.5, 0, 1)

    # the Coulomb kernel of γ = -3, bit for bit, in the operations of the Coulomb formula
    coulomb = kernel_model(-3)
    @test coulomb.γ == -3
    @test Landau(ParticleDistribution(1, 2, 10),
        CollisionEntropy(SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free()))).γ == -3
    for u in us
        u1, u2 = u
        n2 = u1 * u1 + u2 * u2
        inv_n = inv(sqrt(n2))
        inv_n3 = inv_n / n2
        U11 = inv_n - u1 * u1 * inv_n3
        U12 = -u1 * u2 * inv_n3
        U22 = inv_n - u2 * u2 * inv_n3
        @test VM.kernel(u, vβ, coulomb) == SMatrix{2, 2, T}(U11, U12, U12, U22)
    end

    # Maxwell molecules, γ = 0: |u|² I − u uᵀ
    maxwell = kernel_model(0)
    for u in us
        u1, u2 = u
        n2 = u1 * u1 + u2 * u2
        U = VM.kernel(u, vβ, maxwell)
        @test maximum(abs, U - SMatrix{2, 2, T}(n2 - u1 * u1, -u1 * u2, -u1 * u2, n2 -
                                                                              u2 * u2)) ≤
              4 * eps(T) * n2
    end

    for γ in γs
        landau = kernel_model(γ)
        # the zero matrix at u = 0, in the element type of the velocities
        @test VM.kernel(vβ, vβ, landau) === zero(SMatrix{2, 2, T})
        for u in us
            U = @inferred VM.kernel(u + vβ, vβ, landau)
            @test eltype(U) === T
            @test U[1, 2] == U[2, 1]
            # the norm of u, and the expected values, in BigFloat
            n = sqrt(big(u[1])^2 + big(u[2])^2)
            # U u = 0, to the size |u|^{γ+2} |u| of the terms that cancel
            @test maximum(abs, big.(U) * big.(u)) ≤ 8 * eps(T) * n^(γ + 3)
            # the eigenvalues are 0 on u and |u|^{γ+2} on its complement
            λ = maximum(eigvals(Symmetric(Matrix(U))))
            @test abs(λ - n^(γ + 2)) ≤ 8 * eps(T) * n^(γ + 2)
        end
    end

    # radii where |u|³ or 1/|u|³ leaves the range of T: for γ ≠ -3 the kernel is finite and
    # accurate wherever its largest eigenvalue |u|^{γ+2} is a normal number of T
    small, large = cbrt(floatmin(T)) / 16, cbrt(floatmax(T)) * 16
    checked = 0
    for γ in γs, r in (small, large), k in 0:4
        γ == -3 && continue
        u = SVector{2, T}(r * cos(2π * k / 5 + 0.1), r * sin(2π * k / 5 + 0.1))
        λ = sqrt(big(u[1])^2 + big(u[2])^2)^(γ + 2)
        floatmin(T) / eps(T) ≤ λ ≤ floatmax(T) / 4 || continue
        U = VM.kernel(u, vβ, kernel_model(γ))
        @test all(isfinite, U)
        @test abs(maximum(eigvals(Symmetric(Matrix(U)))) - λ) ≤ 8 * eps(T) * λ
        checked += 1
    end
    # γ = -1.5 and γ = 0, each at both radii in five directions
    @test checked == 20
end

@testset "refused exponents" begin
    pd = ParticleDistribution(1, 2, 100)
    ent = CollisionEntropy(SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free()))
    for γ in (-6, -7, -Inf, Inf, NaN)
        err = try
            Landau(pd, ent; γ)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test err isa ArgumentError && occursin("γ > -6", err.msg)
    end
    # the control: one exponent just above the bound, and the most singular one of the tests below
    @test Landau(pd, ent; γ = -5.9).γ == -5.9
    @test Landau(pd, ent; γ = -5.5).γ == -5.5
end

# The kernel |u|^{γ+2} (I − u uᵀ/|u|²), written out apart from the package, with the value `U0`
# at u = 0.
function landau_kernel(u, γ, U0 = zeros(2, 2))
    n = norm(u)
    iszero(n) && return U0
    return n^(γ + 2) * (I - u * u' / n^2)
end

# The grid points, the weights s_a = W_a |f_s(v_a)| and the gradient tables D^c of `sd`'s
# quadrature, the grid enumerated with the first axis fastest, as `compute_L!` enumerates it.
function landau_grid(sd)
    X = quadrature_nodes(sd.quadrature)
    W = quadrature_weights(sd.quadrature)
    Q1, Q2 = length(X[1]), length(X[2])
    pts = [[X[1][a1], X[2][a2]] for a2 in 1:Q2 for a1 in 1:Q1]
    s = [W[1][a1] * W[2][a2] * abs(sd.spline((X[1][a1], X[2][a2])))
         for a2 in 1:Q2 for a1 in 1:Q1]
    return pts, s, map(Matrix, VM.gradient_tabulations(sd))
end

# 𝕃_ij = ½ Σ_ab s_a s_b Σ_cd (∂_c φ_i(v_a) − ∂_c φ_i(v_b)) U^cd_ab (∂_d φ_j(v_a) − ∂_d φ_j(v_b)).
function landau_double_sum(sd, γ)
    pts, s, D = landau_grid(sd)
    M = size(D[1], 1)
    L = zeros(M, M)
    for a in eachindex(pts), b in eachindex(pts)

        U = landau_kernel(pts[a] - pts[b], γ)
        for c in 1:2, d in 1:2

            ga = D[c][:, a] - D[c][:, b]
            gb = D[d][:, a] - D[d][:, b]
            L .+= (s[a] * s[b] * U[c, d] / 2) .* (ga * gb')
        end
    end
    return L
end

# 𝕃 in the expanded form of `compute_L!`: the kernel sums A^cd_a = Σ_b s_b U^cd_ab with a
# diagonal weighting, minus the coupled term. `Usum` and `Ucoupled` give the kernel of the
# grid points a and b in the first and in the second sum.
function landau_expanded(sd, Usum, Ucoupled)
    pts, s, D = landau_grid(sd)
    Q = length(pts)
    M = size(D[1], 1)
    L = zeros(M, M)
    for c in 1:2, d in 1:2

        P = [Usum(pts, a, b)[c, d] for a in 1:Q, b in 1:Q]
        C = [Ucoupled(pts, a, b)[c, d] for a in 1:Q, b in 1:Q]
        L .+= D[c] * Diagonal(s .* (P * s)) * D[d]' -
              (D[c] * Diagonal(s)) * C * (D[d] * Diagonal(s))'
    end
    return L
end

# The fixture of "𝕃 is positive semi-definite on a sign-changing projection" of
# `test/distributions/spline_distribution.jl`: 400 normal particles of seed 0xABBA, 7 knots of
# order 3 (quadratic splines) on [-4, 4]², and the projection of the unit Gaussian, which is
# negative at a node.
function quadratic_fixture(γ)
    Random.seed!(0xABBA)
    pd = ParticleDistribution(1, 2, 400)
    pd.particles.v .= randn(2, 400)
    pd.particles.w[1, :] .= 1 / 400
    sd = SplineDistribution(1, 2, 7, 3, (-4.0, 4.0), 0, Free())
    project_function(v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π), sd)
    return sd, Landau(pd, CollisionEntropy(sd); γ)
end

# The fixture of "L equals the defining double sum" of the same file: 200 normal particles of
# seed 0x51DEA1, scaled by 1/2, projected on 4 knots of order 2 (linear splines) on [-3, 3]².
function linear_fixture(γ)
    sd = SplineDistribution(1, 2, 4, 2, (-3.0, 3.0), 0, Free())
    Random.seed!(0x51DEA1)
    pd = ParticleDistribution(1, 2, 200)
    pd.particles.v .= randn(2, 200) ./ 2
    pd.particles.w[1, :] .= 1 / 200
    projection(pd.particles.v, pd, sd)
    return sd, Landau(pd, CollisionEntropy(sd); γ)
end

@testset "𝕃 for each γ" begin
    for γ in (-5.5, -3, 0), fixture in (quadratic_fixture, linear_fixture)

        sd, landau = fixture(γ)
        M = length(sd)
        L = VM.compute_L!(zeros(M, M), sd, landau)
        scale = maximum(abs, L)
        @test scale > 0
        @test maximum(abs, L - L') ≤ 1e-10 * scale
        @test maximum(abs, L * ones(M)) ≤ 1e-10 * scale
        @test minimum(eigvals(Symmetric(L))) ≥ -1e-12 * scale
        if fixture === linear_fixture
            Lslow = landau_double_sum(sd, γ)
            @test maximum(abs, L - Lslow) < 1e-10 * maximum(abs, Lslow)
        end
    end
    # the projection of the quadratic fixture changes sign at a node
    sd, _ = quadratic_fixture(-3)
    X = quadrature_nodes(sd.quadrature)
    @test minimum(sd.spline((a, b)) for a in X[1], b in X[2]) < 0
end

@testset "the coincident value does not matter" begin
    for γ in (-3, 0)
        sd, landau = linear_fixture(γ)
        M = length(sd)
        L = VM.compute_L!(zeros(M, M), sd, landau)
        scale = maximum(abs, L)
        # the coincident pairs are a = b: the points of one grid are distinct
        zero_at(pts, a, b) = landau_kernel(pts[a] - pts[b], γ)
        seven_at(pts, a, b) = landau_kernel(pts[a] - pts[b], γ, 7.0 * I(2))
        L0 = landau_expanded(sd, zero_at, zero_at)
        L7 = landau_expanded(sd, seven_at, seven_at)
        @test maximum(abs, L0 - L7) ≤ 1e-12 * scale
        @test maximum(abs, L0 - L) ≤ 1e-10 * scale
        @test maximum(abs, L7 - L) ≤ 1e-10 * scale
        # the control: the value 7 I in the kernel sums alone changes 𝕃, so the two sums cancel
        # the coincident terms only together
        @test maximum(abs, landau_expanded(sd, seven_at, zero_at) - L) > 1e-3 * scale
    end
end
