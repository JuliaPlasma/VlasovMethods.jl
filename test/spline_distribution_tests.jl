using VlasovMethods
using LinearAlgebra
using Random
using Test

const VM = VlasovMethods

@testset "SplineDistribution construction" begin
    # 1-D velocity space, cubic (order 4) as both manuscripts specify, on [-10, 10]
    sd = SplineDistribution(1, 1, 41, 4, (-10.0, 10.0), 0, Free())
    @test ndims(sd) == 1
    @test length(sd) == 40 + 3               # n + p = 40 cells + degree 3
    @test polynomial_reproduction(sd.basis) == 3
    @test eltype(sd) == Float64

    # the boundary conditions that break conservation report it
    @test polynomial_reproduction(
        SplineDistribution(1, 1, 41, 4, (-10.0, 10.0), 0, Periodic()).basis) == 0
    @test polynomial_reproduction(
        SplineDistribution(1, 1, 41, 4, (-10.0, 10.0), 0, Dirichlet()).basis) == -1

    # the old symbols still map to the same spaces they used to select
    @test polynomial_reproduction(
        SplineDistribution(1, 1, 21, 3, (-5.0, 5.0), 0, :nothing).basis) == 2
    @test polynomial_reproduction(
        SplineDistribution(1, 1, 21, 3, (-5.0, 5.0), 0, :Dirichlet).basis) == -1
    @test polynomial_reproduction(
        SplineDistribution(1, 1, 21, 3, (-5.0, 5.0), 0, :Periodic).basis) == 0
    # ...and a typo is an error rather than a silent fall-through to the clamped basis
    @test_throws ArgumentError SplineDistribution(1, 1, 21, 3, (-5.0, 5.0), 0, :Perodic)

    # oversized end cells
    sb = SplineDistribution(1, 1, 21, 4, (-5.0, 5.0), 3.0, Free())
    @test VM.domain(sb.basis) == (-8.0 .. 8.0)
    @test ncells(sb.basis) == 22            # 20 interior + one at each end

    # 2-D velocity space: a tensor product
    s2 = SplineDistribution(1, 2, 11, 3, (-4.0, 4.0), 0, Free())
    @test ndims(s2) == 2
    @test size(s2) == (12, 12)
    @test length(s2) == 144
    @test s2.basis isa TensorProductBasis
    @test polynomial_reproduction(s2.basis) == 2

    # the lumped-mass flag is gone rather than silently reinterpreted
    @test_throws MethodError SplineDistribution(1, 1, 21, 3, (-5.0, 5.0), 0, Free(), false)
end

@testset "particle projection reproduces moments" begin
    npart = 20000
    Random.seed!(0xC0FFEE)

    pd = ParticleDistribution(1, 1, npart)
    pd.particles.v[1, :] .= randn(npart)
    pd.particles.x[1, :] .= rand(npart)
    pd.particles.w[1, :] .= 1 / npart

    sd = SplineDistribution(1, 1, 41, 4, (-10.0, 10.0), 0, Free())
    v = pd.particles.v[1, :]
    fs = projection(v, pd, sd)

    b = sd.basis
    q = sd.quadrature
    Φ = basis_values(q, 0)
    w = quadrature_weights(q)

    # ∫ f_s dv, ∫ v f_s dv, ∫ v² f_s dv against the particle moments. This is the property
    # that makes the schemes conserve: 1, v, v² are in the span of a clamped cubic basis, so
    # the L² projection reproduces all three exactly (up to quadrature).
    x = quadrature_nodes(q)
    fvals = fs.(x)
    m0 = sum(w .* fvals)
    m1 = sum(w .* x .* fvals)
    m2 = sum(w .* x .^ 2 .* fvals)

    p0 = sum(pd.particles.w[1, :])
    p1 = sum(pd.particles.w[1, :] .* v)
    p2 = sum(pd.particles.w[1, :] .* v .^ 2)

    @test isapprox(m0, p0; rtol = 1e-10)
    @test isapprox(m1, p1; atol = 1e-10)
    @test isapprox(m2, p2; rtol = 1e-10)

    # A Dirichlet basis cannot represent the constants at all, so the mass is *not* an
    # invariant of the projection. How large the error is depends on how much of the
    # distribution reaches the boundary: for a unit normal on [-10,10] it is only ~1e-8,
    # because f is already negligible where the constraint bites. The loss is structural
    # rather than proportional, so a fair demonstration needs support that reaches the edge.
    npu = 20000
    pu = ParticleDistribution(1, 1, npu)
    pu.particles.v[1, :] .= range(-4.9, 4.9; length = npu)     # fills the domain
    pu.particles.w[1, :] .= 1 / npu
    vu = pu.particles.v[1, :]

    for (bc, tol) in ((Free(), 1e-10), (Dirichlet(), nothing))
        s = SplineDistribution(1, 1, 41, 4, (-5.0, 5.0), 0, bc)
        f = projection(vu, pu, s)
        qq = s.quadrature
        m = sum(quadrature_weights(qq) .* f.(quadrature_nodes(qq)))
        if tol === nothing
            # the Dirichlet projection loses a per-cent-level fraction of the mass here
            @test abs(m - 1) > 1e-3
        else
            @test isapprox(m, 1; rtol = tol)
        end
    end
end

@testset "particles outside the domain are an error" begin
    npart = 100
    pd = ParticleDistribution(1, 1, npart)
    pd.particles.v[1, :] .= range(-1, 1; length = npart)
    pd.particles.v[1, 1] = 50.0                       # outside
    pd.particles.w[1, :] .= 1 / npart

    sd = SplineDistribution(1, 1, 21, 4, (-5.0, 5.0), 0, Free())
    @test_throws ErrorException projection(pd.particles.v[1, :], pd, sd)
end

@testset "project a function and a Maxwellian" begin
    # A smooth function converges at the order of the basis. Checked by refinement rather
    # than against a fixed tolerance: a wrong projection would not converge at all.
    g = v -> exp(-v^2 / 2) / sqrt(2π)
    pts = collect(range(-3, 3; length = 13))
    # h = 0.8, 0.4, 0.2 on a distribution of unit width: coarser than this is not yet in the
    # asymptotic regime and the measured ratio means nothing.
    errs = Float64[]
    for nk in (21, 41, 81)
        s = SplineDistribution(1, 1, nk, 4, (-8.0, 8.0), 0, Free())
        project_function(g, s)
        push!(errs, maximum(abs(s.spline(v) - g(v)) for v in pts))
    end
    @test errs[1] > errs[2] > errs[3]
    @test errs[1] / errs[2] > 6          # cubic: fourth order
    @test errs[2] / errs[3] > 6

    # exact for polynomials up to the degree — the conservation predicate, and the one place a
    # fixed tolerance is right, because the answer really is exact
    sp = SplineDistribution(1, 1, 21, 4, (-3.0, 3.0), 0, Free())
    project_function(v -> 1 + 2v - 0.5v^2, sp)
    for v in (-2.0, 0.0, 1.5)
        @test isapprox(sp.spline(v), 1 + 2v - 0.5v^2; atol = 1e-9)
    end

    # 2-D Maxwellian, which used to exist only for VD == 2 and threw UndefVarError
    G = v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π)
    p2 = [(1.0, -0.5), (0.0, 0.0), (2.0, 2.0), (-1.5, 0.7)]
    e2 = Float64[]
    for nk in (11, 21, 41)
        s2 = SplineDistribution(1, 2, nk, 4, (-6.0, 6.0), 0, Free())
        project_Maxwellian(s2)
        push!(e2, maximum(abs(s2.spline(v) - G(v)) for v in p2))
    end
    @test e2[1] > e2[2] > e2[3]
    @test e2[1] / e2[2] > 6
    @test e2[2] / e2[3] > 6
    @test e2[3] < 1e-4

    # and now also for VD == 1, which had no implementation at all. The normalisation is
    # dimension-dependent: 1/sqrt(2π) in one dimension, not the 1/(2π) that was hard-coded.
    s1 = SplineDistribution(1, 1, 41, 4, (-6.0, 6.0), 0, Free())
    project_Maxwellian(s1)
    @test isapprox(s1.spline(0.0), 1 / sqrt(2π); atol = 1e-4)
    q1 = s1.quadrature
    @test isapprox(sum(quadrature_weights(q1) .* s1.spline.(quadrature_nodes(q1))), 1;
        atol = 1e-6)
end

@testset "conservative Lenard-Bernstein coefficients" begin
    npart = 5000
    Random.seed!(0xBEEF)

    pd = ParticleDistribution(1, 1, npart)
    pd.particles.v[1, :] .= randn(npart)
    pd.particles.w[1, :] .= 1 / npart

    sd = SplineDistribution(1, 1, 31, 4, (-8.0, 8.0), 0, Free())
    v = pd.particles.v[1, :]
    projection(v, pd, sd)

    A1, A2 = VM.compute_coefficients(sd, pd, v)
    @test isfinite(A1) && isfinite(A2)

    # The whole point: with these coefficients the right-hand side annihilates
    # Σ w v̇ (momentum) and Σ w v v̇ (energy).
    fs = sd.spline
    dfs = derivative(fs)
    ν = 1.0
    v̇ = @. -ν * (dfs(v) / fs(v) + A1 + A2 * v)
    w = pd.particles.w[1, :]

    scale = maximum(abs, v̇)
    @test abs(sum(w .* v̇)) < 1e-12 * scale
    @test abs(sum(w .* v .* v̇)) < 1e-11 * scale

    # ...and with NON-uniform weights too. This is what the missing w_α broke: the old code
    # annihilated Σ v̇ rather than Σ w v̇.
    pd2 = ParticleDistribution(1, 1, npart)
    pd2.particles.v[1, :] .= v
    pd2.particles.w[1, :] .= (1 .+ rand(npart)) ./ npart
    sd2 = SplineDistribution(1, 1, 31, 4, (-8.0, 8.0), 0, Free())
    projection(v, pd2, sd2)
    B1, B2 = VM.compute_coefficients(sd2, pd2, v)
    fs2 = sd2.spline
    dfs2 = derivative(fs2)
    v̇2 = @. -ν * (dfs2(v) / fs2(v) + B1 + B2 * v)
    w2 = pd2.particles.w[1, :]
    scale2 = maximum(abs, v̇2)
    @test abs(sum(w2 .* v̇2)) < 1e-12 * scale2
    @test abs(sum(w2 .* v .* v̇2)) < 1e-11 * scale2
end

@testset "entropy" begin
    npart = 5000
    Random.seed!(0x1234)
    pd = ParticleDistribution(1, 1, npart)
    pd.particles.v[1, :] .= randn(npart)
    pd.particles.w[1, :] .= 1 / npart

    sd = SplineDistribution(1, 1, 25, 4, (-6.0, 6.0), 0, Free())
    project_function(v -> exp(-v^2 / 2) / sqrt(2π), sd)

    ent = CollisionEntropy(sd)
    S = ent()
    # ∫ f log f for a unit normal is -½(1 + log 2π) ≈ -1.4189
    @test isapprox(S, -0.5 * (1 + log(2π)); atol = 1e-6)
end

@testset "Landau: K, J and L assemble" begin
    npart = 400
    Random.seed!(0xABBA)
    pd = ParticleDistribution(1, 2, npart)
    pd.particles.v .= randn(2, npart)
    pd.particles.w[1, :] .= 1 / npart

    # small basis so that L is quick
    sd = SplineDistribution(1, 2, 7, 3, (-4.0, 4.0), 0, Free())
    M = length(sd)
    projection(pd.particles.v, pd, sd)

    ent = CollisionEntropy(sd)
    landau = Landau(pd, ent)

    K1 = zeros(M, npart)
    K2 = zeros(M, npart)
    VM.compute_K!(K1, K2, pd.particles.v, sd, landau)

    # K_{k,α} = w_α ∂_c φ_k(v_α): check one particle against direct evaluation
    B = sd.basis
    α = 7
    v = pd.particles.v[:, α]
    for k in 1:M
        I = CartesianIndices(B)[k]
        d1 = VM.evaluate(B, I, (v[1], v[2]), (1, 0))
        d2 = VM.evaluate(B, I, (v[1], v[2]), (0, 1))
        @test isapprox(K1[k, α], pd.particles.w[1, α] * d1; atol = 1e-12)
        @test isapprox(K2[k, α], pd.particles.w[1, α] * d2; atol = 1e-12)
    end

    # calling twice must give the same answer: the arrays are cleared, so stale entries from
    # the previous call cannot survive
    K1b = copy(K1)
    pd.particles.v[:, 1] .+= 2.0     # move a particle to a different cell
    VM.compute_K!(K1, K2, pd.particles.v, sd, landau)
    VM.compute_K!(K1b, K2, pd.particles.v, sd, landau)
    @test K1 == K1b

    # A coarse L² projection of a peaked function undershoots into negative values, and the
    # positivity guard says so instead of returning log of a negative number. This is the
    # problem both manuscripts name and neither solves.
    project_function(v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π), sd)
    Jbad = zeros(M)
    @test_throws ErrorException VM.compute_J!(Jbad, sd, landau)

    # J = M⁻¹ ∫ φ (1 + log f_s), on a distribution that stays positive
    project_function(v -> exp(-(v[1]^2 + v[2]^2) / 2) / (2π) + 0.05, sd)
    J = zeros(M)
    VM.compute_J!(J, sd, landau)
    @test all(isfinite, J)
    # J is the projection of 1 + log f_s, so the spline it defines should match pointwise
    # J defines the spline that approximates 1 + log f_s. The basis here is deliberately tiny
    # (6 cells, quadratic) so that L is quick to assemble, so the projection error is O(0.1);
    # the point of the check is that J is that projection at all.
    JS = Spline(sd.basis, reshape(copy(J), size(sd)))
    for v in ((0.0, 0.0), (1.0, -1.0))
        f = sd.spline(v)
        @test isapprox(JS(v), 1 + log(f); atol = 0.2)
    end

    # L is symmetric
    L = zeros(M, M)
    VM.compute_L!(L, sd, landau)
    @test all(isfinite, L)
    @test maximum(abs, L - L') < 1e-10 * maximum(abs, L)
end

@testset "gradient tabulations match direct evaluation" begin
    sd = SplineDistribution(1, 2, 6, 3, (-2.0, 2.0), 0, Free())
    D1, D2 = VM.gradient_tabulations(sd)
    B = sd.basis
    q = sd.quadrature
    X = quadrature_nodes(q)
    Q1, Q2 = length(X[1]), length(X[2])
    @test size(D1) == (length(sd), Q1 * Q2)

    lin = LinearIndices(B)
    for a1 in (1, 3, Q1), a2 in (2, Q2)

        a = (a2 - 1) * Q1 + a1
        pt = (X[1][a1], X[2][a2])
        for i in (1, 2, size(B, 1)), j in (1, size(B, 2))

            k = lin[i, j]
            @test isapprox(D1[k, a], VM.evaluate(B, (i, j), pt, (1, 0)); atol = 1e-12)
            @test isapprox(D2[k, a], VM.evaluate(B, (i, j), pt, (0, 1)); atol = 1e-12)
        end
    end
end
