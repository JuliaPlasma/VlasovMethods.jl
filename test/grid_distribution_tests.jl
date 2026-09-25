using VlasovMethods
using Test

const VM = VlasovMethods

# ∫ π(v) f_s(v) dv by the Gauß-Legendre rule the spline basis carries, which is exact for a
# cubic spline times a polynomial of degree ≤ 2.
function spline_moment(π, s::SplineDistribution)
    v = quadrature_nodes(s.quadrature)
    w = quadrature_weights(s.quadrature)
    sum(w[r] * π(v[r]) * s.spline(v[r]) for r in eachindex(v))
end

# The three moments of the grid and of the spline for f(x, v) = g(x) h(v), at each x-node.
function moments_on_both(g, h, nx, nv, vdomain; nknots = 81)
    grid = GridDistribution(nx, nv, (0.0, 2π), vdomain)
    grid.values .= g.(grid.x) .* h.(grid.v')
    s = SplineDistribution(1, 1, nknots, 4, vdomain)
    project_function(h, s)
    spline = map(π -> g.(grid.x) .* spline_moment(π, s), (v -> 1, v -> v, v -> v^2 / 2))
    return grid, velocity_moments(grid), spline
end

@testset "GridDistribution construction" begin
    g = GridDistribution(8, 5, (0.0, 2π), (-4.0, 4.0))
    @test size(g) == (8, 5)
    @test length(g) == 40
    # the x-grid is periodic: 2π is the node 0 again and is not stored
    @test first(g.x) == 0.0
    @test step(g.x) ≈ 2π / 8
    @test last(g.x) ≈ 2π - 2π / 8
    # the v-grid is bounded: both ends are nodes
    @test first(g.v) == -4.0
    @test last(g.v) == 4.0
    @test step(g.v) == 2.0

    @test_throws DimensionMismatch GridDistribution(g.x, g.v, zeros(8, 4))
    @test_throws ArgumentError GridDistribution(zeros(8, 1), (0.0, 1.0), (0.0, 1.0))
end

@testset "The three distributions share one interface" begin
    p = ParticleDistribution(1, 1, 10)
    s = SplineDistribution(1, 1, 41, 4, (-8.0, 8.0))
    g = GridDistribution(16, 65, (0.0, 2π), (-8.0, 8.0))

    for d in (p, s, g)
        @test d isa VM.DistributionFunction{Float64, 1, 1}
        @test eltype(d) == Float64
        @test VM.xdim(d) == 1
        @test VM.vdim(d) == 1
    end
    @test size(p) == size(p.particles)
    @test size(s) == (43,)
    @test size(g) == (16, 65)

    # the spline and the grid evaluate through the same call, `d(x, v)`
    h(v) = exp(-v^2 / 2)
    project_function(h, s)
    g.values .= h.(g.v')
    for (x, v) in ((0.3, 0.7), (5.0, -1.2), (1.0, 2.9))
        @test s(x, v) ≈ h(v) atol = 1e-3
        @test g(x, v) ≈ h(v) atol = 2e-2
    end
end

@testset "GridDistribution evaluation" begin
    g = GridDistribution(8, 5, (0.0, 2π), (-2.0, 2.0))
    g.values .= [cos(x) + 3v for x in g.x, v in g.v]

    # the nodes are reproduced, and the right end of the v-range is a node
    @test g(g.x[3], g.v[2]) ≈ g.values[3, 2]
    @test g(g.x[5], 2.0) ≈ g.values[5, 5]
    # between nodes the interpolant is linear in v
    @test g(g.x[3], 0.25) ≈ cos(g.x[3]) + 0.75
    # periodic in x, including the node at the right end of the x-range
    @test g(2π, 0.3) ≈ g(0.0, 0.3)
    @test g(1.1 + 2π, 0.3) ≈ g(1.1, 0.3)
    @test g(1.1 - 4π, 0.3) ≈ g(1.1, 0.3)
    # zero outside the v-range
    @test g(1.0, 2.0 + 1e-12) == 0
    @test g(1.0, -2.5) == 0
    # on [-1, 0.3] with 10 nodes, (v_b - v_a) / h_v rounds above 9, and v_b is still a node
    @test GridDistribution(ones(4, 10), (0.0, 1.0), (-1.0, 0.3))(0.5, 0.3) ≈ 1

    # bilinear interpolation is second order: halving both steps divides the error by four
    f(x, v) = (1 + 0.5 * cos(x)) * exp(-v^2 / 2)
    points = [(x, v)
              for x in range(0.1, 6.1; length = 7), v in range(-2.9, 2.9; length = 7)]
    errors = map((16, 32, 64)) do n
        grid = GridDistribution(n, 2n + 1, (0.0, 2π), (-4.0, 4.0))
        grid.values .= f.(grid.x, grid.v')
        maximum(abs(grid(x, v) - f(x, v)) for (x, v) in points)
    end
    @test all(3.5 .< errors[1:(end - 1)] ./ errors[2:end] .< 4.5)
end

@testset "Grid moments against the spline" begin
    g(x) = 1 + 0.3 * cos(x)
    u = 0.5
    M(v) = exp(-(v - u)^2 / 2) / sqrt(2π)

    # A Maxwellian that vanishes at the ends of the v-range. The L² projection onto the
    # clamped cubic basis keeps ∫ vᵏ h dv for k ≤ 2 exactly, and the rectangle rule is exact up
    # to the end values, so the two agree to round-off, and both give the analytic moments.
    grid, mg, ms = moments_on_both(g, M, 16, 161, (-10.0, 10.0))
    @test mg.density ≈ ms[1] rtol = 1e-12
    @test mg.momentum ≈ ms[2] rtol = 1e-12
    @test mg.energy ≈ ms[3] rtol = 1e-12
    @test mg.density ≈ g.(grid.x) rtol = 1e-12
    @test mg.momentum ≈ u .* g.(grid.x) rtol = 1e-12
    @test mg.energy ≈ (1 + u^2) / 2 .* g.(grid.x) rtol = 1e-12

    # On [-1, 2] the Maxwellian does not vanish at the ends, and the rectangle rule is first
    # order: its error is hᵥ/2 (π(vₐ) h(vₐ) + π(v_b) h(v_b)) plus a second-order remainder.
    va, vb = -1.0, 2.0
    for (k, π) in enumerate((v -> 1, v -> v, v -> v^2 / 2))
        err = Float64[]
        rem = Float64[]
        for nv in (31, 61, 121)
            grid, mg, ms = moments_on_both(g, M, 16, nv, (va, vb))
            e = maximum(abs, (mg.density, mg.momentum, mg.energy)[k] .- ms[k])
            hv = step(grid.v)
            boundary = hv / 2 * (π(va) * M(va) + π(vb) * M(vb))
            r = maximum(abs, (mg.density, mg.momentum, mg.energy)[k] .- ms[k] .-
                             boundary .* g.(grid.x))
            push!(err, e)
            push!(rem, r)
        end
        @test all(1.9 .< err[1:(end - 1)] ./ err[2:end] .< 2.1)
        @test all(3.8 .< rem[1:(end - 1)] ./ rem[2:end] .< 4.2)
    end
end
