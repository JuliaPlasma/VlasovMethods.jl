using ParticleMethods: ParticleList
using SimpleSplines: TensorProductQuadrature
using Test
using VlasovMethods

# No other test calls `Picard_iterate_Landau_nls!`. This call fails on any undefined global in
# its body.
@testset "the Landau Picard solver runs" begin
    n1 = 16
    npart = n1^2
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pdist = ParticleDistribution(1, 2, ParticleList(zeros(4, npart); variables = vars))
    # Cell centres inside the velocity domain, so the projected density stays positive.
    g = [-2.0 + (i - 0.5) * 4.0 / n1 for i in 1:n1]
    pdist.particles.v[1, :] .= repeat(g; outer = n1)
    pdist.particles.v[2, :] .= repeat(g; inner = n1)
    pdist.particles.w[1, :] .= 1 / npart

    axis = BSplineBasis(UniformMesh(8, -2.0 .. 2.0), 3, Free())
    B = TensorProductBasis((axis, axis))
    sdist = SplineDistribution(1, 2, B, TensorProductQuadrature(B), zeros(size(B)...))
    landau = Landau(pdist, CollisionEntropy(sdist))

    v_prev = copy(pdist.particles.v)
    rhs_prev = zeros(2, npart, 2)

    sol = redirect_stdout(devnull) do
        VlasovMethods.Picard_iterate_Landau_nls!(
            landau, 1e-10, 1e-10, 0.5, 1e-6, 1, 0.0, v_prev, copy(v_prev), rhs_prev, 3, 1)
    end

    @test sol.v isa Vector
    @test length(sol.v) == 6
    @test sol.v̇ isa Vector
end
