using ParticleMethods: ParticleList
using SimpleSplines: TensorProductQuadrature
using Test
using VlasovMethods

# No other test calls `Picard_iterate_Landau_nls!`. It is called once on each side of its
# `ti ≥ 4` branch, so an undefined global in either branch throws here.
function landau_solver_setup()
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
    return Landau(pdist, CollisionEntropy(sdist)), copy(pdist.particles.v)
end

function landau_step(landau, v_prev, ti)
    rhs_prev = zeros(2, size(v_prev, 2), 2)
    redirect_stdout(devnull) do
        VlasovMethods.Picard_iterate_Landau_nls!(
            landau, 1e-10, 1e-10, 0.5, 1e-6, ti, 0.0, v_prev, copy(v_prev), rhs_prev, 3, 1)
    end
end

@testset "the Landau Picard solver runs" begin
    landau1, v_prev1 = landau_solver_setup()
    sol1 = landau_step(landau1, v_prev1, 1)
    @test sol1.v isa Vector
    @test length(sol1.v) == 6
    @test sol1.v̇ isa Vector

    landau4, v_prev4 = landau_solver_setup()
    sol4 = landau_step(landau4, v_prev4, 4)
    @test sol4.v isa Vector
    @test length(sol4.v) == 6
    @test sol4.v̇ isa Vector

    # `sol.v[1]` is the extrapolated initial guess. The `ti ≥ 4` call starts from a Hermite guess
    # and the `ti < 4` call from a midpoint one, so the two first entries differ.
    @test sol4.v[1] != sol1.v[1]
end
