# Regenerate the reference outputs of the collision-operator pins, one vector-field evaluation
# of each of the five operators on a StableRNGs-seeded input, plus one Picard step of the
# metriplectic operator (which exercises `f!`).
#
# The pins compare against the tree this runs on, so run it on a tree whose arithmetic is the
# intended reference:
#
#     julia --startup-file=no --project=<worktree> \
#         -e 'using TestEnv; TestEnv.activate("VlasovMethods"); include("test/helpers/generate_collision_reference.jl")'
#
# Each reference is one value per line, `%.17e`, in `test/data/`.
using VlasovMethods
using StableRNGs
using Printf
using ParticleMethods: ParticleList
using SimpleSplines: BSplineBasis, UniformMesh, Free, SplineQuadrature, nbasis,
                     TensorProductBasis, TensorProductQuadrature

const SEED = 20261010
const N = 256          # 1D1V particles; uniform in [-2, 2], so the projection stays positive
const N1 = 24          # Landau grid edge; N1^2 particles, quadratic basis for a positive projection

function write_reference(name, x)
    open(joinpath(@__DIR__, "..", "data", name), "w") do io
        for v in x
            @printf(io, "%.17e\n", v)
        end
    end
end

# --- 1D1V operators ----------------------------------------------------------------
function setup_1d1v(::Type{T}) where {T}
    rng = StableRNG(SEED)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    pdist = ParticleDistribution(1, 1, ParticleList(zeros(T, 3, N); variables = vars))
    pdist.particles.v[1, :] .= (T(4) .* rand(rng, T, N)) .- T(2)
    pdist.particles.w[1, :] .= T(1) / N
    axis = BSplineBasis(UniformMesh(17, T(-2) .. T(2)), 3, Free())
    sdist = SplineDistribution(1, 1, axis, SplineQuadrature(axis), zeros(T, nbasis(axis)))
    return pdist, sdist
end

lb_params(model) = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)

pdist, sdist = setup_1d1v(Float64)

lb = LenardBernstein(pdist, CollisionEntropy(sdist))
v = copy(pdist.particles.v[1, :])
v̇ = similar(v)
VlasovMethods.LB_rhs!(v̇, v, (ν = lb.ν, idist = lb.dist, fdist = lb.ent.dist), 0.0)
write_reference("collision_lb_reference.txt", v̇)

clb = ConservativeLenardBernstein(pdist, CollisionEntropy(sdist))
v = copy(pdist.particles.v[1, :])
v̇ = similar(v)
VlasovMethods.CLB_rhs_GI!(v̇, 0.0, v, lb_params(clb))
write_reference("collision_clb_reference.txt", v̇)

rclb = RescaledConservativeLenardBernstein(pdist, CollisionEntropy(sdist))
v = copy(pdist.particles.v[1, :])
v̇ = similar(v)
VlasovMethods.RCLB_rhs_GI!(v̇, 0.0, v, lb_params(rclb))
write_reference("collision_rclb_reference.txt", v̇)

mlb = MetriplecticLenardBernstein(pdist, CollisionEntropy(sdist))
v = copy(pdist.particles.v[1, :])
v̇ = similar(v)
VlasovMethods.collisional_vectorfield!(v̇, v, nothing, mlb)
write_reference("collision_mlb_reference.txt", v̇)

vstep = Picard_iterate_over_particles(
    zeros(N), v, v, zeros(N, 2), 1, 0.0, 1e-3, 3, 0.5, 3e-16 * sqrt(N), 1e-50, mlb)
write_reference("collision_mlb_step_reference.txt", vstep)

# --- Landau (1D2V) ----------------------------------------------------------------
function setup_1d2v(::Type{T}) where {T}
    rng = StableRNG(SEED)
    npart = N1^2
    vars = (x = 1:1, v = 2:3, z = 1:3, w = 4:4)
    pdist = ParticleDistribution(1, 2, ParticleList(zeros(T, 4, npart); variables = vars))
    pdist.particles.v .= (T(4) .* rand(rng, T, 2, npart)) .- T(2)
    pdist.particles.w[1, :] .= T(1) / npart
    axis = BSplineBasis(UniformMesh(5, T(-2) .. T(2)), 2, Free())
    B = TensorProductBasis((axis, axis))
    sdist = SplineDistribution(1, 2, B, TensorProductQuadrature(B), zeros(T, size(B)...))
    return pdist, sdist
end

lpdist, lsdist = setup_1d2v(Float64)
landau = Landau(lpdist, CollisionEntropy(lsdist))
v = copy(lpdist.particles.v)
v̇ = similar(v)
VlasovMethods.collisional_vectorfield!(v̇, v, nothing, landau)
write_reference("collision_landau_reference.txt", vec(v̇))
