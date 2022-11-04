using OffsetArrays
using PoissonSolvers
using Test
using VlasovMethods
using VlasovMethods.BumpOnTail


nd = 1
np = 100
nt = 10
Δt = 0.1
χ  = 1.0

nh = 16
p = 3
L = 2π

poisson = PoissonSolverPBSplines(p, nh, L)

params = (κ = 0.3, ε = 0.03, a = 0.1, v₀= 4.5, σ = 0.5, χ = 1.0)

particles = BumpOnTail.draw_accept_reject(np, params)


p = PoissonField(poisson)
ep = p(particles.x, particles.w, 0.0)


ϕ = hcat(p.poisson.ϕ, rand(nh,nt))

e1 = ExternalField(poisson, OffsetArray(ϕ, 1:nh, 0:nt), Δt)
e2 = ExternalField(poisson, ϕ, Δt)

ef1 = e1(particles.x, particles.w, 0.0)
ef2 = e2(particles.x, particles.w, 0.0)

@test ep == ef1 == ef2


s1 = ScaledField(p, χ)
es1 = s1(particles.x, particles.w, 0.0)

s2 = ScaledField(e2, χ)
es2 = s2(particles.x, particles.w, 0.0)

@test ep == es1 == es2
