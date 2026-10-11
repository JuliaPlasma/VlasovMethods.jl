using GeometricBase: h5load, h5save
using GeometricIntegrators
using HDF5: h5open, h5read
using ParticleMethods: ParticleList
using PoissonSolvers
using Random
using SimpleSplines: SplineQuadrature, quadrature_nodes, quadrature_weights
using Test
using VlasovMethods
using VlasovMethods: VlasovHamiltonian, projection!

# `ParticleDistribution(xdim, vdim, npart)` always fills `Float64` zeros, so the `Float32` case
# builds its `ParticleList` from a `Float32` matrix, with the same variables.
function particle_distribution(::Type{T}, npart) where {T}
    z = zeros(T, 3, npart)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    return ParticleDistribution(1, 1, ParticleList(z; variables = vars))
end

# A member of element type `T`. The generator is seeded before `initialize!`, or the positions
# differ on every run: `NormalDistribution` draws them before its own `Random.seed!`.
function sample(::Type{T}, seed::Integer; npart = 100) where {T}
    Random.seed!(seed)
    return initialize!(particle_distribution(T, npart), NormalDistribution((T(0), T(1))))
end

function vlasov_poisson(::Type{T}; npart = 100, ncells = 16, order = 3) where {T}
    potential = Potential(PeriodicBasisSpline((T(0), T(1)), order, ncells))
    return VlasovPoisson(sample(T, 1234; npart = npart), potential)
end

# A potential on the model's basis, deposited from the positions `x` with the weights `w`, through
# the public `projection!` path.
function deposit(model, x, w)
    d = particle_distribution(eltype(x), length(x))
    d.particles.x .= x'
    d.particles.w .= w
    p = Potential(PoissonSolvers.basis(model.potential))
    projection!(p, d)
    PoissonSolvers.update!(p)
    return p
end

# `½∫(∂ₓϕ)² dx` of the potential `p`, through the Gauss quadrature of `SimpleSplines`, which is
# exact for a spline. This is the reference of the Hamiltonian's field term, and independent of
# the Galerkin identity `c'rhs` that `VlasovHamiltonian` reads it from.
function field_energy(p)
    q = SplineQuadrature(PoissonSolvers.basis(p))
    return sum(quadrature_weights(q) .* p.(quadrature_nodes(q), 1) .^ 2) / 2
end

# `½Σ w v²`, the kinetic part of the Hamiltonian, at `q`.
function kinetic(w, q)
    e = zero(eltype(q))
    for i in axes(q, 2)
        e += w[1, i] * q[2, i]^2
    end
    return e / 2
end

# `@allocated` through a fixed-arity barrier whose arguments have concrete types, so that a Julia
# 1.11 closure or splat boxing is not counted. It calls the function once before it measures.
function allocations_h(h, t, q, params)
    h(t, q, params)
    return @allocated h(t, q, params)
end

# A control for the barrier: it shows that the machinery sees an allocation, so a zero above is
# not vacuous.
function allocations(f, a)
    f(a)
    return @allocated f(a)
end

@testset "Training data as an ensemble, $T" for T in (Float64, Float32)
    tspan = (T(0), T(1))
    tstep = T(0.1)
    npart = 100
    members = [sample(T, 11; npart = npart), sample(T, 22; npart = npart)]
    model = vlasov_poisson(T; npart = npart)

    @testset "the invariant exists" begin
        ens = sodeensemble(model, members, tspan, tstep)
        inv = invariants(equation(ens))

        # The equation stays a `SODE`, so `hashamiltonian` is false and the clause is met through
        # `invariants`, as for the splitting method of the same model.
        @test !hashamiltonian(equation(ens))
        @test !(inv isa NullInvariants)
        @test haskey(inv, :h)
        @test inv.h isa VlasovHamiltonian
        @test invariants(ens) === inv
    end

    @testset "the invariant reads its argument" begin
        ens = sodeensemble(model, members, tspan, tstep)
        h = invariants(equation(ens)).h
        params = parameters(ens)[1]
        w = params.w

        # A state that is neither the model's own nor a member's, so that a Hamiltonian reading
        # `model.distribution` in place of its argument gives another value.
        q = copy(members[1].particles.z)
        q[1, :] .= T(0.1) .+ T(0.8) .* q[1, :]
        q[2, :] .+= T(0.5)
        @test q != model.distribution.particles.z

        value = h(zero(T), q, params)
        @test value ≈ kinetic(w, q) + field_energy(deposit(model, q[1, :], w)) rtol = npart *
                                                                                      eps(T)

        # One unit of velocity for every particle: the field energy of a fresh deposit of the same
        # positions does not move, so the value rises by `Σ w (v + ½)`.
        q̂ = copy(q)
        q̂[2, :] .+= one(T)
        @test h(zero(T), q̂, params) - value ≈
              sum(w[1, i] * (q[2, i] + T(1) / 2) for i in axes(q, 2)) rtol = npart * eps(T)
    end

    @testset "the members stay independent" begin
        # Two members of one draw, so the same weights, with different initial conditions: only
        # the initial condition can separate their trajectories.
        a = sample(T, 11; npart = npart)
        b = sample(T, 11; npart = npart)
        b.particles.x .= mod.(b.particles.x .+ T(0.25), 1)
        b.particles.v .*= T(0.5)
        @test b.particles.w == a.particles.w
        @test b.particles.z != a.particles.z

        ens = sodeensemble(model, [a, b], tspan, tstep)
        @test initial_conditions(ens)[1].q == a.particles.z
        @test initial_conditions(ens)[2].q == b.particles.z

        sol = integrate(ens, GeometricIntegrators.Strang())
        @test sol[1][:q] != sol[2][:q]

        # Member 2 solved alone is member 2 of the ensemble: every acceleration substep deposits
        # from the state it is given, so no member reads the field of the one solved before it.
        alone = integrate(problem(ens, 2), GeometricIntegrators.Strang())
        @test alone[:q] == sol[2][:q]
    end

    @testset "a member's weights are its own" begin
        # Two members of one draw, so the same state, whose weights differ by a factor of two; the
        # model's own weights are member 1's. A deposit reading `model.distribution.particles.w`
        # would give the two the same trajectory.
        a = sample(T, 11; npart = npart)
        b = sample(T, 11; npart = npart)
        @test b.particles.z == a.particles.z
        b.particles.w .= 2 .* a.particles.w

        model₀ = vlasov_poisson(T; npart = npart)
        model₀.distribution.particles.z .= a.particles.z
        model₀.distribution.particles.w .= a.particles.w

        ens = sodeensemble(model₀, [a, b], tspan, tstep)
        sol = integrate(ens, GeometricIntegrators.Strang())
        @test sol[1][:q] != sol[2][:q]
    end

    @testset "the weights are storable and copied" begin
        a = sample(T, 11; npart = npart)
        ens = sodeensemble(model, [a, sample(T, 22; npart = npart)], tspan, tstep)
        w = parameters(ens)[1].w

        @test w isa Matrix{T}
        @test size(w) == (1, npart)
        @test w == a.particles.w
        @test all(!iszero, w)

        # A later change to the member's distribution does not reach the problem.
        a.particles.w .= 0
        @test parameters(ens)[1].w == w
    end

    @testset "a member whose weights and state disagree" begin
        a = sample(T, 11; npart = npart)
        z = zeros(T, 3, npart)
        # the weights of one particle more than the state, as a parameter of the list, so that
        # `w` is not the view of a row of `z`
        list = ParticleList(z;
            variables = (x = 1:1, v = 2:2, z = 1:2), parameters = (w = zeros(T, 1, npart +
                                                                                   1),))
        bad = ParticleDistribution(1, 1, list)

        @test size(bad.particles.w, 2) == npart + 1
        @test size(bad.particles.z, 2) == npart
        @test_throws ArgumentError sodeensemble(model, [a, bad], tspan, tstep)
    end

    @testset "the round trip" begin
        ens = sodeensemble(model, members, tspan, tstep)
        sol = integrate(ens, GeometricIntegrators.Strang())

        h5file = joinpath(mktempdir(), "ensemble.h5")
        h5open(h5file, "w") do h5
            h5save(h5, sol)
        end
        loaded = h5open(h5file, "r") do h5
            h5load(EnsembleSolution, h5, ens)
        end

        @test datatype(loaded) === T
        for i in 1:nsamples(ens)
            @test eltype(loaded[i][:q]) === T
            @test loaded[i][:q] == sol[i][:q]
            @test parameters(ens)[i].w == members[i].particles.w
        end
        @test loaded.t == sol.t

        # A parameter tuple of the weights alone is what makes the round trip possible; the file
        # holds `(1, npart, nsamples)`.
        raw = h5read(h5file, "parameters/w")
        @test size(raw) == (1, npart, 2)
        @test raw == cat(members[1].particles.w, members[2].particles.w; dims = 3)
    end

    @testset "the Hamiltonian is type-stable and does not allocate" begin
        ens = sodeensemble(model, members, tspan, tstep)
        h = invariants(equation(ens)).h
        q = copy(members[1].particles.z)
        params = parameters(ens)[1]

        @test @inferred(h(zero(T), q, params)) isa T

        # Coverage instrumentation allocates on every counted line, so under it the assertion does
        # not run rather than being weakened.
        if Base.JLOptions().code_coverage == 0
            @test allocations_h(h, zero(T), q, params) == 0
        end
        # a control: the barrier machinery sees an allocation, so the zero above is not vacuous
        @test allocations(collect, q) > 0
    end
end
