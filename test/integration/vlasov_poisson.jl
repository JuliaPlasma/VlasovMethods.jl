using GeometricIntegrators: GeometricIntegrators
using GeometricEquations: ntime
using HDF5: h5read
using ParticleMethods: ParticleList
using PoissonSolvers
using Random
using Test
using VlasovMethods
using VlasovMethods: projection!, update_potential!, v_acceleration!, s_acceleration!

# `ParticleDistribution(xdim, vdim, npart)` always fills `Float64` zeros, so the `Float32` case
# builds its `ParticleList` from a `Float32` matrix, with the same variables.
function particle_distribution(::Type{T}, npart) where {T}
    z = zeros(T, 3, npart)
    vars = (x = 1:1, v = 2:2, z = 1:2, w = 3:3)
    ParticleDistribution(1, 1, ParticleList(z; variables = vars))
end

# `NormalDistribution` draws `x₀` before its own `Random.seed!`, so the generator is seeded here,
# before `initialize!`, or the positions differ on every run.
function vlasov_poisson(::Type{T}; npart = 1000, ncells = 16, order = 3) where {T}
    Random.seed!(1234)
    dist = initialize!(particle_distribution(T, npart), NormalDistribution((T(0), T(1))))
    potential = Potential(PeriodicBasisSpline((T(0), T(1)), order, ncells))
    VlasovPoisson(dist, potential)
end

# A potential on the same basis, deposited from `x` with the weights of `model`, through the
# public `projection!` path.
function deposit(model, x)
    p = Potential(PoissonSolvers.basis(model.potential))
    d = particle_distribution(eltype(x), length(x))
    d.particles.x .= x'
    d.particles.w .= model.distribution.particles.w
    projection!(p, d)
    PoissonSolvers.update!(p)
    return p
end

# `@allocated` through a fixed-arity barrier whose arguments have concrete types, so that a Julia
# 1.11 closure or splat boxing is not counted (`evidence.md`, Allocation assertions). Each calls
# the function once before it measures.
function allocations_s(f, z, t, z̄, t̄, params)
    f(z, t, z̄, t̄, params)
    return @allocated f(z, t, z̄, t̄, params)
end
function allocations_v(f, z, t, z̄, params)
    f(z, t, z̄, params)
    return @allocated f(z, t, z̄, params)
end

@testset "Vlasov–Poisson, $T" for T in (Float64, Float32)
    tspan = (T(0), T(1))
    tstep = T(0.1)

    @testset "integration" begin
        model = vlasov_poisson(T)
        method = SplittingMethod(model, tspan, tstep)
        params = method.equation.parameters

        # ten steps of `tstep` over `tspan`
        @test ntime(method.equation) == 10
        sol = GeometricIntegrators.integrate(method.integrator)
        @test all(isfinite, sol.q[end])
        @test eltype(sol.q[end]) == T
        @test sol.q[end] != sol.q[0]

        z̄ = copy(sol.q[end])
        z = similar(z̄)
        @test @inferred(v_acceleration!(z, tstep, z̄, params)) === nothing
        @test @inferred(s_acceleration!(z, tstep, z̄, zero(T), params)) === nothing
    end

    @testset "the field follows the particles" begin
        model = vlasov_poisson(T)
        method = SplittingMethod(model, tspan, tstep)
        params = method.equation.parameters
        x₀ = vec(copy(model.distribution.particles.x))

        # positions that differ from the model's own
        z̄ = copy(model.distribution.particles.z)
        z̄[1, :] .= mod.(z̄[1, :] .+ T(0.25) .* sinpi.(2 .* z̄[1, :]), 1)
        @test z̄[1, :] != x₀
        z = similar(z̄)
        s_acceleration!(z, tstep, z̄, zero(T), params)

        reference = deposit(model, z̄[1, :])
        @test PoissonSolvers.coefficients(model.potential) ==
              PoissonSolvers.coefficients(reference)

        # The field is evaluated at `z̄`. The reference evaluates the same spline through
        # `Spline`'s own path, which sums in another order, so the two agree to rounding.
        @test z[1, :] == z̄[1, :]
        @test z[2, :] ≈ z̄[2, :] .- tstep .* reference.(z̄[1, :], 1) rtol=√eps(T)

        # the coefficients at step 0 and at step 10 differ
        update_potential!(model)
        c₀ = copy(PoissonSolvers.coefficients(model.potential))
        method = SplittingMethod(vlasov_poisson(T), tspan, tstep)
        GeometricIntegrators.integrate(method.integrator)
        c₁₀ = PoissonSolvers.coefficients(method.model.potential)
        @test c₀ == PoissonSolvers.coefficients(deposit(model, x₀))
        @test c₁₀ != c₀
    end

    @testset "periodic wrap-around" begin
        model = vlasov_poisson(T; ncells = 16)
        x = vec(copy(model.distribution.particles.x))
        w = model.distribution.particles.w
        ρ = copy(PoissonSolvers.rhs(deposit(model, x)))

        # The argument reduction moves `x` by about one ulp, and the slope of a basis function is
        # below `2/h`, so each deposit moves by at most `2·ncells·eps(T)` times its weight.
        tolerance = 2 * 16 * eps(T) * sum(abs, w)
        for shift in (T(1), -T(1), T(3))
            ρ̃ = PoissonSolvers.rhs(deposit(model, x .+ shift))
            @test maximum(abs, ρ̃ .- ρ) ≤ tolerance
        end
    end

    @testset "a Dirichlet basis does not wrap" begin
        # Only a periodic basis wraps. A Dirichlet basis takes each position unchanged, so a
        # particle outside the domain deposits nothing, as it did before the wrap was added.
        b = DirichletBasisSpline((T(0), T(1)), 3, 16)
        x = T[0.3, 0.71, 1.02, -0.01, 1.3]
        w = fill(T(0.2), length(x))
        d = particle_distribution(T, length(x))
        d.particles.x .= x'
        d.particles.w .= w'
        p = Potential(b)
        projection!(p, d)

        # the deposit loop of the tree before this part, which evaluates at `x` as it is
        reference = zeros(T, length(PoissonSolvers.rhs(p)))
        vals = zeros(T, local_width(b))
        for (xᵢ, wᵢ) in zip(x, w)
            j₀ = evaluate_all!(vals, b, xᵢ)
            for (t, value) in pairs(vals)
                iszero(value) && continue
                reference[basis_index(b, j₀ + t - 1)] += wᵢ * value
            end
        end
        @test PoissonSolvers.rhs(p) == reference
    end

    @testset "the HDF5 output" begin
        npart = 100
        model = vlasov_poisson(T; npart)
        z₀ = copy(model.distribution.particles.z)
        method = SplittingMethod(model, tspan, tstep)
        h5file = joinpath(mktempdir(), "vlasov_poisson.hdf5")
        run!(method, h5file)
        z = h5read(h5file, "z")
        @test size(z) == (2, npart, ntime(method.equation) + 1)
        @test z[:, :, 1] == z₀
    end

    @testset "the right-hand side does not allocate" begin
        model = vlasov_poisson(T)
        method = SplittingMethod(model, tspan, tstep)
        params = method.equation.parameters
        z̄ = copy(model.distribution.particles.z)
        z = similar(z̄)
        t, t̄ = tstep, zero(T)

        # Coverage instrumentation allocates on every counted line, so under it the assertion does
        # not run rather than being weakened.
        if Base.JLOptions().code_coverage == 0
            @test allocations_s(s_acceleration!, z, t, z̄, t̄, params) == 0
            @test allocations_v(v_acceleration!, z, t, z̄, params) == 0
        end
    end
end
