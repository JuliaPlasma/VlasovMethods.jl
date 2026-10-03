using GeometricIntegrators: GeometricIntegrators
using GeometricEquations: ntime
using HDF5: h5read
using ParticleMethods: ParticleList
using PoissonSolvers
using Random
using Test
using VlasovMethods
using VlasovMethods: projection!, update_potential!, v_acceleration!, s_acceleration!,
                     _electric_field

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
# 1.11 closure or splat boxing is not counted. Each calls the function once before it measures.
function allocations_s(f, z, t, z̄, t̄, params)
    f(z, t, z̄, t̄, params)
    return @allocated f(z, t, z̄, t̄, params)
end
function allocations_v(f, z, t, z̄, params)
    f(z, t, z̄, params)
    return @allocated f(z, t, z̄, params)
end

# A generic barrier for the control below: it shows that the machinery sees an allocation, so a
# zero from the fixed-arity barriers above is not vacuous.
function allocations(f, a)
    f(a)
    return @allocated f(a)
end

@testset "Vlasov–Poisson, $T" for T in (Float64, Float32)
    tspan = (T(0), T(1))
    tstep = T(0.1)

    @testset "integration" begin
        model = vlasov_poisson(T)
        method = SplittingMethod(model, tspan, tstep)
        params = method.equation.parameters

        # the held buffer carries the element type of the potential's right-hand side
        @test eltype(model.work) === T

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

    @testset "the acceleration field" begin
        model = vlasov_poisson(T)
        method = SplittingMethod(model, tspan, tstep)
        params = method.equation.parameters
        z = copy(model.distribution.particles.z)
        ż = similar(z)
        v_acceleration!(ż, tstep, z, params)

        @test all(iszero, ż[1, :])

        # the field is `-∂ₓϕ`: the derivative of the deposit from the same positions, negated
        reference = deposit(model, z[1, :])
        @test ż[2, :] ≈ -reference.(z[1, :], 1) rtol=√eps(T)
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

    # A pin of the SimpleSplines contract: `evaluate_all!` reduces a periodic argument, so the
    # state may keep an unwrapped position. VlasovMethods adds no reduction of its own.
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

    @testset "a periodic domain that begins away from zero" begin
        # The reduction is modulo the domain length *and* its origin, so a reduction that
        # dropped the origin would shift this deposit. A pin of the SimpleSplines contract.
        b = PeriodicBasisSpline((T(-1), T(1)), 3, 16)
        x = T[-0.7, -0.2, 0.1, 0.55, 0.9]
        w = T[0.1, 0.2, 0.7, 0.05, 0.35]
        p = Potential(b)
        d = particle_distribution(T, length(x))
        d.particles.x .= x'
        d.particles.w .= w'
        projection!(p, d)
        ρ = copy(PoissonSolvers.rhs(p))

        p2 = Potential(b)
        d2 = particle_distribution(T, length(x))
        d2.particles.x .= (x .+ T(2))'   # one domain length of 2
        d2.particles.w .= w'
        projection!(p2, d2)

        tolerance = 2 * 16 * eps(T) * sum(abs, w)
        @test maximum(abs, PoissonSolvers.rhs(p2) .- ρ) ≤ tolerance
    end

    @testset "the field evaluation reduces an unwrapped position" begin
        # As for the deposit, the reduction belongs to `evaluate_all!`. A position one domain
        # length above the domain gives the field of its reduction.
        model = vlasov_poisson(T)
        x = vec(copy(model.distribution.particles.x))
        p = deposit(model, x)
        work = zeros(T, local_width(PoissonSolvers.basis(p)))
        @test _electric_field(p, work, T(1.3)) ≈ _electric_field(p, work, T(0.3)) rtol=√eps(T)
    end

    # A pin of the SimpleSplines contract: only a periodic basis reduces its argument, so a
    # Dirichlet basis takes each position unchanged.
    @testset "a Dirichlet basis does not wrap" begin
        # A particle outside the domain deposits nothing.
        b = DirichletBasisSpline((T(0), T(1)), 3, 16)
        x = T[0.3, 0.71, 1.02, -0.01, 1.3]
        w = fill(T(0.2), length(x))
        d = particle_distribution(T, length(x))
        d.particles.x .= x'
        d.particles.w .= w'
        p = Potential(b)
        projection!(p, d)

        # a reference deposit loop that evaluates at `x` as it is
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

    @testset "the field evaluation on a Dirichlet basis" begin
        # A clamped basis pads the derivative buffer with zeros whose indices run past
        # `nbasis`; the field evaluation skips them rather than indexing past the coefficients.
        b = DirichletBasisSpline((T(0), T(1)), 3, 16)
        p = Potential(b)
        d = particle_distribution(T, 5)
        d.particles.x .= T[0.2, 0.4, 0.6, 0.8, 0.9]'
        d.particles.w .= T[0.2, 0.2, 0.2, 0.2, 0.2]'
        projection!(p, d)
        PoissonSolvers.update!(p)

        # 0.999 lies in the last cell, where `evaluate_all!` leaves the tail of the buffer zero
        work = zeros(T, local_width(b))
        @test isfinite(_electric_field(p, work, T(0.999)))
    end

    @testset "the deposit reads each particle's weight" begin
        # Distinct weights: a deposit that read only the first weight would give another result.
        b = PeriodicBasisSpline((T(0), T(1)), 3, 16)
        x = T[0.05, 0.3, 0.62, 0.91, 0.44]
        w = T[0.1, 0.2, 0.7, 0.05, 0.35]
        d = particle_distribution(T, length(x))
        d.particles.x .= x'
        d.particles.w .= w'
        p = Potential(b)
        projection!(p, d)

        # a reference deposit loop, each particle with its own weight
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
        # the run is written out step by step, not shifted by one from the initial slice
        @test z[:, :, 2] != z[:, :, 1]
        # the last slice is the final integration state, written back into the model's particles
        @test z[:, :, end] == model.distribution.particles.z
        # each slice is the state of the same step of `integrate` on the same problem
        sol = GeometricIntegrators.integrate(
            SplittingMethod(vlasov_poisson(T; npart), tspan, tstep).integrator)
        same = [z[:, :, n + 1] == sol.q[n] for n in 0:ntime(method.equation)]
        @test all(same)
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
        # a control: the barrier machinery sees an allocation, so the zeros above are not vacuous
        @test allocations(collect, z) > 0
    end
end
