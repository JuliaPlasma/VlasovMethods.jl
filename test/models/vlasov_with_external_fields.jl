using GeometricIntegrators: GeometricIntegrators
using HDF5: h5read
using LinearAlgebra: ⋅
using StaticArrays
using Test
using VlasovMethods
using VlasovMethods: electric_field, magnetic_field, s_drift!, s_lorentz!

include("../helpers/external_fields.jl")

# A constant `B = B₀ ẑ` and no electric field. The arguments are the ones every test below uses;
# `q = ±2`, `m = 3` and `B₀ = 1.7` are not powers of two, so a `q/m ↔ m/q` swap or an `abs(q)`
# changes the numbers.
function uniform_b_field(::Type{T}, B₀) where {T}
    UniformField(zero(SVector{3, T}), SVector{3, T}(0, 0, B₀))
end

# The analytic solution of `ẍ = (q/m) ẋ × B` for a constant `B = B₀ ẑ`: a circle of radius
# `m v⊥ / (|q| B₀)` about a fixed centre, and a constant velocity along `ẑ`. The sign of `q` is the
# sign of `ω`, so the direction of rotation is part of the solution the push is compared with.
function gyro_orbit(::Type{T}, x₀, v₀, q, m, B₀, t) where {T}
    ω = q * B₀ / m
    return SVector{3, T}(
        x₀[1] + (v₀[1] * sin(ω * t) + v₀[2] * (1 - cos(ω * t))) / ω,
        x₀[2] + (v₀[2] * sin(ω * t) - v₀[1] * (1 - cos(ω * t))) / ω,
        x₀[3] + v₀[3] * t)
end

# The solution of one charged particle in a field, as the integrator produces it.
function trajectory(::Type{T}, field, x₀, v₀, q, m, tspan, tstep) where {T}
    dist = one_particle(T, x₀, v₀)
    model = VlasovWithExternalFields(dist, field; charge = q, mass = m)
    method = SplittingMethod(model, tspan, tstep)
    return GeometricIntegrators.integrate(method.integrator)
end

# Every field of `T`, and every field of a field, has a concrete type. The recursion stops at a
# bitstype, an array and a function: an array's own storage is not the model's, and a function's
# body is opaque (`FieldFunction` holds a `RuntimeGeneratedFunction` whose body is an `Expr`).
function concrete_fields(T::Type)
    isconcretetype(T) || return false
    (isbitstype(T) || T <: AbstractArray || T <: Function) && return true
    fieldcount(T) == 0 && return true
    return all(concrete_fields, fieldtypes(T))
end

@testset "VlasovWithExternalFields, $T" for T in (Float64, Float32)
    m, B₀ = T(3), T(1.7)
    x₀ = SVector{3, T}(0.3, -0.2, 0.1)
    v₀ = SVector{3, T}(0.7, 0.4, 0.25)
    field = uniform_b_field(T, B₀)

    @testset "the constructor" begin
        dist = particle_distribution(T, 4)
        @test concrete_fields(typeof(VlasovWithExternalFields(dist, field)))

        # one species, with the defaults
        model = VlasovWithExternalFields(dist, field)
        @test model.charge === one(T)
        @test model.mass === one(T)

        # the keywords are converted to the particles' element type, `Int32` and `Int` included
        model = VlasovWithExternalFields(dist, field; charge = Int32(-2), mass = 3)
        @test model.charge === T(-2)
        @test model.mass === T(3)

        # a two-dimensional distribution is not a `ParticleDistribution{T,3,3}`
        @test_throws MethodError VlasovWithExternalFields(
            ParticleDistribution(1, 1, 4), field)
    end

    @testset "the mass guard" begin
        dist = particle_distribution(T, 2)
        @test_throws ArgumentError VlasovWithExternalFields(dist, field; mass = 0)
        @test_throws ArgumentError VlasovWithExternalFields(dist, field; mass = -1)
        # `NaN > 0` is false, so the guard rejects it; `Inf` passes and would give `qm = 0`
        @test_throws ArgumentError VlasovWithExternalFields(dist, field; mass = NaN)
        @test VlasovWithExternalFields(dist, field; mass = Inf).mass === T(Inf)

        # the guard reads the converted mass: `1e-50` is `0f0` in `Float32` and rejected there,
        # and positive in `Float64`
        if T === Float32
            @test_throws ArgumentError VlasovWithExternalFields(dist, field; mass = 1e-50)
        else
            @test VlasovWithExternalFields(dist, field; mass = 1e-50).mass == 1e-50
        end
    end

    @testset "the gyro-orbit" begin
        for q in (T(2), T(-2))
            # ten gyro-periods; `T_g` is positive for either sign of `q`
            T_g = 2 * T(π) * m / (abs(q) * B₀)
            err = Dict{T, T}()
            fin = Dict{T, T}()
            for (i, Δt) in enumerate((T_g / 32, T_g / 64))
                sol = trajectory(T, field, x₀, v₀, q, m, (zero(T), 10 * T_g), Δt)
                # The step count of the equation rounds up (`div(t₁ - t₀, Δt, RoundUp)`), while
                # the solution's time series is the range `t₀:Δt:t₁`, so the two disagree by one
                # whenever the span is not an exact multiple. Both errors are therefore taken at
                # the times the solution carries, never at an assumed count.
                N = lastindex(sol.q)
                err[i] = zero(T)
                for n in 0:N
                    x = sol.q[n]
                    d = SVector{3, T}(x[1, 1], x[2, 1], x[3, 1]) -
                        gyro_orbit(T, x₀, v₀, q, m, B₀, sol.t[n])
                    err[i] = max(err[i], sqrt(sum(abs2, d)))
                end
                x = sol.q[N]
                fin[i] = sqrt(sum(abs2,
                    SVector{3, T}(x[1, 1], x[2, 1], x[3, 1]) -
                    gyro_orbit(T, x₀, v₀, q, m, B₀, sol.t[N])))
            end
            # The push is second order, so halving the step divides the error by four. The window
            # is ±10% about 4. A reflected rotation is off by up to the diameter whatever the step,
            # so its ratio is near 1; a whole number of periods brings the wrong direction back to
            # its start in `x`–`y`, which is why the ratio uses the maximum over the trajectory.
            @test 3.6 ≤ err[1] / err[2] ≤ 4.4
            @test 3.6 ≤ fin[1] / fin[2] ≤ 4.4
        end
    end

    @testset "the energy" begin
        for q in (T(2), T(-2))
            T_g = 2 * T(π) * m / (abs(q) * B₀)
            test_energy(trajectory(T, field, x₀, v₀, q, m, (zero(T), 10 * T_g), T_g / 32), v₀)
        end
    end

    @testset "a field with E ∥ B" begin
        # With `E = E₀ ẑ` and `B = B₀ ẑ`, the rotation leaves `v_z` alone, so the motion along
        # `ẑ` is `z₀ + v∥ t + ½(q/m)E₀t²` exactly for this splitting: the two half drifts carry
        # the kicks' velocity change in full.
        q, E₀, tstep = T(2), T(0.3), T(4) / T(37)
        local field = UniformField(SVector{3, T}(0, 0, E₀), SVector{3, T}(0, 0, B₀))
        sol = trajectory(T, field, x₀, v₀, q, m, (zero(T), T(4)), tstep)
        N = lastindex(sol.q)
        t = sol.t[N]
        z = sol.q[N]
        @test z[3, 1] ≈ x₀[3] + v₀[3] * t + q * E₀ * t^2 / (2 * m) rtol = N * eps(T)
        # the kick is the only thing that moves `v_z`, and it does so by `(q/m)E₀` per unit time
        @test z[6, 1] ≈ v₀[3] + q * E₀ * t / m rtol = N * eps(T)
    end

    @testset "free streaming" begin
        # A zero field leaves `v` bitwise unchanged: `τ = 0` makes `s = 0` and `qh·E = 0`
        # exactly. Two `Int` coordinates make a `UniformField{Int}`, which the push converts to
        # the particles' element type.
        local field = UniformField(SVector(0, 0, 0), SVector(0, 0, 0))
        @test field isa UniformField{Int}
        @test all(iszero, field.E) && all(iszero, field.B)
        q, tstep = T(-2), T(0.1)
        sol = trajectory(T, field, x₀, v₀, q, m, (zero(T), T(1)), tstep)
        N = lastindex(sol.q)
        t = sol.t[N]
        z = sol.q[N]
        @test all(z[4:6, 1] .=== v₀)
        @test z[1:3, 1] ≈ x₀ + v₀ * t rtol = N * eps(T)
    end

    @testset "zero and one particle" begin
        for npart in (0, 1)
            dist = particle_distribution(T, npart)
            npart == 1 && (dist.particles.v .= v₀)
            model = VlasovWithExternalFields(dist, field; charge = T(-2), mass = m)
            method = SplittingMethod(model, (zero(T), T(1)), T(0.1))
            params = method.equation.parameters
            z̄ = copy(model.distribution.particles.z)
            z = similar(z̄)
            s_drift!(z, T(0.1), z̄, zero(T), params)
            s_lorentz!(z, T(0.1), z̄, zero(T), params)
            @test size(z) == (6, npart)
            @test all(isfinite, z)
            if npart == 1
                @test (z[4, 1]^2 + z[5, 1]^2 + z[6, 1]^2) / 2 ≈ v₀ ⋅ v₀ / 2 rtol = √eps(T)
            end

            # The empty list is a vacuum run, not an error: `GeometricBase` 0.15.1 keeps the
            # vector field of an empty state variable, so `integrate` returns an empty state at
            # the time points the solution carries.
            sol = GeometricIntegrators.integrate(method.integrator)
            @test size(sol.q[end]) == (6, npart)
            @test all(isfinite, sol.q[end])
            @test length(sol.t) == 11

            # `run!` serves the empty list too: its chunk has no zero extent, and the dataset
            # holds one slice per time point, the last of them the model's own state.
            h5file = tempname() * ".h5"
            try
                @test run!(method, h5file) === model.distribution
                h5z = h5read(h5file, "z")
                @test size(h5z) == (6, npart, 11)
                @test h5z[:, :, end] == model.distribution.particles.z
            finally
                rm(h5file; force = true)
            end
        end
    end

    @testset "the two vector fields" begin
        # `E` and `B` both nonzero, and no value a power of two
        local field = UniformField(SVector{3, T}(0.2, -0.1, 0.3), SVector{3, T}(0, 0, B₀))
        test_vector_fields(T, field, SVector{3, T}(0.3, -0.2, 0.1), v₀)
    end

    @testset "the push does not allocate" begin
        test_push_does_not_allocate(T, field)
    end

    @testset "the field of another precision" begin
        # `electric_field` returns the element type of the coordinates, whatever the field's own
        # precision: the model converts, so a `Float32` particle sees a `Float32` field.
        x32 = SVector{3, Float32}(0.1, 0.2, 0.3)
        f64 = UniformField(SVector(0.0, 0.0, 0.3), SVector(0.0, 0.0, 1.7))
        @test electric_field(f64, 0.0f0, x32) isa SVector{3, Float32}
        @test magnetic_field(f64, 0.0f0, x32) isa SVector{3, Float32}
        @test electric_field(f64, 0.0f0, x32) == SVector{3, Float32}(0, 0, 0.3)
        @test magnetic_field(f64, 0.0f0, x32) == SVector{3, Float32}(0, 0, 1.7)

        # and the model with a `Float64` field and `Float32` particles runs
        dist = particle_distribution(Float32, 1)
        model = VlasovWithExternalFields(dist, f64; charge = 2, mass = 3)
        method = SplittingMethod(model, (0.0f0, 1.0f0), 0.1f0)
        @test eltype(GeometricIntegrators.integrate(method.integrator).q[end]) === Float32
    end

    @testset "a field with no method" begin
        # A length other than three is the `DimensionMismatch` of `SVector{3}`, and a field of
        # no method raises its `MethodError` at the first kick — the constructor does not check
        # the field, because the model puts no bound on it.
        @test_throws DimensionMismatch UniformField(zeros(2), zeros(3))
        @test_throws DimensionMismatch UniformField(zeros(3), zeros(4))
        model = VlasovWithExternalFields(
            particle_distribution(T, 1), nothing; charge = T(2), mass = m)
        method = SplittingMethod(model, (zero(T), T(1)), T(0.1))
        @test_throws MethodError GeometricIntegrators.integrate(method.integrator)
    end
end
