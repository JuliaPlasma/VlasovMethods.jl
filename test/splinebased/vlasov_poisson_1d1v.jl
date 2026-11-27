using Test
import HDF5
using HDF5: h5open, h5read
using Random: Xoshiro
using SimpleSplines: BSplineBasis, PeriodicBSplineBasis, UniformMesh, Dirichlet, (..)
using SimpleSplines: GeneralMesh
using PoissonSolvers: FFTWBasis, Potential
using VlasovMethods
using VlasovMethods: _galerkin_form, _averaged_form, _spline_1d1v_rhs!, _spline_energy,
                     _spline_mass

include(joinpath(@__DIR__, "..", "helpers", "spline_galerkin.jl"))

const GS = VlasovMethods.GeometricSolvers

median(x) = sort(x)[(length(x) + 1) ÷ 2]

function spline_model(T, f; nv = 100)
    xb, vb = spline_bases(T; nv)
    dist = TensorSplineDistribution(xb, vb, f)
    return VlasovPoisson(dist, Potential(xb)), spline_tables(xb, vb)
end

field(model, A) = _spline_1d1v_rhs!(similar(A), model, A)

# `run!` into a new file, and the states and the times it writes
function run_states(method)
    file = tempname(; suffix = ".hdf5")
    dist = run!(method, file)
    return dist, h5read(file, "z"), h5read(file, "t")
end

# the number of HDF5 files that the process holds open
open_hdf5_files() = HDF5.API.h5f_get_obj_count(HDF5.API.H5F_OBJ_ALL, HDF5.API.H5F_OBJ_FILE)

# the test's own implicit midpoint step of the plain field, by Newton–GMRES; the largest relative
# energy change over `nsteps` steps
function plain_midpoint_drift(t, A, nsteps, τ, tol)
    T = eltype(A)
    nx, nv = size(A)
    H₀ = spline_H(t, A)
    drift = zero(T)
    for _ in 1:nsteps
        A₀ = A
        function F!(r, y, p)
            f̄ = (A₀ .+ reshape(y, nx, nv)) ./ 2
            r .= y .- vec(A₀) .- τ .* vec(plain_bracket(t, spline_phi(t, f̄), true, f̄))
            return r
        end
        function jvp!(Jv, y, w, p)
            f̄ = (A₀ .+ reshape(y, nx, nv)) ./ 2
            W = reshape(w, nx, nv)
            J = plain_bracket(t, spline_phi(t, W), false, f̄) .+
                plain_bracket(t, spline_phi(t, f̄), true, W)
            Jv .= w .- (τ / 2) .* vec(J)
            return Jv
        end
        method = GS.Newton(;
            linsolve = GS.GMRES(; maxiter = 100, forcing = GS.FixedForcing(1e-3)),
            jacobian = GS.JVPFunction(jvp!))
        options = GS.Options(T; f_abstol = tol * norm(A₀), f_reltol = 0, x_reltol = 0,
            max_iterations = 10)
        y, status = GS.solve(GS.NonlinearProblem(F!), copy(vec(A₀)), method; options)
        @test status.code === GS.SUCCESS
        A = reshape(y, nx, nv)
        drift = max(drift, abs(spline_H(t, A) - H₀) / abs(H₀))
    end
    return drift
end

@testset "construction ($T)" for T in (Float64, Float32)
    xb, vb = spline_bases(T)
    dist = TensorSplineDistribution(xb, vb, landau_state(T))

    # the potential is a spline potential on the x-basis of the distribution
    @test_throws r"^ArgumentError: the potential of the spline Vlasov–Poisson model solves -φ'' = ρ by a PoissonSolverSpline" VlasovPoisson(
        dist, Potential(FFTWBasis((zero(T), T(4π)), 16)))
    other_basis = r"^ArgumentError: the potential of the spline Vlasov–Poisson model is on the x-basis"
    for b in (PeriodicBSplineBasis(UniformMesh(16, zero(T) .. T(2π)), 3),
        PeriodicBSplineBasis(UniformMesh(16, zero(T) .. T(4π)), 2),
        PeriodicBSplineBasis(UniformMesh(32, zero(T) .. T(4π)), 3))
        @test_throws other_basis VlasovPoisson(dist, Potential(b))
    end
    # breakpoints that both element types hold exactly, in the other element type
    S = T === Float64 ? Float32 : Float64
    on8 = TensorSplineDistribution(
        PeriodicBSplineBasis(UniformMesh(16, zero(T) .. T(8)), 3), vb,
        landau_state(T))
    @test_throws other_basis VlasovPoisson(on8,
        Potential(PeriodicBSplineBasis(UniformMesh(16, zero(S) .. S(8)), 3)))
    # a quadratic Dirichlet basis on the breakpoints of a quadratic periodic basis has as many
    # functions, and PoissonSolverSpline accepts it
    xq2 = PeriodicBSplineBasis(UniformMesh(16, zero(T) .. T(4π)), 2)
    dirichlet = BSplineBasis(UniformMesh(16, zero(T) .. T(4π)), 2, Dirichlet())
    @test nbasis(dirichlet) == nbasis(xq2)
    @test_throws other_basis VlasovPoisson(
        TensorSplineDistribution(xq2, vb, landau_state(T)),
        Potential(dirichlet))
    model = VlasovPoisson(dist, Potential(xb))
    @test model.distribution === dist

    # the coefficient matrices of the forms, the field, the energy and the mass have the size of
    # the distribution's
    wrong = r"^DimensionMismatch: a coefficient matrix of the distribution has the size"
    A, B = dist.coefficients, zeros(T, 16, 99)
    @test_throws wrong _spline_1d1v_rhs!(similar(B), model, A)
    @test_throws wrong _spline_1d1v_rhs!(similar(A), model, B)
    @test_throws wrong _galerkin_form(dist, A, A, B)
    @test_throws wrong _galerkin_form(dist, B, A, A)
    @test_throws wrong _averaged_form(dist, A, B, A)
    @test_throws wrong _spline_energy(model, B)
    @test_throws wrong _spline_mass(model, B)

    # the method converts the span and the step, and copies the coefficients
    method = GalerkinMidpoint(model, (0.0, 1.0), 0.1)
    @test method isa VlasovMethods.VlasovMethod
    @test method.tspan === (zero(T), one(T))
    @test method.tstep === T(0.1)
    @test method.tol === 4eps(T)
    @test method.max_iterations == 10
    @test method.initial == dist.coefficients
    @test method.initial !== dist.coefficients
    m2 = GalerkinMidpoint(model, (0, 2), 1 // 4; tol = 1e-3, max_iterations = 0)
    @test m2.tspan === (zero(T), T(2)) && m2.tstep === T(0.25)
    @test m2.tol === T(1e-3) && m2.max_iterations == 0

    # non-finite coefficients
    for bad in (T(NaN), T(Inf), T(-Inf))
        C = copy(dist.coefficients)
        C[3, 7] = bad
        bad_model = VlasovPoisson(TensorSplineDistribution(xb, vb, C), Potential(xb))
        @test_throws r"^ArgumentError: .*finite" GalerkinMidpoint(bad_model, (0.0, 1.0), 0.1)
    end

    # the step, the tolerance and the iteration count
    for τ in (-0.1, 0, NaN, Inf, -Inf)
        @test_throws r"^ArgumentError: the step of GalerkinMidpoint" GalerkinMidpoint(
            model, (0.0, 1.0), τ)
    end
    # a step that overflows the element type, named as it was given
    T === Float32 &&
        @test_throws r"^ArgumentError: the step of GalerkinMidpoint .* got tstep = 1.0e39$" GalerkinMidpoint(
            model, (0.0, 1.0), 1e39)
    for tol in (-1, -eps(T), NaN)
        @test_throws r"^ArgumentError: the tolerance of GalerkinMidpoint" GalerkinMidpoint(
            model, (0.0, 1.0), 0.1; tol)
    end
    @test GalerkinMidpoint(model, (0.0, 1.0), 0.1; tol = 0).tol === zero(T)
    @test_throws r"^ArgumentError: the iteration count of GalerkinMidpoint" GalerkinMidpoint(
        model, (0.0, 1.0), 0.1; max_iterations = -1)
end

@testset "form ($T)" for T in (Float64, Float32)
    xb, vb = spline_bases(T)
    dist = TensorSplineDistribution(xb, vb, zeros(T, 16, 100))
    t8 = spline_tables(xb, vb; nq = 8)

    # the default quadrature is exact: it agrees with eight points per cell
    rng = Xoshiro(1)
    for _ in 1:10
        a, b, c = (randn(rng, T, 16, 100) for _ in 1:3)
        value = _galerkin_form(dist, a, b, c)
        @test value isa T
        ref = plain_form(t8, a, b, c)
        @test abs(value - ref) ≤ 16384eps(T) * abs(ref)
    end

    # the averaged form is totally antisymmetric, the plain form is not
    rng = Xoshiro(2)
    plain = T[]
    for _ in 1:20
        a, b, c = (randn(rng, T, 16, 100) for _ in 1:3)
        Tabc, Tbca, Tcab = _galerkin_form(dist, a, b, c), _galerkin_form(dist, b, c, a),
        _galerkin_form(dist, c, a, b)
        S = abs(Tabc) + abs(Tbca) + abs(Tcab)
        t̄ = _averaged_form(dist, a, b, c)
        @test t̄ isa T
        @test abs(t̄ + _averaged_form(dist, b, a, c)) ≤ 16eps(T) * S
        @test abs(t̄ + _averaged_form(dist, a, c, b)) ≤ 16eps(T) * S
        @test abs(t̄ + _averaged_form(dist, c, b, a)) ≤ 16eps(T) * S
        @test abs(t̄ - _averaged_form(dist, b, c, a)) ≤ 16eps(T) * S
        @test abs(t̄ - (Tabc + Tbca + Tcab) / 3) ≤ 16eps(T) * S
        push!(plain, abs(Tabc + _galerkin_form(dist, b, a, c)) / S)

        # on functions that vanish at ±V the plain form is antisymmetric in its first two slots
        for g in (a, b, c)
            g[:, 1] .= 0
            g[:, end] .= 0
        end
        Tabc, Tbac = _galerkin_form(dist, a, b, c), _galerkin_form(dist, b, a, c)
        @test abs(Tabc + Tbac) ≤ 1024eps(T) * (abs(Tabc) + abs(Tbac))
    end
    @test median(plain) ≥ 0.02
end

@testset "field ($T)" for T in (Float64, Float32)
    # the field is the averaged form tested against Ψ, for an even and an odd count of v-splines
    for nv in (100, 101)
        model, t = spline_model(T, pedestal_state(T); nv)
        dist = model.distribution
        A = dist.coefficients
        Ȧ = similar(A)
        @test _spline_1d1v_rhs!(Ȧ, model, A) === Ȧ
        @test eltype(Ȧ) === T
        ĥ = spline_hhat(t, A)
        for s in 1:5
            Ψ = randn(Xoshiro(s), T, 16, nv)
            sc = abs(_galerkin_form(dist, Ψ, ĥ, A)) + abs(_galerkin_form(dist, ĥ, A, Ψ)) +
                 abs(_galerkin_form(dist, A, Ψ, ĥ))
            @test abs(dot(Ψ, t.Mx * Ȧ * t.Mv) - _averaged_form(dist, Ψ, ĥ, A)) ≤
                  64eps(T) * sc
        end

        # the field is equivariant under the reflection v ↦ -v
        C = A .+ T(0.01) .* randn(Xoshiro(9), T, 16, nv)
        r = field(model, C)
        @test norm(field(model, reflect(C)) .+ reflect(r)) ≤ 256eps(T) * norm(r)
    end

    for (name, f) in (("Landau", landau_state(T)), ("pedestal", pedestal_state(T)))
        model, t = spline_model(T, f)
        A = model.distribution.coefficients
        Ȧ = field(model, A)

        # the energy rate vanishes; the plain field's does not where f reaches ±V
        ∇H = energy_gradient(t, A)
        sc = norm(∇H) * norm(t.Mx * Ȧ * t.Mv)
        @test abs(dot(∇H, Ȧ)) ≤ 8eps(T) * sc
        if name == "pedestal"
            @test abs(dot(∇H, plain_field(t, A))) ≥ 1e-6 * sc
        end

        # the energy and the mass are the test's own
        H = _spline_energy(model, A)
        @test H isa T
        @test abs(H - spline_H(t, A)) ≤ 64eps(T) * abs(H)
        M = _spline_mass(model, A)
        @test M isa T
        @test abs(M - spline_M(t, A)) ≤ 64eps(T) * abs(M)

        # the mass rate is the boundary flux (1/3) ∫ φ' (f(x, V) − f(x, −V)) dx
        rate = dot(t.ix, Ȧ * t.iv)
        flux = spline_flux(t, A)
        @test abs(rate - flux) ≤ 64eps(T) * abs(M) * norm(Ȧ) / norm(A)
        name == "pedestal" && @test abs(flux) ≥ 1e-3
    end

    # a homogeneous Maxwellian is stationary
    landau, _ = spline_model(T, landau_state(T))
    homogeneous, _ = spline_model(T, (x, v) -> exp(-v^2 / 2) / sqrt(T(2π)))
    @test norm(field(homogeneous, homogeneous.distribution.coefficients)) ≤
          1024eps(T) * norm(field(landau, landau.distribution.coefficients))
end

@testset "nonuniform x-mesh ($T)" for T in (Float64, Float32)
    # 16 cells of [0, 4π) of widths between 0.52 and 1.48 times the mean
    L = T(4π)
    bp = collect(range(zero(T), L; length = 17))
    bp[2:16] .+= (L / 16) .* T(1.25) .* sin.(T(2π) .* (1:15) ./ 16)
    xb = PeriodicBSplineBasis(GeneralMesh(bp), 3)
    _, vb = spline_bases(T)
    t = spline_tables(xb, vb)
    K = t.X1 * (t.wx .* t.X1')
    nonuniform(f) = VlasovPoisson(TensorSplineDistribution(xb, vb, f), Potential(xb))

    # a homogeneous Maxwellian has no field and is stationary
    landau = nonuniform(landau_state(T))
    homogeneous = nonuniform((x, v) -> exp(-v^2 / 2) / sqrt(T(2π)))
    H₀, A₀ = homogeneous.distribution.coefficients, landau.distribution.coefficients
    φ̂ = VlasovMethods._spline_potential(homogeneous, H₀)
    @test eltype(φ̂) === T
    @test norm(φ̂) ≤ 1024eps(T) * norm(VlasovMethods._spline_potential(landau, A₀))
    @test norm(field(homogeneous, H₀)) ≤ 1024eps(T) * norm(field(landau, A₀))

    for f in (landau_state(T), pedestal_state(T))
        model = nonuniform(f)
        A = model.distribution.coefficients
        # φ solves -φ'' = ρ - ρ̄ for the mean ρ̄ = ∫ρ dx / L, and ∫φ dx = 0
        b = t.Mx * (A * t.iv)
        b̃ = b .- (sum(b) / sum(t.ix)) .* t.ix
        φ̂ = VlasovMethods._spline_potential(model, A)
        @test norm(K * φ̂ - b̃) ≤ 64eps(T) * norm(b)
        @test abs(dot(t.ix, φ̂)) ≤ 64eps(T) * sum(t.ix) * norm(φ̂, Inf)

        # the energy is the kinetic energy and ½ ∫ (φ')² dx of the test's own potential
        φx = t.X1' * solve(t.solver, b̃)
        H = _spline_energy(model, A)
        @test H isa T
        @test abs(H - (dot(t.ix, A * t.kv) + sum(t.wx .* φx .^ 2) / 2)) ≤ 64eps(T) * abs(H)

        # the energy is quadratic in A, so its central difference along the field is its rate,
        # which vanishes
        Ȧ = field(model, A)
        ε = norm(A) / norm(Ȧ)
        H⁺, H⁻ = _spline_energy(model, A .+ ε .* Ȧ), _spline_energy(model, A .- ε .* Ȧ)
        @test abs(H⁺ - H⁻) ≤ 64eps(T) * (abs(H⁺) + abs(H⁻))
    end
end

@testset "residual and its derivative ($T)" for T in (Float64, Float32)
    model, _ = spline_model(T, pedestal_state(T))
    A = model.distribution.coefficients
    τ = T(0.05)
    y = vec(A .+ T(0.1) .* randn(Xoshiro(4), T, 16, 100) .* A)
    w = vec(randn(Xoshiro(5), T, 16, 100) .* A)
    F(y) = VlasovMethods._midpoint_residual!(similar(y), model, A, τ, y)

    # F(Y) = Y − A − τ Mx⁻¹ r(f̄) Mv⁻¹, with the field of the averaged form
    f̄ = (A .+ reshape(y, 16, 100)) ./ 2
    ref = y .- vec(A) .- τ .* vec(field(model, f̄))
    @test eltype(F(y)) === T
    @test norm(F(y) .- ref) ≤ 64eps(T) * norm(y)
    @test F(vec(A)) != zero(vec(A))

    # F is quadratic in y, so its central difference is its derivative up to round-off
    Jw = VlasovMethods._midpoint_jvp!(similar(w), model, A, τ, y, w)
    @test eltype(Jw) === T
    difference = (F(y .+ w) .- F(y .- w)) ./ 2
    @test norm(Jw .- difference) ≤ 64eps(T) * (norm(F(y .+ w)) + norm(F(y .- w)))
    @test norm(Jw .- w) ≥ 1e-4 * norm(w)
end

@testset "step ($T)" for T in (Float64, Float32)
    τ = T(0.05)
    bound = T === Float64 ? 1e-12 : 1e-4

    # the pedestal state reaches ±V, so its mass moves by the boundary flux
    model, t = spline_model(T, pedestal_state(T))
    A = copy(model.distribution.coefficients)
    method = GalerkinMidpoint(model, (0, 40τ), τ; tol = 64eps(T))
    _, z, times = run_states(method)
    @test size(z, 2) == 41 && length(times) == 41
    states = [reshape(z[:, n], 16, 100) for n in axes(z, 2)]
    @test states[1] == A
    @test states[2] != A
    H₀, M₀ = spline_H(t, A), spline_M(t, A)
    @test maximum(abs(spline_H(t, S) - H₀) / abs(H₀) for S in states) ≤ bound
    for n in 1:40
        Δ = spline_M(t, states[n + 1]) - spline_M(t, states[n])
        flux = spline_flux(t, (states[n] .+ states[n + 1]) ./ 2)
        @test abs(Δ - τ * flux) ≤ 64eps(T) * abs(M₀)
    end
    @test maximum(abs(spline_M(t, S) - M₀) / abs(M₀) for S in states) ≥ 1e-7

    # the plain field's own implicit midpoint step drifts in energy
    @test plain_midpoint_drift(t, A, 40, τ, 64eps(T)) ≥ 1e-4

    # one step and one step of the reflected state return to the start
    model, _ = spline_model(T, landau_state(T))
    A = copy(model.distribution.coefficients)
    A₁ = run!(GalerkinMidpoint(model, (0, τ), τ), tempname(; suffix = ".hdf5")).coefficients
    @test A₁ != A
    back = VlasovPoisson(
        TensorSplineDistribution(model.distribution.xbasis,
            model.distribution.vbasis, reflect(A₁)),
        model.potential)
    A₂ = run!(GalerkinMidpoint(back, (0, τ), τ), tempname(; suffix = ".hdf5")).coefficients
    @test norm(reflect(A₂) .- A) ≤ 64eps(T) * norm(A)
end

@testset "stopping test ($T)" for T in (Float64, Float32)
    τ = T(0.05)
    model, _ = spline_model(T, landau_state(T))
    A = copy(model.distribution.coefficients)
    F(y) = VlasovMethods._midpoint_residual!(similar(y), model, A, τ, y)

    # the solve stops at ‖F‖₂ ≤ tol ‖A‖₂: at half the residual of the start the target binds
    # after the first Newton iteration, and a target looser by a factor of 2 accepts the start
    tol = norm(F(vec(A))) / (2 * norm(A))
    Y = run!(GalerkinMidpoint(model, (0, τ), τ; tol), tempname(; suffix = ".hdf5")).coefficients
    @test eltype(Y) === T
    @test Y != A
    @test norm(F(vec(Y))) ≤ tol * norm(A)

    # with tol = 0 only the iteration count ends the solve, because the step test is off
    model, _ = spline_model(T, landau_state(T))
    err = try
        run!(GalerkinMidpoint(model, (0, τ), τ; tol = 0, max_iterations = 8),
            tempname(; suffix = ".hdf5"))
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("returned MAXITERS", err.msg) && occursin("after 8 iterations", err.msg)
end

@testset "zero state ($T)" for T in (Float64, Float32)
    model, _ = spline_model(T, (x, v) -> zero(x))
    @test iszero(model.distribution.coefficients)
    _, z, _ = run_states(GalerkinMidpoint(model, (0.0, 0.15), 0.05))
    @test size(z, 2) == 4
    @test all(iszero, z)
end

@testset "failure ($T)" for T in (Float64, Float32)
    model, _ = spline_model(T, landau_state(T))
    A = copy(model.distribution.coefficients)
    method = GalerkinMidpoint(model, (0.0, 1.0), 0.05; max_iterations = 0)
    file = tempname(; suffix = ".hdf5")
    # the garbage collector is held off, so that no finalizer closes a handle that run! leaves
    # open; HDF5 refuses to truncate a file that is open
    GC.gc()
    before = open_hdf5_files()
    GC.enable(false)
    err = try
        run!(method, file)
    catch e
        e
    end
    after = open_hdf5_files()
    initial = h5open(file, "r") do h5
        read(h5["z"])[:, 1]
    end
    reopened = try
        close(h5open(file, "w"))
        true
    catch
        false
    end
    GC.enable(true)
    @test err isa ErrorException
    @test occursin("MAXITERS", err.msg)
    @test occursin("‖F‖", err.msg) && occursin("iterations", err.msg)
    # the file is closed and holds the initial state, and the distribution keeps its coefficients
    @test after == before
    @test reopened
    @test initial == vec(A)
    @test model.distribution.coefficients == A
end

@testset "run ($T)" for T in (Float64, Float32)
    model, _ = spline_model(T, landau_state(T))
    A = copy(model.distribution.coefficients)
    # the residual floor of a step of 0.3 is about 4 eps(T) ‖A‖ in Float64, the default tolerance
    dist, z, times = run_states(GalerkinMidpoint(model, (0.0, 1.0), 0.3; tol = 64eps(T)))
    @test size(z) == (16 * 100, 4)
    @test times == collect(T(0):T(0.3):T(1))
    @test eltype(z) === T && eltype(times) === T && eltype(dist.coefficients) === T
    @test z[:, 1] == vec(A)
    @test z[:, 2] != z[:, 1]
    @test dist === model.distribution
    @test vec(dist.coefficients) == z[:, end]

    # two runs of one method start from the same state
    method = GalerkinMidpoint(model, (0.0, 0.25), 0.05)
    _, z₁, _ = run_states(method)
    _, z₂, _ = run_states(method)
    @test size(z₁, 2) == 6
    @test z₁ == z₂
end
