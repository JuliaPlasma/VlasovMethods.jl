import GeometricEquations
using GeometricEquations: ntime
using GeometricIntegrators: GeometricIntegrators
using HDF5: h5read
using LinearAlgebra: dot, issymmetric
using PoissonSolvers: FFTWBasis, FiniteDifferenceBasis, PeriodicBasisSpline, Potential,
                      PoissonSolverFFT, PoissonSolverMatrixFree, solve
using Test
using VlasovMethods

const VM = VlasovMethods

# The scheme's own stepping loop, without the HDF5 file, so that every saved step can be read:
# `run!(::GeometricIntegrator, h5file)` advances the same solution step the same way. The two
# `march(model, …)` forms build the integrator, so that the three-argument one exercises the
# method the constructor defaults to.
#
# `H` is the discrete Hamiltonian of the scheme, `hx*hv*Σ f v²/2 + (hx/2)*dot(ρ, ϕ)` with `ϕ` the
# field the backend returns for `ρ = hv*Σ f`. It is computed here and not in `src`, so that a sign
# or a scale flip in the right-hand side cannot be matched by an equally flipped energy. `Σ` and
# `L¹` are the two sums check 2 separates, `edge` the velocity edge of check 3 in units of
# `eps(T)*maximum(abs, f)`, and `negmass` the negative mass of the L¹ bound.
function march(int)
    prob = int.equation
    dist = int.model.distribution
    T = eltype(dist.values)
    nx, nv = size(dist.values)
    hx, hv = step(dist.x), step(dist.v)
    v = dist.v
    solver = int.model.potential.solver

    function energy(f)
        F = reshape(f, nx, nv)
        ρ = hv .* vec(sum(F; dims = 2))
        return hx * hv * sum(F .* (v' .^ 2 ./ 2)) + (hx / 2) * dot(ρ, solve(solver, ρ))
    end

    times = GeometricEquations.initialtime(prob):GeometricEquations.timestep(prob):GeometricEquations.finaltime(prob)
    n = ntime(prob)
    solstep = GeometricIntegrators.GeometricIntegratorsBase.solutionstep(
        int.integrator, GeometricEquations.initialstate(prob))

    H = Vector{T}(undef, n + 1)
    Σ = Vector{T}(undef, n + 1)
    L¹ = Vector{T}(undef, n + 1)
    negmass = Vector{T}(undef, n + 1)
    minf = Vector{T}(undef, n + 1)
    edge = Vector{T}(undef, n + 1)

    function record!(s, f)
        F = reshape(f, nx, nv)
        H[s] = energy(f)
        Σ[s] = hx * hv * sum(F)
        L¹[s] = hx * hv * sum(abs, F)
        negmass[s] = hx * hv * sum(max.(-F, zero(T)))
        minf[s] = minimum(F)
        edge[s] = max(maximum(abs, F[:, 1]), maximum(abs, F[:, nv])) /
                  (eps(T) * maximum(abs, F))
    end

    record!(1, solstep.q)
    for s in 1:n
        GeometricIntegrators.reset!(solstep, times[s + 1])
        GeometricIntegrators.integrate!(solstep, int.integrator)
        record!(s + 1, solstep.q)
    end

    return (H = H, Σ = Σ, L¹ = L¹, negmass = negmass, minf = minf, edge = edge)
end

march(model, tspan, tstep, method) = march(GeometricIntegrator(model, tspan, tstep, method))
march(model, tspan, tstep) = march(GeometricIntegrator(model, tspan, tstep))

# The initial state of the conservation runs, on the grid of check 4: `x ∈ [0, 4π)` with `nx`
# nodes and the periodic endpoint not stored, `v ∈ [-vlim, vlim]` with `nv` nodes and both ends
# nodes. `α` and `vlim` are the ones of the run's row.
function maxwellian_grid(T, nx, nv, vlim, α; k = T(0.5))
    xd = (zero(T), T(4π))
    vd = (T(-vlim), T(vlim))
    dist = GridDistribution(zeros(T, nx, nv), xd, vd)
    dist.values .= (1 .+ T(α) .* cos.(k .* dist.x)) .* exp.(-dist.v' .^ 2 ./ 2) ./
                   sqrt(T(2π))
    return dist
end

@testset "Grid 1D1V Vlasov–Poisson" begin
    @testset "the exact right-hand side, $T" for T in (Float32, Float64)
        # check 4's grid facts, and the stencil they are read through. `f = cos(k x_i) v_j` has a
        # vanishing density, so `ϕ = 0` and `h = v²/2` exactly, and the expected right-hand side is
        # the transcribed Arakawa stencil at every row, rows 1 and nv included.
        nx, nv, k = 8, 12, T(0.5)
        xd = (zero(T), T(4π))
        vd = (T(-5), T(5))
        x = range(xd[1], xd[2]; length = nx + 1)[1:nx]
        v = range(vd[1], vd[2]; length = nv)
        hx, hv = step(x), step(v)

        # the grid of `GridDistribution`, not the driver's `range(0, 1, length = nx)` with
        # `h₁ = x[2] - x[1]`, which stores the periodic endpoint twice
        @test hx == T(4π) / nx
        @test hv == (vd[2] - vd[1]) / (nv - 1)
        @test first(x) == zero(T)
        @test last(x) == T(4π) - hx
        @test first(v) == vd[1]
        @test last(v) == vd[2]

        dist = GridDistribution([cos(k * x[i]) * v[j] for i in 1:nx, j in 1:nv], xd, vd)
        model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))

        int = GeometricIntegrator(model, (zero(T), one(T)), one(T))
        prob = int.equation
        q = copy(prob.ics.q)
        ḟ = similar(q)
        prob.equation.v(ḟ, zero(T), q, prob.parameters)

        c = cos.(k .* x)
        H = v .^ 2 ./ 2
        expected = zeros(T, nx, nv)
        for i in 1:nx, j in 1:nv

            i₋, i₊ = mod1(i - 1, nx), mod1(i + 1, nx)
            j₋, j₊ = mod1(j - 1, nv), mod1(j + 1, nv)
            Δc = c[i₊] - c[i₋]
            expected[i, j] = -Δc * (2 * v[j] * (H[j₊] - H[j₋]) +
                              v[j₊] * (H[j₊] - H[j]) +
                              v[j₋] * (H[j] - H[j₋])) / (12 * hx * hv)
        end

        atol = 8nv * eps(T) * maximum(abs, expected)
        @test maximum(abs, reshape(ḟ, nx, nv) .- expected) ≤ atol
        # the interior rows reduce to −Δcᵢ (vⱼ² + hv²/6) / (2hx)
        @test expected[4, 6] ≈ -T(1) * (c[5] - c[3]) * (v[6]^2 + hv^2 / 6) / (2 * hx) atol = atol
    end

    @testset "the Green's matrix, $T" for T in (Float32, Float64)
        nx, nv, k = 16, 12, T(0.5)
        xd = (zero(T), T(4π))
        vd = (T(-5), T(5))
        dist = GridDistribution(zeros(T, nx, nv), xd, vd)
        potential = Potential(FFTWBasis(xd, nx))
        model = VlasovPoisson(dist, potential)

        @test model.work isa Matrix{T}
        @test issymmetric(model.work)

        x = dist.x
        ρ = 1 .+ cos.(k .* x) .+ T(0.3) * sin.(2k .* x)
        @test model.work * ρ ≈ solve(potential.solver, ρ) rtol = 16nx * eps(T)
    end

    @testset "the fourth-order matrix-free backend" begin
        # Float64 only. The error of G·cos(kx) against cos(kx)/k² falls by 16 from nx = 16 to 32
        # under the (−1, 16, −30, 16, −1)/12h² stencil, and by 4 under the old
        # (5, −32, 54, −32, 5)/12h² one, which P5 replaced.
        k = 0.5
        errs = Float64[]
        for nx in (16, 32)
            xd = (0.0, 4π)
            vd = (-5.0, 5.0)
            dist = GridDistribution(zeros(nx, 4), xd, vd)
            model = VlasovPoisson(dist, Potential(FiniteDifferenceBasis(xd, nx; order = 4)))
            x = dist.x
            push!(errs, maximum(abs, model.work * cos.(k .* x) .- cos.(k .* x) ./ k^2))
        end
        @test errs[1] / errs[2] ≥ 12
    end

    @testset "checks 1 to 3: conservation, $T" for (T, vlim, α, dt, nsteps) in ((
        Float64, 10.0, 0.01, 0.05, 40), (Float32, 7.5, 0.5, 0.2, 14))
        # one conserving run, read by all three checks: the horizon `T` and the step `dt` are
        # Q5's, the grid is check 4's, and the velocity range keeps `f` below `eps(T)·maximum(f)`
        # at both ends (check 3, Q7)
        nx, nv = 16, 64
        dist = maxwellian_grid(T, nx, nv, vlim, α)
        model = VlasovPoisson(dist, Potential(FFTWBasis((zero(T), T(4π)), nx)))
        tspan = (zero(T), T(nsteps) * T(dt))

        run = march(model, tspan, T(dt))

        # check 1: the discrete Hamiltonian, and the explicit midpoint of the same order at the
        # same step, whose drift must be at least ten times the bound
        drift = abs(run.H[end] - run.H[1]) / abs(run.H[1])
        bound = nsteps * 2eps(T)
        explicit = march(model, tspan, T(dt), GeometricIntegrators.ExplicitMidpoint())
        driftₑ = abs(explicit.H[end] - explicit.H[1]) / abs(explicit.H[1])

        @test drift ≤ bound
        @test driftₑ ≥ 10 * bound
        @test driftₑ ≥ 10 * drift

        # check 2: the mass, and the L¹ norm `hx*hv*Σ|f|` rather than `hx*hv*Σf`. Arakawa keeps
        # L¹ only up to twice the negative mass, which is the second term; a test that used `Σf`
        # here would compare the mass drift against this bound and PINS.
        mass_bound = 10 * nsteps * eps(T)
        @test abs(run.Σ[end] - run.Σ[1]) / abs(run.Σ[1]) ≤ mass_bound
        @test abs(run.L¹[end] - run.L¹[1]) / abs(run.L¹[1]) ≤
              mass_bound + 2 * run.negmass[end] / abs(run.L¹[1])

        # check 3: the torus identities of the first two hold only while `f` vanishes at `v_a`
        # and `v_b`, which is asserted at every saved step
        @test maximum(run.edge) ≤ 1

        diagnostics = (; T_horizon = tspan[2], n_steps = nsteps, drift, driftₑ, bound,
            minf = minimum(run.minf), negmass = run.negmass[end], edge = maximum(run.edge))
        @info "checks 1 to 3" T diagnostics
    end

    @testset "the constructor's guards" begin
        nx, nv = 8, 12
        xd = (0.0, 4π)
        vd = (-5.0, 5.0)
        dist = GridDistribution(zeros(nx, nv), xd, vd)

        # a spline basis samples ρ at the nodes rather than forming a Galerkin load vector
        @test_throws ArgumentError VlasovPoisson(dist, Potential(PeriodicBasisSpline(xd, 4, nx)))
        # one degree of freedom per x-node, and no more
        @test_throws DimensionMismatch VlasovPoisson(dist, Potential(FFTWBasis(xd, nx + 1)))
        # the potential's own x-nodes, compared exactly
        @test_throws ArgumentError VlasovPoisson(dist,
            Potential(FFTWBasis((0.0, 4π + 0.1), nx)))
        # fewer than three velocity nodes is Arakawa's own check, raised at construction
        @test_throws ArgumentError VlasovPoisson(GridDistribution(zeros(nx, 2), xd, vd),
            Potential(FFTWBasis(xd, nx)))
        @test_throws ArgumentError VlasovPoisson(GridDistribution(zeros(2, nv), xd, vd),
            Potential(FFTWBasis(xd, 2)))
    end

    @testset "run! writes the grid state back" begin
        nx, nv = 8, 12
        xd = (0.0, 4π)
        dist = maxwellian_grid(Float64, nx, nv, 5.0, 0.01)
        f₀ = copy(vec(dist.values))
        model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))

        method = GeometricIntegrator(model, (0.0, 2 * 0.05), 0.05)
        h5file = joinpath(mktempdir(), "grid_vlasov_poisson.hdf5")
        @test run!(method, h5file) === model.distribution

        z = h5read(h5file, "z")
        @test size(z) == (nx * nv, ntime(method.equation) + 1)
        @test z[:, 1] == f₀
        @test z[:, end] != z[:, 1]
        @test model.distribution.values == reshape(z[:, end], nx, nv)
    end
end
