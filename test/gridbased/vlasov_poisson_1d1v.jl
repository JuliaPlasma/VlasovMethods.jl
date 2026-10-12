using GeometricEquations: ntime
using GeometricIntegrators: GeometricIntegrators
using HDF5: h5read
using LinearAlgebra: dot, issymmetric
using PoissonSolvers: FFTWBasis, FiniteDifferenceBasis, PeriodicBasisSpline, Potential,
                      PoissonSolverFFT, PoissonSolverMatrixFree, solve
using Test
using VlasovMethods

const VM = VlasovMethods

# The diagnostics of every step of a run of `int`, which builds the method the constructor
# defaults to when `march(model, tspan, tstep)` is called. `H` is the discrete Hamiltonian
# `hx*hv*Σ f v²/2 + (hx/2)*dot(ρ, ϕ)` with `ϕ` the backend's field for `ρ = hv*Σ f`, computed here
# and not in `src`, so that a sign or a scale flip in the right-hand side cannot be matched by an
# equally flipped energy. `Σ` is `hx*hv*Σf`, `L¹` is `hx*hv*Σ|f|`, `negmass` the negative mass of
# the L¹ bound, and `edge` the velocity edge in units of `eps(T)*maximum(abs, f)`.
function march(int)
    dist = int.model.distribution
    T = eltype(dist.values)
    nx, nv = size(dist.values)
    hx, hv = step(dist.x), step(dist.v)
    v = dist.v
    solver = int.model.potential.solver

    function energy(F)
        ρ = hv .* vec(sum(F; dims = 2))
        return hx * hv * sum(F[i, j] * v[j]^2 / 2 for j in 1:nv, i in 1:nx) +
               (hx / 2) * dot(ρ, solve(solver, ρ))
    end
    edge(F) = max(maximum(abs, @view(F[:, 1])), maximum(abs, @view(F[:, nv]))) /
              (eps(T) * maximum(abs, F))

    sol = GeometricIntegrators.integrate(int.integrator)
    Fs = [reshape(sol.q[n], nx, nv) for n in 0:ntime(int.equation)]

    return (H = energy.(Fs), Σ = [hx * hv * sum(F) for F in Fs],
        L¹ = [hx * hv * sum(abs, F) for F in Fs],
        negmass = [hx * hv * sum(fᵢ -> max(-fᵢ, zero(T)), F) for F in Fs],
        minf = minimum.(Fs), edge = edge.(Fs))
end

march(model, tspan, tstep, method) = march(GeometricIntegrator(model, tspan, tstep, method))
march(model, tspan, tstep) = march(GeometricIntegrator(model, tspan, tstep))

# The initial state of the conservation runs, on the scheme's own grid: `x ∈ [0, 4π)` with `nx`
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
        # the grid facts, and the stencil they are read through. `f = cos(k x_i) v_j` has a
        # vanishing density, so `ϕ = 0` and `h = v²/2` exactly, and the expected right-hand side is
        # the transcribed Arakawa stencil at every row, rows 1 and nv included.
        nx, nv, k = 8, 12, T(0.5)
        xd = (zero(T), T(4π))
        vd = (T(-5), T(5))
        dist = GridDistribution(zeros(T, nx, nv), xd, vd)
        x, v = dist.x, dist.v
        hx, hv = step(x), step(v)

        # the grid of `GridDistribution`, not the driver's `range(0, 1, length = nx)` with
        # `h₁ = x[2] - x[1]`, which stores the periodic endpoint twice
        @test hx == T(4π) / nx
        @test hv == (vd[2] - vd[1]) / (nv - 1)
        @test first(x) == zero(T)
        @test last(x) == T(4π) - hx
        @test first(v) == vd[1]
        @test last(v) == vd[2]

        dist.values .= cos.(k .* x) .* v'
        model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))

        int = GeometricIntegrator(model, (zero(T), one(T)), one(T))
        prob = int.equation
        q = copy(prob.ics.q)
        ḟ = similar(q)
        prob.equation.v(ḟ, zero(T), q, prob.parameters)
        # the vector field of every Newton iteration is type-stable
        @test @inferred(VM._grid_1d1v_rhs!(ḟ, zero(T), q, prob.parameters)) === ḟ

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
        # under the (−1, 16, −30, 16, −1)/12h² stencil, and by 4 under the second-order
        # (5, −32, 54, −32, 5)/12h² one.
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
        # one conserving run, read by all three checks: the horizon `T` and the step `dt` are the
        # row's, the grid is the scheme's, and the velocity range keeps `f` below
        # `eps(T)·maximum(f)` at both ends
        nx, nv = 16, 64
        dist = maxwellian_grid(T, nx, nv, vlim, α)
        model = VlasovPoisson(dist, Potential(FFTWBasis((zero(T), T(4π)), nx)))
        tspan = (zero(T), T(nsteps) * T(dt))

        res = march(model, tspan, T(dt))

        # the discrete Hamiltonian, and the explicit midpoint of the same order at the same step,
        # whose drift must be at least ten times the bound
        drift = abs(res.H[end] - res.H[1]) / abs(res.H[1])
        bound = nsteps * 2eps(T)
        explicit = march(model, tspan, T(dt), GeometricIntegrators.ExplicitMidpoint())
        driftₑ = abs(explicit.H[end] - explicit.H[1]) / abs(explicit.H[1])

        @test drift ≤ bound
        @test driftₑ ≥ 10 * bound

        # the mass, and the L¹ norm `hx*hv*Σ|f|` rather than `hx*hv*Σf`. Arakawa keeps L¹ only up
        # to twice the negative mass, which is the second term. The identity below is what
        # separates the two sums: `|f| = f + 2*max(-f, 0)` at every node, so a run whose L¹ was
        # `Σf` would miss its right-hand side by twice this negative mass.
        mass_bound = 10 * nsteps * eps(T)
        @test abs(res.Σ[end] - res.Σ[1]) / abs(res.Σ[1]) ≤ mass_bound
        @test abs(res.L¹[end] - res.L¹[1]) / abs(res.L¹[1]) ≤
              mass_bound + 2 * res.negmass[end] / abs(res.L¹[1])
        @test isapprox(res.L¹[end], res.Σ[end] + 2 * res.negmass[end];
            rtol = 4 * nx * nv * eps(T))

        # the scheme is periodic in `v`, so a run is valid only while `f` vanishes at `v_a` and `v_b`;
        # the bound is asserted at every saved step
        @test maximum(res.edge) ≤ 1

        diagnostics = (; T_horizon = tspan[2], n_steps = nsteps, drift, driftₑ, bound,
            minf = minimum(res.minf), negmass = res.negmass[end], edge = maximum(res.edge))
        @info "checks 1 to 3" T diagnostics
    end

    @testset "the constructor's guards" for T in (Float32, Float64)
        nx, nv = 8, 12
        xd = (zero(T), T(4π))
        vd = (T(-5), T(5))
        dist = GridDistribution(zeros(T, nx, nv), xd, vd)

        # a spline basis samples ρ at the nodes rather than forming a Galerkin load vector
        @test_throws ArgumentError VlasovPoisson(dist, Potential(PeriodicBasisSpline(xd, 4, nx)))
        # one degree of freedom per x-node, and no more
        @test_throws DimensionMismatch VlasovPoisson(dist, Potential(FFTWBasis(xd, nx + 1)))
        # the potential's own x-nodes, compared exactly
        @test_throws ArgumentError VlasovPoisson(dist,
            Potential(FFTWBasis((zero(T), T(4π) + T(0.1)), nx)))
        # fewer than three nodes per direction is Arakawa's own check, raised at construction
        @test_throws ArgumentError VlasovPoisson(GridDistribution(zeros(T, nx, 2), xd, vd),
            Potential(FFTWBasis(xd, nx)))
        @test_throws ArgumentError VlasovPoisson(GridDistribution(zeros(T, 2, nv), xd, vd),
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
