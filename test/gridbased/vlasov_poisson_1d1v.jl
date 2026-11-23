using ForwardDiff: ForwardDiff
using GeometricBrackets: poisson_apply
using GeometricEquations: ntime
using GeometricIntegrators: GeometricIntegrators
using HDF5: h5read
using LinearAlgebra: dot, issymmetric, norm
using PoissonSolvers: FFTWBasis, FiniteDifferenceBasis, PeriodicBasisSpline, Potential,
                      PoissonSolverFFT, PoissonSolverMatrixFree, solve
using Random: Xoshiro
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

# The right-hand side written out with the allocating bracket `poisson_apply`, on the bracket `a`,
# the velocity nodes `v` and the Green's matrix `G`: `ρ = hv Σⱼ f`, `ϕ = G ρ`, `h = ϕ .+ v²/2`
# and `ḟ = poisson_apply(a, f, h) / (hx hv)`.
function allocating_rhs(f, a, v, G)
    ϕ = G * (a.hv .* vec(sum(reshape(f, a.nx, a.nv); dims = 2)))
    h = vec(ϕ .+ v' .^ 2 ./ 2)
    return poisson_apply(a, f, h) ./ (a.hx * a.hv)
end

# The bytes of one call of `rhs!` after a warm call, through a barrier with concrete arguments.
function rhs_allocations(rhs!::R, ḟ::D, t::T, f::F, params::P) where {R, D, T, F, P}
    rhs!(ḟ, t, f, params)
    return @allocated rhs!(ḟ, t, f, params)
end

# A `Dual{Nothing, T, 12}` copy of `q`, the chunk of the Jacobian of `ImplicitMidpoint`.
function dual12(q::AbstractVector{T}) where {T}
    ForwardDiff.Dual{Nothing}.(q, ntuple(_ -> one(T), 12)...)
end

# The Jacobian of `allocating_rhs` at `f`, column by column. The bracket is bilinear and `h` is
# affine in `f`, so the column of the unit vector `e` is
# `(poisson_apply(a, e, h(f)) + poisson_apply(a, f, δh(e))) / (hx hv)`, with
# `h(f) = vec(G ρ(f) .+ v'.^2 / 2)`, `δh(e) = vec(G ρ(e) .+ 0 .* v')` and `ρ(f) = hv Σⱼ f`.
function reference_jacobian(f, a, v, G)
    ρ(g) = a.hv .* vec(sum(reshape(g, a.nx, a.nv); dims = 2))
    h = vec(G * ρ(f) .+ v' .^ 2 ./ 2)
    J = similar(f, length(f), length(f))
    for k in eachindex(f)
        e = zero(f)
        e[k] = one(eltype(f))
        δh = vec(G * ρ(e) .+ zero(eltype(f)) .* v')
        J[:, k] .= (poisson_apply(a, e, h) .+ poisson_apply(a, f, δh)) ./ (a.hx * a.hv)
    end
    return J
end

# The model and the problem of the seeded state, at `nx, nv = 16, 64` unless named: the perturbed
# Maxwellian of the `Float32` conservation run, whose density does not vanish, so `ϕ ≠ 0`.
function seeded_problem(T; nx = 16, nv = 64)
    dist = maxwellian_grid(T, nx, nv, 7.5, 0.5)
    model = VlasovPoisson(dist, Potential(FFTWBasis((zero(T), T(4π)), nx)))
    return model, GeometricIntegrator(model, (zero(T), one(T)), T(0.1)).equation
end

# The bounded velocity closure, written out on the test's side. The unknowns are the nodes
# `j = 2 … nv-1`; `f` and the Hamiltonian are read as zero on and beyond the end nodes, and the
# field is `ḟ = χ ∘ [Ph, χ ∘ f]` with `Ph = χ ∘ (ϕ ⊗ 1 + 1 ⊗ v²/2)`, `ϕ = G ρ` and
# `ρ = hv Σ_{j=2}^{nv-1} f`. The oracle is the field of the cyclic average `T̄` of Arakawa's form,
# built column by column from a transcription of the stencil of its own.

# The tolerance `1e-12` in `Float64` is 4504 eps(Float64); the same budget in units of eps(T),
# rounded down to a power of 2, gives `4096 eps(Float32)` in `Float32`.
tolerance(::Type{Float64}) = 1e-12
tolerance(::Type{Float32}) = 4096 * eps(Float32)

rel(a, b) = norm(a .- b) / norm(b)

# standard-normal draws in Float64, rounded to T
draw(rng, T, dims...) = T.(randn(rng, dims...))

# `a` with its columns 1 and nv set to zero, by assignment, so that a NaN there does not survive
function cut(a::AbstractMatrix)
    b = copy(a)
    b[:, [1, end]] .= 0
    return b
end

# a subtype of `VelocityClosure` that is neither closure, which `GeometricIntegrator` refuses
struct OtherClosure <: VM.VelocityClosure end

# the model and the problem of the bounded closure on [0, 4π) × [-6, 6]
function bounded_problem(T, nx, nv)
    xd = (zero(T), T(4π))
    dist = GridDistribution(zeros(T, nx, nv), xd, (T(-6), T(6)))
    model = VlasovPoisson(dist, Potential(FFTWBasis(xd, nx)))
    int = GeometricIntegrator(model, (zero(T), one(T)), T(0.1);
        closure = BoundedVelocityClosure())
    return model, int.equation
end

# the vector field of the problem `prob` applied to the matrix `f`, as a matrix
function applied(prob, f)
    a = prob.parameters.arakawa
    ḟ = similar(vec(f))
    @test prob.equation.v(ḟ, zero(eltype(f)), vec(copy(f)), prob.parameters) === ḟ
    return reshape(ḟ, a.nx, a.nv)
end

# [a, b] on the periodic grid, from GeometricBrackets, which gives hx hv [c, û]
function bracket(p, a, b)
    A = p.arakawa
    return reshape(poisson_apply(A, vec(b), vec(a)) ./ (A.hx * A.hv), A.nx, A.nv)
end

bounded_potential(p, f) = p.G * (p.arakawa.hv .* vec(sum(cut(f); dims = 2)))
cut_h(p, f) = cut(bounded_potential(p, f) .+ p.kinetic')
uncut_h(p, f) = bounded_potential(p, f) .+ p.kinetic'
bounded_field(p, f) = cut(bracket(p, cut_h(p, f), cut(f)))
# the boundary flux −hx hv Σ_i Σ_{j ∈ {1, nv}} [Ph, χ ∘ f]_ij
function boundary_flux(p, f)
    B = bracket(p, cut_h(p, f), cut(f))
    return -p.arakawa.hx * p.arakawa.hv * (sum(B[:, 1]) + sum(B[:, end]))
end
# ⟨a, b⟩ over the unknowns
ip(p, a, b) = p.arakawa.hx * p.arakawa.hv * sum((a .* b)[:, 2:(end - 1)])

# Arakawa's stencil [F, H] at (i, j), on index readers F and H
function arakawa_at(F, H, i, j, hx, hv)
    J1 = (F(i + 1, j) - F(i - 1, j)) * (H(i, j + 1) - H(i, j - 1)) -
         (F(i, j + 1) - F(i, j - 1)) * (H(i + 1, j) - H(i - 1, j))
    J2 = F(i + 1, j) * (H(i + 1, j + 1) - H(i + 1, j - 1)) -
         F(i - 1, j) * (H(i - 1, j + 1) - H(i - 1, j - 1)) -
         F(i, j + 1) * (H(i + 1, j + 1) - H(i - 1, j + 1)) +
         F(i, j - 1) * (H(i + 1, j - 1) - H(i - 1, j - 1))
    J3 = F(i + 1, j + 1) * (H(i, j + 1) - H(i + 1, j)) -
         F(i - 1, j - 1) * (H(i - 1, j) - H(i, j - 1)) -
         F(i - 1, j + 1) * (H(i, j + 1) - H(i - 1, j)) +
         F(i + 1, j - 1) * (H(i + 1, j) - H(i, j - 1))
    return (J1 + J2 + J3) / (12hx * hv)
end
# a reader that wraps in x and returns zero outside the unknowns
reader(a) = (i, j) -> (2 ≤ j ≤ size(a, 2) - 1) ? a[mod1(i, size(a, 1)), j] : zero(eltype(a))
function stencil_bracket(p, a, b)
    A = p.arakawa
    out = zeros(promote_type(eltype(a), eltype(b)), A.nx, A.nv)
    for i in 1:(A.nx), j in 2:(A.nv - 1)

        out[i, j] = arakawa_at(reader(a), reader(b), i, j, A.hx, A.hv)
    end
    return out
end
T3(p, a, b, c) = ip(p, a, stencil_bracket(p, b, c))
Tbar(p, a, b, c) = (T3(p, a, b, c) + T3(p, b, c, a) + T3(p, c, a, b)) / 3

# the field of T̄: ⟨ψ, ḟ⟩ = T̄(ψ, Ph, f) for the unit vectors ψ of the unknowns
function tbar_field(p, f)
    A = p.arakawa
    Ph = cut_h(p, f)
    out = zeros(eltype(f), A.nx, A.nv)
    for i in 1:(A.nx), j in 2:(A.nv - 1)

        e = zeros(eltype(f), A.nx, A.nv)
        e[i, j] = 1
        out[i, j] = Tbar(p, e, Ph, f) / (A.hx * A.hv)
    end
    return out
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

    @testset "the buffers of the right-hand side, $T" for T in (Float32, Float64)
        model, prob = seeded_problem(T)
        dist, params = model.distribution, prob.parameters
        nx, nv = size(dist.values)

        @test params.kinetic isa Vector{T}
        @test params.kinetic == collect(dist.v) .^ 2 ./ 2
        @test params.ρ isa Vector{T} && length(params.ρ) == nx
        @test params.ϕ isa Vector{T} && length(params.ϕ) == nx
        @test params.h isa Vector{T} && length(params.h) == nx * nv
    end

    @testset "the in-place right-hand side agrees with the allocating one, $T" for T in (
        Float32, Float64)
        # two states of different densities in sequence on one `params`, so that a buffer that
        # keeps the value of the previous call is seen. The bound is the one of the exact testset.
        model, prob = seeded_problem(T)
        params = prob.parameters
        nv = size(model.distribution.values, 2)
        q = Vector(prob.ics.q)
        for f in (q, q .* (1 .+ T(0.25) .* sin.(T.(1:length(q)))))
            ḟ = similar(f)
            @test VM._grid_1d1v_rhs!(ḟ, zero(T), f, params) === ḟ
            ref = allocating_rhs(f, params.arakawa, params.v, params.G)
            @test maximum(abs, ḟ .- ref) ≤ 8nv * eps(T) * maximum(abs, ref)
        end
    end

    @testset "a state of another length, $T" for T in (Float32, Float64)
        # the guard refuses before any buffer is written; `reshape` and `poisson_apply!` raise
        # the same type on the same input, so the phrase is the evidence
        _, prob = seeded_problem(T)
        params = prob.parameters
        q = Vector(prob.ics.q)
        ḟ = similar(q)
        VM._grid_1d1v_rhs!(ḟ, zero(T), q, params)
        before = [copy(x) for x in values(params) if x isa Array]
        ḟ₀ = copy(ḟ)

        short = q[1:(15 * 64)]
        @test_throws DimensionMismatch("_grid_1d1v_rhs! on a 16 × 64 grid takes a state f " *
                                       "and a derivative ḟ of nx * nv values each; got 960 and 1024") VM._grid_1d1v_rhs!(
            ḟ, zero(T), short, params)
        @test [x for x in values(params) if x isa Array] == before
        @test ḟ == ḟ₀

        long = similar(q, 17 * 64)
        @test_throws DimensionMismatch("_grid_1d1v_rhs! on a 16 × 64 grid takes a state f " *
                                       "and a derivative ḟ of nx * nv values each; got 1024 and 1088") VM._grid_1d1v_rhs!(
            long, zero(T), q, params)
        @test [x for x in values(params) if x isa Array] == before
    end

    @testset "params of another form take the allocating method, $T" for T in (Float32, Float64)
        # `params` without the buffers, as `(arakawa, v, G)`, has a bracket of element type `T`
        # but no buffers to write into
        _, prob = seeded_problem(T)
        params = prob.parameters
        base = (; params.arakawa, params.v, params.G)
        q = Vector(prob.ics.q)
        ḟ = similar(q)
        @test VM._grid_1d1v_rhs!(ḟ, zero(T), q, base) === ḟ
        @test ḟ == allocating_rhs(q, params.arakawa, params.v, params.G)
    end

    @testset "a Dual state takes the buffers of its own cache" begin
        # the Jacobian of `ImplicitMidpoint`'s Newton solver passes `ForwardDiff.Dual` states,
        # which the `Float64` buffers cannot hold
        model, prob = seeded_problem(Float64)
        params = prob.parameters
        nv = size(model.distribution.values, 2)
        q = Vector(prob.ics.q)
        qd = ForwardDiff.Dual.(q, 1.0)
        ḟ, ḟd = similar(q), similar(qd)

        @test which(VM._grid_1d1v_rhs!, typeof((ḟ, 0.0, q, params))) !=
              which(VM._grid_1d1v_rhs!, typeof((ḟd, 0.0, qd, params)))
        @test !haskey(params.cache.caches, eltype(qd))
        @test @inferred(VM._grid_1d1v_rhs!(ḟd, 0.0, qd, params)) === ḟd
        # the first call of the type builds the cache of that type, in that type
        cache = params.cache[eltype(qd)]
        @test cache.ρ isa Vector{eltype(qd)} && length(cache.ρ) == length(params.ρ)
        @test cache.ϕ isa Vector{eltype(qd)} && length(cache.ϕ) == length(params.ϕ)
        @test cache.h isa Vector{eltype(qd)} && length(cache.h) == length(params.h)
        # the state's own type holds the buffers of `params`
        @test params.cache[Float64].ρ === params.ρ

        # `G ρ` takes BLAS in `Float64` and the generic product for `Dual`, and the allocating
        # form sums `ρ` in another order, so the parts agree to the bound of the exact testset,
        # not bit for bit
        VM._grid_1d1v_rhs!(ḟ, 0.0, q, params)
        @test maximum(abs, ForwardDiff.value.(ḟd) .- ḟ) ≤
              8nv * eps(Float64) * maximum(abs, ḟ)
        ref = ForwardDiff.partials.(allocating_rhs(qd, params.arakawa, params.v, params.G), 1)
        @test maximum(abs, ForwardDiff.partials.(ḟd, 1) .- ref) ≤
              8nv * eps(Float64) * maximum(abs, ref)
    end

    @testset "a Dual state allocates nothing, $T, $nx × $nv" for T in (Float32, Float64),
        (nx, nv) in ((16, 64), (64, 256))

        _, prob = seeded_problem(T; nx, nv)
        params = prob.parameters
        qd = dual12(Vector(prob.ics.q))
        @test eltype(qd) === ForwardDiff.Dual{Nothing, T, 12}
        ḟd = similar(qd)
        if Base.JLOptions().code_coverage == 0
            @test rhs_allocations(VM._grid_1d1v_rhs!, ḟd, zero(T), qd, params) == 0
        end
        # the control: `params` without the cache take the allocating method
        base = (; params.arakawa, params.v, params.G)
        @test rhs_allocations(VM._grid_1d1v_rhs!, similar(qd), zero(T), qd, base) > 0
    end

    @testset "the Jacobian of the right-hand side, $T" for T in (Float32, Float64)
        _, prob = seeded_problem(T)
        params = prob.parameters
        q = Vector(prob.ics.q)
        J = ForwardDiff.jacobian((ḟ, f) -> VM._grid_1d1v_rhs!(ḟ, zero(T), f, params),
            similar(q), q)
        Jref = reference_jacobian(q, params.arakawa, params.v, params.G)
        @test eltype(J) === T
        @test norm(J - Jref) ≤ 8 * eps(T) * norm(Jref)
    end

    @testset "the right-hand side allocates nothing, $T" for T in (Float32, Float64)
        _, prob = seeded_problem(T)
        params = prob.parameters
        q = Vector(prob.ics.q)
        ḟ = similar(q)
        # the problem's own state is a wrapper, not a `Vector`, as the integrators' states are
        qs = copy(prob.ics.q)
        @test !(qs isa Vector)

        # coverage instrumentation allocates on every counted line, so under it the assertion
        # does not run rather than being weakened
        qd = ForwardDiff.Dual.(q, one(T))
        if Base.JLOptions().code_coverage == 0
            @test rhs_allocations(VM._grid_1d1v_rhs!, ḟ, zero(T), q, params) == 0
            @test rhs_allocations(VM._grid_1d1v_rhs!, similar(qs), zero(T), qs, params) == 0
            @test rhs_allocations(VM._grid_1d1v_rhs!, similar(qd), zero(T), qd, params) == 0
        end
        # the control: the barrier sees the allocations of the method of `params` without
        # buffers
        base = (; params.arakawa, params.v, params.G)
        @test rhs_allocations(VM._grid_1d1v_rhs!, similar(q), zero(T), q, base) > 0
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

    @testset "the closure keyword, $T" for T in (Float32, Float64)
        # the seeded problem, applied to a state whose end columns are standard normal, where
        # the two closures differ
        nx, nv = 16, 64
        dist = maxwellian_grid(T, nx, nv, 7.5, 0.5)
        model = VlasovPoisson(dist, Potential(FFTWBasis((zero(T), T(4π)), nx)))
        tspan, dt = (zero(T), one(T)), T(0.1)
        default = GeometricIntegrator(model, tspan, dt).equation
        periodic = GeometricIntegrator(model, tspan, dt;
            closure = PeriodicVelocityClosure()).equation
        bounded = GeometricIntegrator(model, tspan, dt;
            closure = BoundedVelocityClosure()).equation

        f = copy(dist.values)
        f[:, [1, nv]] .= draw(Xoshiro(11), T, nx, 2)
        ḟ = applied(default, f)
        @test eltype(ḟ) === T
        @test ḟ == applied(periodic, f)
        @test default.equation.v === periodic.equation.v === VM._grid_1d1v_rhs!
        @test keys(bounded.parameters) == keys(default.parameters)

        ref = bounded_field(bounded.parameters, f)
        ḟb = applied(bounded, f)
        @test eltype(ḟb) === T
        @test maximum(abs, ḟb .- ref) ≤ 8nv * eps(T) * maximum(abs, ref)
        @test maximum(abs, ḟ .- ref) > 8nv * eps(T) * maximum(abs, ref)

        @test_throws r"^TypeError: in keyword argument closure, expected VlasovMethods\.VelocityClosure, got a value of type Symbol$" GeometricIntegrator(
            model, tspan, dt; closure = :bounded)
        @test_throws r"^TypeError: in GeometricIntegrator, in keyword closure, expected Union\{.*\}, got a value of type .*OtherClosure$" GeometricIntegrator(
            model, tspan, dt; closure = OtherClosure())
    end

    @testset "the bounded field is the T̄ field, $T" for T in (Float32, Float64)
        # random `f`, its end columns drawn too, reaches the end nodes, where the closure
        # differs from the periodic one
        for (nx, nv) in ((8, 12), (16, 24)), seed in 1:3

            _, prob = bounded_problem(T, nx, nv)
            p = prob.parameters
            f = draw(Xoshiro(seed), T, nx, nv)
            ḟ = applied(prob, f)
            @test eltype(ḟ) === T
            @test rel(ḟ, tbar_field(p, f)) ≤ tolerance(T)
            @test all(==(0), ḟ[:, 1]) && all(==(0), ḟ[:, nv])
        end
    end

    @testset "the end columns are not read, $T" for T in (Float32, Float64)
        nx, nv = 16, 24
        _, prob = bounded_problem(T, nx, nv)
        f = draw(Xoshiro(4), T, nx, nv)
        f[:, [1, nv]] .= 0
        normal = copy(f)
        normal[:, [1, nv]] .= draw(Xoshiro(5), T, nx, 2)
        nan, inf = copy(f), copy(f)
        nan[3, 1] = nan[6, nv] = T(NaN)
        inf[3, 1] = inf[6, nv] = T(Inf)

        ḟ = applied(prob, f)
        @test all(isfinite, ḟ)
        for g in (normal, nan, inf)
            ḟg = applied(prob, g)
            @test ḟg == ḟ
            @test all(isfinite, ḟg)
        end
    end

    @testset "energy, reflection and mass flux, $T" for T in (Float32, Float64)
        nx, nv = 16, 32
        _, prob = bounded_problem(T, nx, nv)
        p = prob.parameters
        tol = tolerance(T)
        reflect(a) = a[:, end:-1:1]
        @test collect(p.v) ≈ -reverse(collect(p.v)) atol = 4 * eps(T) * 6
        for seed in 1:5
            f = cut(draw(Xoshiro(seed), T, nx, nv))
            ḟ = applied(prob, f)
            Ph = cut_h(p, f)
            scale(u) = sqrt(ip(p, u, u) * ip(p, ḟ, ḟ))

            @test abs(ip(p, Ph, ḟ)) / scale(Ph) ≤ tol
            @test rel(applied(prob, reflect(f)), -reflect(ḟ)) ≤ tol
            # the mass rate is the boundary flux, and need not vanish
            one_ = one.(f)
            @test abs(ip(p, one_, ḟ) - boundary_flux(p, f)) / scale(one_) ≤ tol

            # the control: the uncut `h` in the bracket does not keep the energy
            ḟu = cut(bracket(p, uncut_h(p, f), f))
            @test abs(ip(p, Ph, ḟu)) / sqrt(ip(p, Ph, Ph) * ip(p, ḟu, ḟu)) ≥ 1e-3
        end
    end

    @testset "the edges of the bounded right-hand side, $T" for T in (Float32, Float64)
        nx, nv = 16, 24
        _, prob = bounded_problem(T, nx, nv)
        params = prob.parameters
        rhs! = prob.equation.v
        q = vec(draw(Xoshiro(6), T, nx, nv))
        sentinel = T(12345)

        # a state or a derivative of another length: the message of the periodic closure, and
        # the output untouched
        msg(nf, nḟ) = "_grid_1d1v_rhs! on a 16 × 24 grid takes a state f and a derivative ḟ " *
                      "of nx * nv values each; got $nf and $nḟ"
        ḟ = fill(sentinel, nx * nv)
        @test_throws DimensionMismatch(msg(nx * nv + 1, nx * nv)) rhs!(
            ḟ, zero(T), [q; one(T)], params)
        @test all(==(sentinel), ḟ)
        long = fill(sentinel, nx * nv + 1)
        @test_throws DimensionMismatch(msg(nx * nv, nx * nv + 1)) rhs!(
            long, zero(T), q, params)
        @test all(==(sentinel), long)
        # a derivative that is the state: the message of the periodic closure's bracket, and the
        # state untouched
        qa = copy(q)
        @test_throws ArgumentError("Pc shares memory with û or c, which the stencil reads " *
                                   "while it writes Pc") rhs!(qa, zero(T), qa, params)
        @test qa == q

        # a non-finite value of an unknown propagates
        for bad in (T(NaN), T(Inf))
            fb = reshape(copy(q), nx, nv)
            fb[3, 5] = bad
            @test any(!isfinite, applied(prob, fb))
        end

        # nv = 3: the one unknown row has no velocity neighbour among the unknowns
        _, prob3 = bounded_problem(T, 5, 3)
        @test all(==(0), applied(prob3, draw(Xoshiro(7), T, 5, 3)))

        # a Dual state takes the buffers of its own cache
        qd = dual12(q)
        ḟd = similar(qd)
        @test rhs!(ḟd, zero(T), qd, params) === ḟd
        ḟq = vec(applied(prob, reshape(q, nx, nv)))
        @test maximum(abs, ForwardDiff.value.(ḟd) .- ḟq) ≤ 8nv * eps(T) * maximum(abs, ḟq)

        # the allocations of the state's own element type, on the seeded 16 × 64 fixture;
        # skipped where a run forces `--check-bounds=yes`, and under coverage instrumentation
        model64, _ = seeded_problem(T)
        prob64 = GeometricIntegrator(model64, (zero(T), one(T)), T(0.1);
            closure = BoundedVelocityClosure()).equation
        q64 = Vector(prob64.ics.q)
        if Base.JLOptions().check_bounds != 1 && Base.JLOptions().code_coverage == 0
            @test rhs_allocations(prob64.equation.v, similar(q64), zero(T), q64,
                prob64.parameters) == 0
        end
        # the control: the barrier sees the allocations of the allocating periodic method
        base = (; prob64.parameters.arakawa, prob64.parameters.v, prob64.parameters.G)
        @test rhs_allocations(VM._grid_1d1v_rhs!, similar(q64), zero(T), q64, base) > 0
    end

    @testset "the initial state of a bounded run, $T" for T in (Float32, Float64)
        nx, nv = 16, 64
        model, _ = seeded_problem(T)
        @test all(!=(0), model.distribution.values[:, [1, nv]])
        before = copy(model.distribution.values)

        dt = T(0.1)
        int = GeometricIntegrator(model, (zero(T), 3dt), dt, GeometricIntegrators.RK4();
            closure = BoundedVelocityClosure())
        q₀ = reshape(Vector(int.equation.ics.q), nx, nv)
        @test all(==(0), q₀[:, [1, nv]])
        @test q₀[:, 2:(nv - 1)] == before[:, 2:(nv - 1)]
        @test model.distribution.values == before

        h5file = joinpath(mktempdir(), "grid_vlasov_poisson_bounded.hdf5")
        @test run!(int, h5file) === model.distribution
        @test all(==(0), model.distribution.values[:, [1, nv]])
        @test model.distribution.values[:, 2:(nv - 1)] != before[:, 2:(nv - 1)]
    end
end
