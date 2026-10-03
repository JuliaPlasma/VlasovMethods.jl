using VlasovMethods
using LinearAlgebra
using Test

const multiindex = VlasovMethods.multiindex
const linearindex = VlasovMethods.linearindex

# --- helpers for the invariant test (check 5) ---------------------------------

# the cubic contraction Q₃[o] = Σ MC̃_cubic[o, a, b, c] f[a] f[b] f[c]
function _contract_cubic!(Q, MC, f)
    N = length(f)
    for o in 1:N
        s = zero(eltype(Q))
        for a in 1:N, b in 1:N, c in 1:N
            s += MC[o, a, b, c] * f[a] * f[b] * f[c]
        end
        Q[o] = s
    end
    return Q
end

# the quadratic contraction Q₂[o] = Σ MC̃_quadratic[o, a, c] f[a] f[c]
function _contract_quadratic!(Q, MC, f)
    N = length(f)
    for o in 1:N
        s = zero(eltype(Q))
        for a in 1:N, c in 1:N

            s += MC[o, a, c] * f[a] * f[c]
        end
        Q[o] = s
    end
    return Q
end

# ten RK4 steps of df/dt = Q(f), returning the grid
function _rk4(f_grid, MC, nx, nv, Δt, n_steps, cubic::Bool)
    T = eltype(f_grid)
    N = nx * nv
    f = vec(copy(f_grid))
    Q1 = zeros(T, N)
    Q2 = zeros(T, N)
    Q3 = zeros(T, N)
    Q4 = zeros(T, N)
    for _ in 1:n_steps
        cubic ? _contract_cubic!(Q1, MC, f) : _contract_quadratic!(Q1, MC, f)
        c2 = f .+ (Δt / 2) .* Q1
        cubic ? _contract_cubic!(Q2, MC, c2) : _contract_quadratic!(Q2, MC, c2)
        c3 = f .+ (Δt / 2) .* Q2
        cubic ? _contract_cubic!(Q3, MC, c3) : _contract_quadratic!(Q3, MC, c3)
        c4 = f .+ Δt .* Q3
        cubic ? _contract_cubic!(Q4, MC, c4) : _contract_quadratic!(Q4, MC, c4)
        f = f .+ (Δt / 6) .* (Q1 .+ 2 .* Q2 .+ 2 .* Q3 .+ Q4)
    end
    return reshape(f, nx, nv)
end

# the three moments ∫dv·f, ∫vdv·f, ∫v²dv·f over the x nodes
_moments(∫dv, ∫vdv, ∫v²dv, f) = (∫dv * f, ∫vdv * f, ∫v²dv * f)

# the matrices that go into the assembler, by the rectangle rule of _apply_∫dv!
function _moment_matrices(::Type{T}, nx, nv, v, h₂) where {T}
    li = LinearIndices((nx, nv))
    N = nx * nv
    ∫dv = zeros(T, nx, N)
    ∫vdv = zeros(T, nx, N)
    ∫v²dv = zeros(T, nx, N)
    for i in 1:nx, j in 1:nv

        k = li[i, j]
        ∫dv[i, k] = h₂
        ∫vdv[i, k] = h₂ * v[j]
        ∫v²dv[i, k] = h₂ * v[j]^2
    end
    return ∫dv, ∫vdv, ∫v²dv
end

function _shifted_maxwellian(::Type{T}, nx, nv, v, A) where {T}
    return [A * (1 + T(0.2) * i) / sqrt(T(2π) * (T(0.5) + T(0.1) * i)) *
            exp(-(v[j] - (T(0.5) + T(0.2) * i))^2 / (2 * (T(0.5) + T(0.1) * i)))
            for i in 1:nx, j in 1:nv]
end

_relative_drift(a, b) = norm(b .- a) / norm(a)

# Function barriers for the allocation assertions of check 4: the argument is concrete, and the
# call runs once before `@allocated` measures it. For the function argument `f::F`, the `where`
# clause forces specialization: without it the call would dispatch dynamically and allocate on
# its own.
_getindex3(x, i, j, k) = x[i, j, k]
function _alloc_getindex(x::X, i::Int, j::Int, k::Int) where {X}
    (_getindex3(x, i, j, k); @allocated _getindex3(x, i, j, k))
end
function _alloc_call3(f::F, a::A, b::B, c::C) where {F, A, B, C}
    (f(a, b, c); @allocated f(a, b, c))
end
_alloc_call1(f::F, a::A) where {F, A} = (f(a); @allocated f(a))

@testset "Collisions" begin
    # every test below runs in Float32 and Float64 with nx ≠ nv and nv ≥ 3
    for T in (Float32, Float64)
        @testset "check 1: CollisionTensor indexing and the hand projection" begin
            nx, nv = 3, 5
            v = collect(range(T(-3), T(3); length = nv))
            hx = one(T) / nx
            hv = (v[end] - v[1]) / (nv - 1)
            qc = VlasovMethods.QuadraticCollisions(nx, nv, hx, hv, v)
            ct = VlasovMethods.CollisionTensor(T, nx, nv, qc)

            # the three-Int getindex returns ct.f at the same three multi-indices, for every
            # entry of the grid: a permutation or a wrong index of any of the three fails
            N = nx * nv
            @test all(ct[i, j, k] == qc(multiindex(i, nx, nv), multiindex(j, nx, nv),
                          multiindex(k, nx, nv)) for i in 1:N, j in 1:N, k in 1:N)

            # the hand projection with three different projection matrices
            # (the one-argument constructor would set all three equal and hide a swap)
            Pi = T[sin(T(0.3) * (a + b)) for a in 1:N, b in 1:2]
            Pj = T[cos(T(0.7) * (a + 2b)) for a in 1:N, b in 1:3]
            Pk = T[exp(T(-0.1) * a * b) for a in 1:N, b in 1:2]
            rt = VlasovMethods.ReducedCollisionTensor(ct, Pi, Pj, Pk)
            @test size(rt) == (2, 3, 2)
            for i in axes(rt, 1), j in axes(rt, 2), k in axes(rt, 3)
                hand = zero(T)
                for m1 in 1:nx, m2 in 1:nv

                    m2₋ = mod1(m2 - 1, nv)
                    m2₊ = mod1(m2 + 1, nv)
                    for o2 in (m2₋, m2, m2₊), n2 in 1:nv

                        M = CartesianIndex(m1, m2)
                        Nn = CartesianIndex(m1, n2)
                        O = CartesianIndex(m1, o2)
                        m = linearindex(M, nx, nv)
                        n = linearindex(Nn, nx, nv)
                        o = linearindex(O, nx, nv)
                        hand += qc(M, Nn, O) * Pi[m, i] * Pj[n, j] * Pk[o, k]
                    end
                end
                @test rt[i, j, k] ≈ hand rtol = sqrt(eps(T))
            end
        end

        @testset "check: the stencil matches the assembler on the same grid" begin
            nx, nv = 3, 5
            v = collect(range(T(-3), T(3); length = nv))
            hv = (v[end] - v[1]) / (nv - 1)
            hx = one(T) / nx
            qc = VlasovMethods.QuadraticCollisions(nx, nv, hx, hv, v)
            ct = VlasovMethods.CollisionTensor(T, nx, nv, qc)
            N = nx * nv
            ci = CartesianIndices((nx, nv))
            li = LinearIndices((nx, nv))
            ∫dv, ∫vdv, ∫v²dv = _moment_matrices(T, nx, nv, v, hv)
            MC2 = VlasovMethods._get_MC̃_quadratic(Matrix{T}(I, N, N), ∫dv, ∫vdv, ∫v²dv,
                v, ci, li, hx, hv)
            # With V = I and the rectangle-rule moments, MC̃_quadratic = qc / hv, and
            # ct[o, a, c] is qc at the same three multi-indices: the stencil and the assembler
            # agree entry by entry, through two independent implementations.
            @test all(ct[o, a, c] ≈ hv * MC2[o, a, c] for o in 1:N, a in 1:N, c in 1:N)
        end

        @testset "check 2: the fields of QuadraticCollisions" begin
            nx, nv = 3, 5
            v = collect(range(T(-3), T(3); length = nv))
            @test fieldnames(VlasovMethods.QuadraticCollisions) == (:nx, :nv, :hx, :hv, :v)
            qc = VlasovMethods.QuadraticCollisions(nx, nv, one(T) / nx,
                (v[end] - v[1]) / (nv - 1), v)
            @test qc isa VlasovMethods.QuadraticCollisions{T}
            @test qc.nx == nx && qc.nv == nv
        end

        @testset "check 3: the v argument of the two assemblers" begin
            nx, nv = 4, 5
            v = collect(range(T(-3), T(3); length = nv))
            N = nx * nv
            hx = one(T) / nx
            hv = (v[end] - v[1]) / (nv - 1)
            ci = CartesianIndices((nx, nv))
            li = LinearIndices((nx, nv))
            ∫dv, ∫vdv, ∫v²dv = _moment_matrices(T, nx, nv, v, hv)
            V = Matrix{T}(I, N, N)

            MC3 = VlasovMethods._get_MC̃_cubic(V, ∫dv, ∫vdv, ∫v²dv, v, ci, li, hx, hv)
            @test size(MC3) == (N, N, N, N)
            @test eltype(MC3) == T
            @test_throws DimensionMismatch VlasovMethods._get_MC̃_cubic(
                V, ∫dv, ∫vdv, ∫v²dv, v[1:(end - 1)], ci, li, hx, hv)
            @test_throws DimensionMismatch VlasovMethods._get_MC̃_cubic(
                V, ∫dv, ∫vdv, ∫v²dv, vcat(v, v[end]), ci, li, hx, hv)

            MC2 = VlasovMethods._get_MC̃_quadratic(V, ∫dv, ∫vdv, ∫v²dv, v, ci, li, hx, hv)
            @test size(MC2) == (N, N, N)
            @test eltype(MC2) == T
            @test_throws DimensionMismatch VlasovMethods._get_MC̃_quadratic(
                V, ∫dv, ∫vdv, ∫v²dv, v[1:(end - 1)], ci, li, hx, hv)
            @test_throws DimensionMismatch VlasovMethods._get_MC̃_quadratic(
                V, ∫dv, ∫vdv, ∫v²dv, vcat(v, v[end]), ci, li, hx, hv)
        end

        @testset "check 4: the tensor paths are concrete, inferred and allocation-free" begin
            nx, nv = 3, 5
            v = collect(range(T(-3), T(3); length = nv))
            hv = (v[end] - v[1]) / (nv - 1)
            I = CartesianIndex(1, 1)
            for vv in (v, range(T(-3), T(3); length = nv))
                qc = VlasovMethods.QuadraticCollisions(nx, nv, one(T) / nx, hv, vv)
                @test isconcretetype(fieldtype(typeof(qc), :v))
                @test @inferred(qc(I, I, I)) isa T
            end

            # the two index paths that the reduction and the collision step call
            qc = VlasovMethods.QuadraticCollisions(nx, nv, one(T) / nx, hv, v)
            ct = VlasovMethods.CollisionTensor(T, nx, nv, qc)
            N = nx * nv
            Pi = T[sin(T(0.3) * (a + b)) for a in 1:N, b in 1:2]
            Pj = T[cos(T(0.7) * (a + 2b)) for a in 1:N, b in 1:3]
            Pk = T[exp(T(-0.1) * a * b) for a in 1:N, b in 1:2]
            rt = VlasovMethods.ReducedCollisionTensor(ct, Pi, Pj, Pk)

            @test @inferred(ct[1, 2, 3]) isa T
            @test @inferred(rt[1, 2, 1]) isa T

            # each hot path asserts type stability and an empty allocation, measured
            # through a function barrier whose arguments are concrete
            @test _alloc_call3(qc, I, I, I) == 0
            @test _alloc_getindex(ct, 1, 2, 3) == 0
            @test _alloc_getindex(rt, 1, 2, 1) == 0
            @test _alloc_call1(collect, v) > 0   # the control: the barrier sees an allocation
        end

        @testset "check 5: the invariants of the two forms" begin
            nx, nv = 2, 25
            v = collect(range(T(-10), T(10); length = nv))
            hv = (v[end] - v[1]) / (nv - 1)
            hx = one(T) / nx
            ci = CartesianIndices((nx, nv))
            li = LinearIndices((nx, nv))
            N = nx * nv
            ∫dv, ∫vdv, ∫v²dv = _moment_matrices(T, nx, nv, v, hv)
            V = Matrix{T}(I, N, N)
            MC3 = VlasovMethods._get_MC̃_cubic(V, ∫dv, ∫vdv, ∫v²dv, v, ci, li, hx, hv)
            MC2 = VlasovMethods._get_MC̃_quadratic(V, ∫dv, ∫vdv, ∫v²dv, v, ci, li, hx, hv)

            n_steps = 10
            atol = 10 * n_steps * eps(T)
            for A in (one(T), T(3e-8), T(3e8))
                f_grid = _shifted_maxwellian(T, nx, nv, v, A)
                # the invariant case requires f to vanish at the ends of the v-grid
                @test maximum(abs, f_grid[:, [1, nv]]) ≤ eps(T) * maximum(abs, f_grid)

                f3 = _rk4(f_grid, MC3, nx, nv, T(1e-3) / A^2, n_steps, true)
                f2 = _rk4(f_grid, MC2, nx, nv, T(1e-3) / A, n_steps, false)
                @test maximum(abs, f3[:, [1, nv]]) ≤ eps(T) * maximum(abs, f3)
                @test maximum(abs, f2[:, [1, nv]]) ≤ eps(T) * maximum(abs, f2)

                # The invariants are asserted on the trajectory, and the zero operator is a
                # fixed point of it: both operators must move the state, or the assertions hold
                # for an operator that is identically zero.
                @test maximum(abs, f3 - f_grid) > 0
                @test maximum(abs, f2 - f_grid) > 0

                m0, m1, m2 = _moments(∫dv, ∫vdv, ∫v²dv, vec(f_grid))
                M0, M1, M2 = _moments(∫dv, ∫vdv, ∫v²dv, vec(f3))
                n0, n1, n2 = _moments(∫dv, ∫vdv, ∫v²dv, vec(f2))

                # the cubic form conserves mass, momentum and energy
                @test _relative_drift(m0, M0) ≤ atol
                @test _relative_drift(m1, M1) ≤ atol
                @test _relative_drift(m2, M2) ≤ atol

                # the quadratic form conserves mass and energy, not momentum
                @test _relative_drift(m0, n0) ≤ atol
                @test _relative_drift(m2, n2) ≤ atol
                @test _relative_drift(m1, n1) > 1e-3
            end
        end
    end

    @testset "check 5: the docstrings name the invariants" begin
        cubic = string(@doc VlasovMethods._get_MC̃_cubic)
        quadratic = string(@doc VlasovMethods._get_MC̃_quadratic)
        @test occursin("mass", cubic) && occursin("momentum", cubic) &&
              occursin("energy", cubic)
        @test occursin("mass", quadratic) && occursin("energy", quadratic)
        @test occursin("momentum", quadratic)
    end

    @testset "check 6: the docstrings state the v-end requirement" begin
        cubic = string(@doc VlasovMethods._get_MC̃_cubic)
        quadratic = string(@doc VlasovMethods._get_MC̃_quadratic)
        @test occursin("vanish", cubic) && occursin("ends", cubic)
        @test occursin("vanish", quadratic) && occursin("ends", quadratic)
    end

    @testset "edges: nv < 3 is rejected" begin
        for T in (Float32, Float64)
            v = collect(range(T(-3), T(3); length = 3))
            for nv in (1, 2)
                @test_throws ArgumentError VlasovMethods.QuadraticCollisions(
                    3, nv, one(T) / 3, (v[end] - v[1]) / 2, v[1:nv])
            end
        end
    end
end
