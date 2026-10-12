using VlasovMethods
using HDF5
using LinearAlgebra
using ParticleMethods
using PoissonSolvers
using Random
using ReducedBasisMethods: ComplexSVD, CotangentLiftEVD, CotangentLiftSVD, PODSVD, PODEVD,
                           ReducedBasis, projection_error
using ReducedComplexityModeling: Parameter, ParameterSpace, h5load, h5save
using Test

# The snapshot sweep: two position and two velocity components, one Poisson solver of `NH` degrees
# of freedom, `NS` saved times and the six training `β`. No two extents are equal, so a swapped
# axis changes a size.
const ND, NP, NH, NS, NPARAM = 2, 100, 16, 11, 6

"The six training `β` of the sweep, on `[5.95, 6.05]`, as bounds of the element type `T`."
parameter_space(T) = ParameterSpace(Parameter(:β, T(5.95), T(6.05), NPARAM))

"""
A `Snapshots{T}` of the size above, with every entry of the seven arrays drawn from a seeded
`Xoshiro`. No entry keeps the constructor's zero, and no two of the seven arrays are equal.
"""
function snapshot_fixture(T; seed = 20261012)
    s = Snapshots(T, ND, NP, NH, NS, parameter_space(T))
    rng = Xoshiro(seed)
    s.X .= randn(rng, T, size(s.X))
    s.V .= randn(rng, T, size(s.V))
    s.A .= randn(rng, T, size(s.A))
    s.Φ .= randn(rng, T, size(s.Φ))
    s.W .= randn(rng, T, size(s.W))
    s.K .= randn(rng, T, size(s.K))
    s.M .= randn(rng, T, size(s.M))
    return s
end

"`s` with one entry of `field` raised by one, or another `ParameterSpace`, and nothing else."
function perturbed(s::Snapshots{T}, field::Symbol) where {T}
    args = (copy(s.X), copy(s.V), copy(s.A), copy(s.Φ), copy(s.W), copy(s.K), copy(s.M))
    field === :paramspace &&
        return Snapshots(args..., ParameterSpace(Parameter(:β, T(5.9), T(6.1), NPARAM)))
    a = args[findfirst(isequal(field), (:X, :V, :A, :Φ, :W, :K, :M))]
    a[1] += one(T)
    return Snapshots(args..., s.paramspace)
end

@testset "the snapshot layout" begin
    @testset "check 1 — the `nd == 2` sweep survives" begin
        for T in (Float32, Float64)
            s = snapshot_fixture(T)
            M = VlasovMethods._snapshot_matrix(s.X)

            @test size(M) == (ND * NP, NS * NPARAM)
            for p in 1:NPARAM
                for t in 1:NS
                    c = (p - 1) * NS + t
                    @test M[:, c] == vec(s.X[:, :, t, p])
                end
            end

            # a `reshape` of the array: it shares the memory of the snapshots and copies nothing
            M[1] += one(T)
            @test s.X[1, 1, 1, 1] == M[1]
        end
    end

    @testset "the particle matrix is the stacked `[Q; P]`" begin
        for T in (Float32, Float64)
            s = snapshot_fixture(T)
            M = VlasovMethods._particle_snapshots(s)
            F = VlasovMethods._field_snapshots(s)

            # `2·nd·nₚ` rows, the positions first and the velocities second, and no interleaving:
            # RBM reduces any even number of rows, so an interleaved matrix raises no error
            @test size(M) == (2 * ND * NP, NS * NPARAM)
            @test size(F) == (ND * NP, NS * NPARAM)
            for p in 1:NPARAM
                for t in 1:NS
                    c = (p - 1) * NS + t
                    @test M[:, c] == vcat(vec(s.X[:, :, t, p]), vec(s.V[:, :, t, p]))
                    @test F[:, c] == vec(s.A[:, :, t, p])
                end
            end
        end
    end

    @testset "check 2 — no column mixes two times or two parameters" begin
        for T in (Float32, Float64)
            s = snapshot_fixture(T)
            s.X[ND, :, :, :] .= zero(T)

            σ₂ = svdvals(VlasovMethods._snapshot_matrix(s.X))
            # the `nd = 1` reference is assembled by hand, column `(p - 1) * NS + t` =
            # `vec(X[1, :, t, p])`, and not with `_snapshot_matrix`: a layout that does not depend
            # on `nd` would be applied to both sides and cancel from the comparison
            M₁ = reduce(hcat, [s.X[1, :, t, p] for p in 1:NPARAM for t in 1:NS])
            σ₁ = svdvals(M₁)

            # a swapped row keeps the singular values; a column that mixes two times does not
            @test length(σ₂) == length(σ₁) == NS * NPARAM
            @test isapprox(σ₂, σ₁; atol = √eps(T) * first(σ₁), rtol = 0)
        end
    end

    @testset "check 3 — `nd` comes from the list, `nh` from the solver" begin
        for T in (Float32, Float64)
            rng = Xoshiro(20261012)
            poisson = PoissonSolverFFT(FFTWBasis((zero(T), one(T)), NH))
            ps = parameter_space(T)
            pl = ParticleList(randn(rng, T, ND, NP), randn(rng, T, ND, NP), ones(T, 1, NP))
            s = Snapshots(pl, poisson, NS, ps)

            @test size(s.X) == (ND, NP, NS, NPARAM)
            @test size(s.V) == size(s.X)
            @test size(s.A) == size(s.X)
            @test size(s.Φ) == (NH, NS, NPARAM)
            @test size(s.W) == (NS, NPARAM)
            @test size(s.K) == (NS, NPARAM)
            @test size(s.M) == (NS, NPARAM)
            @test eltype(s.X) === T
            @test s.paramspace == ps

            # `q` and `p` of one size: a list whose `v` holds another `nd` is refused
            pl₁ = ParticleList(randn(rng, T, 1, NP), randn(rng, T, ND, NP), ones(T, 1, NP))
            @test_throws ArgumentError Snapshots(pl₁, poisson, NS, ps)

            # the array form, and a model whose solver carries no potential
            s₀ = Snapshots(T, ND, NP, 0, NS, ps)
            @test size(s₀.X) == (ND, NP, NS, NPARAM)
            @test size(s₀.Φ) == (0, NS, NPARAM)
        end
    end

    @testset "check 4 — `A` feeds the field basis, and only it" begin
        for T in (Float32, Float64)
            s = snapshot_fixture(T)
            s′ = Snapshots(copy(s.X), copy(s.V), s.A .+ one(T), copy(s.Φ), copy(s.W),
                copy(s.K), copy(s.M), s.paramspace)

            # the particle basis is built from `X` and `V` only
            for alg in (PODSVD(), PODEVD(), CotangentLiftSVD(), CotangentLiftEVD(), ComplexSVD())
                @test Matrix(particle_basis(alg, s′; k = 5)) ==
                      Matrix(particle_basis(alg, s; k = 5))
            end

            # the field basis has the per-snapshot length of `A`, not the particle count of `X`
            fb = field_basis(PODSVD(), s; k = 5)
            @test size(Matrix(fb), 1) == length(s.A[:, :, 1, 1]) == ND * NP

            # at full rank each POD tag reproduces `A`, whose columns span another subspace than
            # `X`'s, so a field basis built from `X` does not pass
            F = VlasovMethods._field_snapshots(s)
            for alg in (PODSVD(), PODEVD())
                full = field_basis(alg, s; k = NS * NPARAM)
                @test projection_error(full, F) ≤ √eps(T) * norm(F)
            end
        end
    end

    @testset "check 5 — the HDF5 round trip" begin
        for T in (Float32, Float64)
            s = snapshot_fixture(T)

            # a missing or swapped key is only visible if the arrays differ from each other
            @test length(unique([s.X, s.V, s.A])) == 3
            @test length(unique([s.W, s.K, s.M])) == 3

            mktempdir() do dir
                fpath = joinpath(dir, "snapshots.h5")
                h5open(fpath, "w") do file
                    h5save(file, s; path = "snapshots")
                end

                s₂ = h5open(fpath, "r") do file
                    h5load(Snapshots, file; path = "snapshots")
                end

                @test s₂ isa Snapshots{T}
                @test eltype(s₂.X) === T
                @test size(s₂.X) == (ND, NP, NS, NPARAM)
                @test size(s₂.Φ) == (NH, NS, NPARAM)
                @test size(s₂.W) == (NS, NPARAM)
                @test s₂ == s
            end
        end
    end
end

@testset "the inner constructor refuses a mismatched array" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)

        # one array alone, off by one along one axis, and one `ParameterSpace` of five samples
        V′ = zeros(T, ND, NP + 1, NS, NPARAM)
        A′ = zeros(T, ND, NP + 1, NS, NPARAM)
        Φ′ = zeros(T, NH, NS + 1, NPARAM)
        W′ = zeros(T, NS + 1, NPARAM)
        K′ = zeros(T, NS + 1, NPARAM)
        M′ = zeros(T, NS + 1, NPARAM)
        ps′ = ParameterSpace(Parameter(:β, T(5.95), T(6.05), NPARAM - 1))

        @test_throws ArgumentError Snapshots(s.X, V′, s.A, s.Φ, s.W, s.K, s.M, s.paramspace)
        @test_throws ArgumentError Snapshots(s.X, s.V, A′, s.Φ, s.W, s.K, s.M, s.paramspace)
        @test_throws ArgumentError Snapshots(s.X, s.V, s.A, Φ′, s.W, s.K, s.M, s.paramspace)
        @test_throws ArgumentError Snapshots(s.X, s.V, s.A, s.Φ, W′, s.K, s.M, s.paramspace)
        @test_throws ArgumentError Snapshots(s.X, s.V, s.A, s.Φ, s.W, K′, s.M, s.paramspace)
        @test_throws ArgumentError Snapshots(s.X, s.V, s.A, s.Φ, s.W, s.K, M′, s.paramspace)
        @test_throws ArgumentError Snapshots(s.X, s.V, s.A, s.Φ, s.W, s.K, s.M, ps′)
    end
end

@testset "`==` over the eight fields" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)
        @test s == deepcopy(s)
        for field in (:X, :V, :A, :Φ, :W, :K, :M, :paramspace)
            @test perturbed(s, field) != s
        end
    end
end

@testset "the two builders" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)
        pars = (μ = one(T),)

        for alg in (PODSVD(), PODEVD(), CotangentLiftSVD(), CotangentLiftEVD(), ComplexSVD())
            b = particle_basis(alg, s; k = 5)
            @test b.k == 5
            @test b.parameters == NamedTuple()
            @test b.paramspace == s.paramspace

            # a `parameters` keyword reaches the returned basis, so the builder passes it on
            @test particle_basis(alg, s; parameters = pars, k = 5).parameters == pars

            # `k`, `tolerance` and the parameter space reach `ReducedBasis` unchanged
            @test particle_basis(alg, s; tolerance = 1e-3).k ==
                  ReducedBasis(alg, VlasovMethods._particle_snapshots(s), NamedTuple(),
                s.paramspace; tolerance = 1e-3).k
        end

        for alg in (PODSVD(), PODEVD())
            b = field_basis(alg, s; k = 5)
            @test b.k == 5
            @test b.parameters == NamedTuple()
            @test b.paramspace == s.paramspace

            @test field_basis(alg, s; parameters = pars, k = 5).parameters == pars

            @test field_basis(alg, s; tolerance = 1e-3).k ==
                  ReducedBasis(alg, VlasovMethods._field_snapshots(s), NamedTuple(),
                s.paramspace; tolerance = 1e-3).k
        end

        # the field basis is built from the `nd·nₚ` rows of `A`, which no PSD tag splits
        @test_throws MethodError field_basis(CotangentLiftSVD(), s; k = 5)
    end
end
