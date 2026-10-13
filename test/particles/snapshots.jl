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
    for f in (:X, :V, :A, :Φ, :W, :K, :M)
        getfield(s, f) .= randn(rng, T, size(getfield(s, f)))
    end
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

"The snapshot matrix of `Y` by hand: column `(p - 1) * NS + t` is `vec(Y[:, :, t, p])`."
columns(Y) = reduce(hcat, [vec(Y[:, :, t, p]) for p in 1:NPARAM for t in 1:NS])

@testset "the snapshot matrix is a column-per-snapshot reshape" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)
        M = VlasovMethods._snapshot_matrix(s.X)
        @test M == columns(s.X)

        # it shares the memory of the snapshots and copies nothing
        M[1] += one(T)
        @test s.X[1, 1, 1, 1] == M[1]

        # `Φ` has no component axis, so it has no snapshot matrix of this layout
        @test_throws MethodError VlasovMethods._snapshot_matrix(s.Φ)
    end
end

@testset "`nd` comes from the particle list, `nh` from the Poisson solver" begin
    for T in (Float32, Float64)
        rng = Xoshiro(20261012)
        poisson = PoissonSolverFFT(FFTWBasis((zero(T), one(T)), NH))
        ps = parameter_space(T)
        pl = ParticleList(randn(rng, T, ND, NP), randn(rng, T, ND, NP), ones(T, 1, NP))
        s = Snapshots(pl, poisson, NS, ps)

        @test size(s.X) == size(s.V) == size(s.A) == (ND, NP, NS, NPARAM)
        @test size(s.Φ) == (NH, NS, NPARAM)
        @test size(s.W) == size(s.K) == size(s.M) == (NS, NPARAM)
        @test eltype(s.X) === T
        @test s.paramspace == ps

        # `q` and `p` of one size: a list whose `v` holds another `nd` is refused
        pl₁ = ParticleList(randn(rng, T, 1, NP), randn(rng, T, ND, NP), ones(T, 1, NP))
        @test_throws ArgumentError Snapshots(pl₁, poisson, NS, ps)
    end
end

@testset "the inner constructor refuses a mismatched array" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)

        args = (s.X, s.V, s.A, s.Φ, s.W, s.K, s.M, s.paramspace)

        # one argument alone replaced: an array off by one along one axis, or a `ParameterSpace`
        # of five samples
        Y′ = zeros(T, ND, NP + 1, NS, NPARAM)
        Φ′ = zeros(T, NH, NS + 1, NPARAM)
        D′ = zeros(T, NS + 1, NPARAM)
        ps′ = ParameterSpace(Parameter(:β, T(5.95), T(6.05), NPARAM - 1))
        for (i, a) in (2 => Y′, 3 => Y′, 4 => Φ′, 5 => D′, 6 => D′, 7 => D′, 8 => ps′)
            @test_throws ArgumentError Snapshots(Base.setindex(args, a, i)...)
        end

        # the fields are `Array`s, and another array type is not converted in silence
        X′ = view(s.X, :, :, :, :)
        @test_throws MethodError Snapshots(Base.setindex(args, X′, 1)...)
    end
end

@testset "`==`, `isequal` and `hash` over the eight fields" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)
        s′ = deepcopy(s)
        @test s == s′
        @test isequal(s, s′)
        @test hash(s) == hash(s′)
        @test length(unique([s, s′])) == 1
        for field in (:X, :V, :A, :Φ, :W, :K, :M, :paramspace)
            @test perturbed(s, field) != s
        end

        # a signed zero is `==` and not `isequal`, as it is in the arrays
        s′.X[1] = zero(T)
        s₋ = deepcopy(s′)
        s₋.X[1] = -zero(T)
        @test s₋ == s′
        @test !isequal(s₋, s′)
    end
end

@testset "the particle basis reduces `[X; V]`, the field basis `A`" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)
        ps = s.paramspace
        pars = (μ = one(T),)
        QP, F = [columns(s.X); columns(s.V)], columns(s.A)

        for alg in (PODSVD(), PODEVD(), CotangentLiftSVD(), CotangentLiftEVD(), ComplexSVD())
            b = particle_basis(alg, s; parameters = pars, k = 5)
            @test Matrix(b) == Matrix(ReducedBasis(alg, QP, pars, ps; k = 5))
            @test b.parameters == pars
            @test b.paramspace == ps
            @test particle_basis(alg, s; tolerance = 1e-3).k ==
                  ReducedBasis(alg, QP, NamedTuple(), ps; tolerance = 1e-3).k
        end

        for alg in (PODSVD(), PODEVD())
            b = field_basis(alg, s; parameters = pars, k = 5)
            @test Matrix(b) == Matrix(ReducedBasis(alg, F, pars, ps; k = 5))
            @test b.parameters == pars
            @test b.paramspace == ps

            # at full rank each POD tag reproduces `A`, whose columns span another subspace than
            # those of `X`
            full = field_basis(alg, s; k = NS * NPARAM)
            @test full.parameters == NamedTuple()
            @test projection_error(full, F) ≤ √eps(T) * norm(F)
        end

        # the rows of `A` are not a `[Q; P]` stacking, so no PSD tag splits them
        @test_throws MethodError field_basis(CotangentLiftSVD(), s; k = 5)
    end
end

@testset "the HDF5 round trip" begin
    for T in (Float32, Float64)
        s = snapshot_fixture(T)
        s₀ = Snapshots(T, ND, NP, 0, NS, s.paramspace)
        @test size(s₀.Φ) == (0, NS, NPARAM)

        # a missing or swapped key is only visible if the arrays differ from each other
        @test allunique([s.X, s.V, s.A])
        @test allunique([s.W, s.K, s.M])

        mktempdir() do dir
            fpath = joinpath(dir, "snapshots.h5")
            h5open(fpath, "w") do file
                h5save(file, s; path = "snapshots")
            end
            s₂ = h5open(fpath, "r") do file
                h5load(Snapshots, file; path = "snapshots")
            end
            @test s₂ isa Snapshots{T}
            @test s₂ == s

            # the file-path forms of `ReducedComplexityModeling`, on a model with no potential
            fpath₀ = joinpath(dir, "no_potential.h5")
            h5save(fpath₀, s₀; path = "snapshots")
            @test h5load(Snapshots, fpath₀; path = "snapshots") == s₀
        end
    end
end
