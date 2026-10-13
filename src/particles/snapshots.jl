"""
    Snapshots{T, PS}

The snapshots of a parameter sweep: the positions `X`, the velocities `V`, the accelerations `A`,
the potentials `Φ`, and the diagnostics `W`, `K` and `M`, on the sweep's `ParameterSpace`.

`X`, `V` and `A` are indexed `[component, particle, time, sample]` and have the size
`(nd, nₚ, nₛ, nparam)`, with `nparam = length(paramspace)`. `X[:, :, t, p]` is the `(nd, nₚ)`
matrix of the positions saved at the `t`-th of the `nₛ` saved times of the sample
`paramspace(p)`, `X[i, j, t, p]` is component `i` of particle `j`, and `V` and `A` hold the
velocities and the accelerations of that state the same way. So `vec(X[:, :, t, p])` is the
component-fastest order of `vec(particles.x)` of that state, and a snapshot matrix built from
these arrays has one column per saved time and sample.

`Φ` has the size `(nh, nₛ, nparam)`: the potential is a scalar on the one dimension of the
Poisson problem, so it carries no component axis, and `nh` is the number of degrees of freedom of
the sweep's Poisson solver. `W`, `K` and `M` have the size `(nₛ, nparam)`, one value each
per saved time and sample. `W` is the field energy, not a particle weight: the snapshots store no
weights.

The component count `nd` of `X`, `V` and `A` and the degree count `nh` of `Φ` are the leading
extents of those arrays, so they are read from the arrays rather than stored beside them.

`==` compares all eight fields, `paramspace` included.
"""
struct Snapshots{T <: Number, PS <: ParameterSpace}
    X::Array{T, 4}
    V::Array{T, 4}
    A::Array{T, 4}
    Φ::Array{T, 3}
    W::Array{T, 2}
    K::Array{T, 2}
    M::Array{T, 2}
    paramspace::PS

    function Snapshots(X::Array{T, 4}, V::Array{T, 4}, A::Array{T, 4}, Φ::Array{T, 3},
            W::Array{T, 2}, K::Array{T, 2}, M::Array{T, 2},
            paramspace::PS) where {T, PS <: ParameterSpace}
        nts = size(X)[3:4]
        size(V) == size(A) == size(X) || throw(ArgumentError(
            "X, V and A need one size, got $(size(X)), $(size(V)) and $(size(A))"))
        size(Φ)[2:3] == size(W) == size(K) == size(M) == nts || throw(ArgumentError(
            "Φ, W, K and M need the time and sample extents $(nts) of X"))
        nts[2] == length(paramspace) || throw(ArgumentError(
            "the snapshots hold $(nts[2]) samples, the parameter space holds " *
            "$(length(paramspace))"))
        new{T, PS}(X, V, A, Φ, W, K, M, paramspace)
    end
end

"""
    Snapshots(T, nd, nₚ, nh, nₛ, paramspace)

The zero `Snapshots{T}` of `nd` components, `nₚ` particles, `nh` potential degrees of freedom,
`nₛ` saved times and `length(paramspace)` samples, which a sweep fills as it runs. `nh = 0` is the
sweep of a model with no potential, such as one of `VlasovWithExternalFields`, and its empty `Φ`
round-trips through HDF5 like any other.
"""
function Snapshots(T, nd, nₚ, nh, nₛ, paramspace::ParameterSpace)
    nparam = length(paramspace)
    X = zeros(T, nd, nₚ, nₛ, nparam)
    V = zeros(T, nd, nₚ, nₛ, nparam)
    A = zeros(T, nd, nₚ, nₛ, nparam)
    Φ = zeros(T, nh, nₛ, nparam)
    W = zeros(T, nₛ, nparam)
    K = zeros(T, nₛ, nparam)
    M = zeros(T, nₛ, nparam)
    Snapshots(X, V, A, Φ, W, K, M, paramspace)
end

"""
    Snapshots(particles::ParticleList, poisson::PoissonSolver, nₛ, paramspace)

The zero `Snapshots{T}` of a sweep of `particles` on the Poisson solver `poisson`, with `nₛ` saved
times and `length(paramspace)` samples: `nd = size(particles.x, 1)` components,
`nₚ = length(particles)` particles and `nh = length(poisson)` potential degrees of freedom. Both
element types are `T`.

`particles.v` needs the `nd` rows of `particles.x`, because `X`, `V` and `A` share one component
count and the symplectic structure of the reduced model pairs positions and velocities of one
size. A list that holds another count throws an `ArgumentError`.
"""
function Snapshots(particles::ParticleList{T}, poisson::PoissonSolver{T}, nₛ,
        paramspace::ParameterSpace) where {T}
    nd = size(particles.x, 1)
    size(particles.v, 1) == nd || throw(ArgumentError(
        "the particle list has $(nd) position components and $(size(particles.v, 1)) velocity " *
        "components, which the snapshots cannot pair"))
    Snapshots(T, nd, length(particles), length(poisson), nₛ, paramspace)
end

function Base.:(==)(a::Snapshots, b::Snapshots)
    all(f -> getfield(a, f) == getfield(b, f), fieldnames(Snapshots))
end
function Base.isequal(a::Snapshots, b::Snapshots)
    all(f -> isequal(getfield(a, f), getfield(b, f)), fieldnames(Snapshots))
end
Base.hash(s::Snapshots, h::UInt) = hash(s.X, hash(s.V, hash(s.A, hash(:Snapshots, h))))

### the reduced bases of a sweep

"The `(n₁·n₂) × (n₃·n₄)` reshape of `Y`; its column `(p - 1) * n₃ + t` is `vec(Y[:, :, t, p])`."
_snapshot_matrix(Y::AbstractArray{<:Any, 4}) = reshape(Y, size(Y, 1) * size(Y, 2), :)

"""
    particle_basis(algorithm, s::Snapshots; parameters = NamedTuple(), kwargs...)

The `ReducedBasis` of the stacked `[Q; P]` snapshot matrix of `s`: the `nd * nₚ` rows of the
positions `X` above the `nd * nₚ` rows of the velocities `V`, one column per saved time and
sample. `algorithm` is one of the tags `PODSVD`, `PODEVD`, `CotangentLiftSVD`, `CotangentLiftEVD`
and `ComplexSVD` of `ReducedBasisMethods`, which VlasovMethods does not export.

`parameters`, `s.paramspace` and `kwargs`, `k` or `tolerance`, go to `ReducedBasis` unchanged.
"""
function particle_basis(alg, s::Snapshots; parameters = NamedTuple(), kwargs...)
    Q = _snapshot_matrix(s.X)
    P = _snapshot_matrix(s.V)
    ReducedBasis(alg, [Q; P], parameters, s.paramspace; kwargs...)
end

"""
    field_basis(algorithm::Union{PODSVD, PODEVD}, s::Snapshots; parameters = NamedTuple(),
                kwargs...)

The POD `ReducedBasis` of the snapshot matrix of the accelerations `s.A`, whose `nd * nₚ` rows
match the position block of `particle_basis`. `algorithm` is the tag `PODSVD` or `PODEVD` of
`ReducedBasisMethods`, which VlasovMethods does not export. A PSD or complex tag would split the
rows of `A` into positions and velocities, so it throws a `MethodError`.

`parameters`, `s.paramspace` and `kwargs`, `k` or `tolerance`, go to `ReducedBasis` unchanged.
"""
function field_basis(alg::Union{PODSVD, PODEVD}, s::Snapshots; parameters = NamedTuple(),
        kwargs...)
    ReducedBasis(alg, _snapshot_matrix(s.A), parameters, s.paramspace; kwargs...)
end

### HDF5

const _SNAPSHOT_KEYS = (:X, :V, :A, :Φ, :W, :K, :M)

function Snapshots(h5::HDF5.H5DataStore, path::AbstractString = "/")
    group = h5[path]
    Snapshots((read(group[string(k)]) for k in _SNAPSHOT_KEYS)...,
        ParameterSpace(group, "parameterspace"))
end

function h5save(h5::HDF5.H5DataStore, s::Snapshots; path::AbstractString = "/")
    group = haskey(h5, path) ? h5[path] : create_group(h5, path)
    for k in _SNAPSHOT_KEYS
        group[string(k)] = getfield(s, k)
    end
    h5save(group, s.paramspace; path = "parameterspace")
end

function h5load(::Type{Snapshots}, h5::HDF5.H5DataStore; path::AbstractString = "/")
    Snapshots(h5, path)
end
