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

    function Snapshots(X::AbstractArray{T, 4}, V::AbstractArray{T, 4},
            A::AbstractArray{T, 4}, Φ::AbstractArray{T, 3}, W::AbstractArray{T, 2},
            K::AbstractArray{T, 2}, M::AbstractArray{T, 2},
            paramspace::PS) where {T, PS <: ParameterSpace}
        size(V) == size(X) || throw(ArgumentError(
            "V needs the size $(size(X)) of X, got $(size(V))"))
        size(A) == size(X) || throw(ArgumentError(
            "A needs the size $(size(X)) of X, got $(size(A))"))
        size(Φ)[2:3] == size(X)[3:4] || throw(ArgumentError(
            "Φ needs the time and sample extents $(size(X)[3:4]) of X, got $(size(Φ)[2:3])"))
        for (name, D) in (("W", W), ("K", K), ("M", M))
            size(D) == size(X)[3:4] || throw(ArgumentError(
                "$name needs the time and sample extents $(size(X)[3:4]) of X, got $(size(D))"))
        end
        size(X, 4) == length(paramspace) || throw(ArgumentError(
            "the snapshots hold $(size(X, 4)) samples, the parameter space holds " *
            "$(length(paramspace))"))
        new{T, PS}(X, V, A, Φ, W, K, M, paramspace)
    end
end

"""
    Snapshots(T, nd, nₚ, nh, nₛ, paramspace)

The zero `Snapshots{T}` of `nd` components, `nₚ` particles, `nh` potential degrees of freedom,
`nₛ` saved times and `length(paramspace)` samples, which a sweep fills as it runs. `nh = 0` is the
sweep of a model with no potential, such as one of `VlasovWithExternalFields`; the HDF5 round trip
of a `Snapshots` with no potential is not claimed.
"""
function Snapshots(T, nd, np, nh, nt, paramspace::ParameterSpace)
    nparam = length(paramspace)
    X = zeros(T, nd, np, nt, nparam)
    V = zeros(T, nd, np, nt, nparam)
    A = zeros(T, nd, np, nt, nparam)
    Φ = zeros(T, nh, nt, nparam)
    W = zeros(T, nt, nparam)
    K = zeros(T, nt, nparam)
    M = zeros(T, nt, nparam)
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
function Snapshots(particles::ParticleList{T}, poisson::PoissonSolver{T}, nt,
        paramspace::ParameterSpace) where {T}
    nd = size(particles.x, 1)
    size(particles.v, 1) == nd || throw(ArgumentError(
        "the particle list has $(nd) position components and $(size(particles.v, 1)) velocity " *
        "components, which the snapshots cannot pair"))
    Snapshots(T, nd, length(particles), length(poisson), nt, paramspace)
end

function Base.:(==)(s1::Snapshots, s2::Snapshots)
    (
        s1.X == s2.X
        && s1.V == s2.V
        && s1.A == s2.A
        && s1.Φ == s2.Φ
        && s1.W == s2.W
        && s1.K == s2.K
        && s1.M == s2.M
        && s1.paramspace == s2.paramspace)
end

### the snapshot matrices the reduced bases are built from

"""
    _snapshot_matrix(Y)

The `size(Y, 1) * size(Y, 2) × size(Y, 3) * size(Y, 4)` snapshot matrix of the array `Y` of one
saved quantity, whose column `(p - 1) * size(Y, 3) + t` is `vec(Y[:, :, t, p])`: the
`(size(Y, 1), size(Y, 2))` matrix of the saved time `t` of the sample `p`, in the
component-fastest order of `vec`.

It is a `reshape` of `Y`: it shares the memory of `Y` and copies nothing. A reshape to
`(size(Y, 2), :)` instead is valid only for `size(Y, 1) == 1`.
"""
function _snapshot_matrix(Y::AbstractArray)
    return reshape(Y, size(Y, 1) * size(Y, 2), size(Y, 3) * size(Y, 4))
end

"""
    _particle_snapshots(s::Snapshots)

The `2 * nd * nₚ × nₛ * nparam` particle snapshot matrix `[Q; P]` of `s`: the `nd * nₚ` rows of
`_snapshot_matrix(s.X)` followed by the `nd * nₚ` rows of `_snapshot_matrix(s.V)`.

The cotangent lift and the complex SVD read that stacking and split it at the row `nd * nₚ`, and
the POD routes reduce the whole matrix; neither reads an interleaved `[q₁; p₁; q₂; p₂; …]`. Both
accept any even number of rows, so an interleaved matrix would be reduced without an error.
"""
_particle_snapshots(s::Snapshots) = vcat(_snapshot_matrix(s.X), _snapshot_matrix(s.V))

"""
    _field_snapshots(s::Snapshots)

The `nd * nₚ × nₛ * nparam` field snapshot matrix of `s`: the acceleration `A` of each saved time
and sample, in the component-fastest order of `vec`, one column each.

The field basis uses `A`, not the potential: the reduced field couples it to the position
block of the particle basis, which has the `nd * nₚ` rows of `X`.
"""
_field_snapshots(s::Snapshots) = _snapshot_matrix(s.A)

### the reduced bases of a sweep

"""
    particle_basis(algorithm, s::Snapshots; parameters = NamedTuple(), k = nothing,
                   tolerance = nothing)

The `ReducedBasis` of the particle snapshots of `s`, built by `algorithm`, which is any of
`PODSVD`, `PODEVD`, `CotangentLiftSVD`, `CotangentLiftEVD` and `ComplexSVD`.

The basis is built from the stacked `[Q; P]` matrix `_particle_snapshots(s)` of `2 * nd * nₚ`
rows, whose first half holds the positions `X` and whose second half the velocities `V`: exactly
the snapshot the two PSD routes and the complex SVD split, and the whole matrix the two POD routes
reduce.

`parameters`, `s.paramspace` and exactly one of `k` and `tolerance` are passed to `ReducedBasis`
unchanged, so its `ArgumentError`s are the guards of this function.
"""
function particle_basis(alg, s::Snapshots; parameters = NamedTuple(), k = nothing,
        tolerance = nothing)
    return ReducedBasis(alg, _particle_snapshots(s), parameters, s.paramspace; k, tolerance)
end

"""
    field_basis(algorithm::Union{PODSVD, PODEVD}, s::Snapshots; parameters = NamedTuple(),
                k = nothing, tolerance = nothing)

The `ReducedBasis` of the field snapshots of `s`, built by `algorithm`.

The basis is built from `_field_snapshots(s)`, the `nd * nₚ × nₛ * nparam` matrix of the
accelerations, and not from the potential `Φ`, whose `nh` degrees of freedom are another space.
Its `nd * nₚ` rows match the position block of the stacked particle basis.

`algorithm` is one of the two POD tags: a PSD or complex tag would split the `nd * nₚ` rows of
`A` in the middle and reduce the halves as positions and velocities, so it is refused with a
`MethodError` rather than reduced. `parameters`, `s.paramspace` and exactly one of `k` and
`tolerance` are passed to `ReducedBasis` unchanged.
"""
function field_basis(alg::Union{PODSVD, PODEVD}, s::Snapshots; parameters = NamedTuple(),
        k = nothing, tolerance = nothing)
    return ReducedBasis(alg, _field_snapshots(s), parameters, s.paramspace; k, tolerance)
end

### HDF5

function Snapshots(h5::HDF5.H5DataStore, path::AbstractString = "/")
    group = h5[path]

    X = read(group["X"])
    V = read(group["V"])
    A = read(group["A"])
    Φ = read(group["Φ"])
    W = read(group["W"])
    K = read(group["K"])
    M = read(group["M"])
    paramspace = ParameterSpace(group, "parameterspace")

    Snapshots(X, V, A, Φ, W, K, M, paramspace)
end

function Snapshots(fpath::AbstractString, path::AbstractString = "/")
    h5open(fpath, "r") do file
        Snapshots(file, path)
    end
end

function h5save(h5::HDF5.H5DataStore, s::Snapshots; path::AbstractString = "/")
    group = haskey(h5, path) ? h5[path] : create_group(h5, path)

    group["X"] = s.X
    group["V"] = s.V
    group["A"] = s.A
    group["Φ"] = s.Φ
    group["W"] = s.W
    group["K"] = s.K
    group["M"] = s.M
    h5save(group, s.paramspace; path = "parameterspace")
end

function h5load(::Type{Snapshots}, h5::HDF5.H5DataStore; path::AbstractString = "/")
    Snapshots(h5, path)
end
