
struct Snapshots{T <: Number}
    X::Array{T, 4}
    V::Array{T, 4}
    A::Array{T, 4}
    Φ::Array{T, 4}
    W::Array{T, 2}
    K::Array{T, 2}
    M::Array{T, 2}

    function Snapshots(X::AbstractArray{T}, V::AbstractArray{T}, A::AbstractArray{T},
            Φ::AbstractArray{T}, W::AbstractArray{T},
            K::AbstractArray{T}, M::AbstractArray{T}) where {T}
        new{T}(X, V, A, Φ, W, K, M)
    end
end

function Snapshots(T, nd, np, nh, nt, nparam)
    X = zeros(T, nd, np, nt, nparam)
    V = zeros(T, nd, np, nt, nparam)
    A = zeros(T, nd, np, nt, nparam)
    Φ = zeros(T, nd, nh, nt, nparam)
    W = zeros(T, nt, nparam)
    K = zeros(T, nt, nparam)
    M = zeros(T, nt, nparam)
    Snapshots(X, V, A, Φ, W, K, M)
end

function Snapshots(T, ip::IntegratorParameters)
    Snapshots(T, 1, ip.nₚ, ip.nₕ, ip.nₛ, ip.nparam)
end

function Snapshots(particles::ParticleList{T}, poisson::PoissonSolver{T}, nt, pspace::ParameterSpace) where {T}
    Snapshots(T, 1, length(particles), length(poisson), nt, length(pspace))
end

function Base.:(==)(s1::Snapshots, s2::Snapshots)
    (
        s1.X == s2.X
        && s1.V == s2.V
        && s1.A == s2.A
        && s1.Φ == s2.Φ
        && s1.W == s2.W
        && s1.K == s2.K
        && s1.M == s2.M)
end

function Snapshots(h5::H5DataStore, path::AbstractString = "/")
    group = h5[path]

    X = read(group["X"])
    V = read(group["V"])
    A = read(group["A"])
    Φ = read(group["Φ"])
    W = read(group["W"])
    K = read(group["K"])
    M = read(group["M"])

    Snapshots(X, V, A, Φ, W, K, M)
end

function Snapshots(fpath::AbstractString, path::AbstractString = "/")
    h5open(fpath, "r") do file
        Snapshots(file, path)
    end
end

function h5save(h5::H5DataStore, TS::Snapshots; path::AbstractString = "/")
    group = _create_group(h5, path)

    group["X"] = TS.X
    group["V"] = TS.V
    group["A"] = TS.A
    group["Φ"] = TS.Φ
    group["W"] = TS.W
    group["K"] = TS.K
    group["M"] = TS.M
end

function h5load(::Type{Snapshots}, h5::H5DataStore; path::AbstractString = "/")
    Snapshots(h5, path)
end
