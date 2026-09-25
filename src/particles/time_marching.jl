
using LinearAlgebra: lu, ldiv!, norm, Diagonal
using PoissonSolvers: PBSpline, stiffnessmatrix, eval_deriv_PBSBasis, rhs_particles_PBSBasis
using ParticleMethods: ParticleList

struct IntegratorParameters{T}
    dt::T          # time step
    nₜ::Int        # number of time steps
    nₛ::Int        # number of saved time steps
    nₕ::Int        # number of basis functions
    nₚ::Int        # number of particles
    nparam::Int    # number of equation parameters sampled

    t::Vector{T}

    function IntegratorParameters(
            dt::T, nₜ::Int, nₛ::Int, nₕ::Int, nₚ::Int, nparam::Int) where {T}
        t = collect(range(0, stop = dt*nₜ, length = nₛ))
        new{T}(dt, nₜ, nₛ, nₕ, nₚ, nparam, t)
    end
end

function IntegratorParameters(h5::H5DataStore, path::AbstractString = "/")
    group = h5[path]
    IntegratorParameters(
        read(attributes(group)["dt"]),
        read(attributes(group)["nt"]),
        read(attributes(group)["ns"]),
        read(attributes(group)["nh"]),
        read(attributes(group)["np"]),
        read(attributes(group)["nparam"])
    )
end

function IntegratorParameters(fpath::AbstractString, path::AbstractString = "/")
    h5open(fpath, "r") do file
        IntegratorParameters(file, path)
    end
end

function Base.:(==)(ip1::IntegratorParameters{T1}, ip2::IntegratorParameters{T2}) where {
        T1, T2}
    T1 == T2 &&
        ip1.dt == ip2.dt &&
        ip1.nₜ == ip2.nₜ &&
        ip1.nₛ == ip2.nₛ &&
        ip1.nₕ == ip2.nₕ &&
        ip1.nₚ == ip2.nₚ &&
        ip1.nparam == ip2.nparam &&
        ip1.t == ip2.t
end

"""
save integrator parameters
"""
function h5save(h5::H5DataStore, IP::IntegratorParameters; path = "/")
    g = _create_group(h5, path)
    attributes(g)["dt"] = IP.dt
    attributes(g)["nt"] = IP.nₜ
    attributes(g)["ns"] = IP.nₛ
    attributes(g)["nh"] = IP.nₕ
    attributes(g)["np"] = IP.nₚ
    attributes(g)["nparam"] = IP.nparam
end

mutable struct ReducedIntegratorCache{T}
    zₓ::Vector{T}
    zᵥ::Vector{T}
    zₐ::Vector{T}
    w::Vector{T}
end

function ReducedIntegratorCache(IP::IntegratorParameters{T}, k::Int) where {T}
    ReducedIntegratorCache(zeros(T, k), # zₓ
        zeros(T, k), # zᵥ
        zeros(T, k), # zₐ
        zeros(T, IP.nₚ) # w
    )
end

function save_solution(SS, IC, Ψ, efield, w, p, t, ts, nsave, save = true)
    if save && ts % nsave == 0
        # advance time step index
        ts = ts+1

        # create views
        x = @view SS.X[1, :, ts, p]
        v = @view SS.V[1, :, ts, p]

        # reconstruct high fidelity solution
        mul!(x, Ψ, IC.zₓ)
        mul!(v, Ψ, IC.zᵥ)

        # solve for potential and copy efield coefficients
        update!(efield, IC.zₓ, w, t)
        SS.Φ[1, :, ts, p] .= coefficients(efield)

        # diagnostics
        SS.W[ts, p] = energy(efield)
        SS.K[ts, p] = dot(v, Diagonal(w), v) / 2
        SS.M[ts, p] = dot(w, v)
    end

    return ts
end

function reduced_integrate_vp(P₀::ParticleList{T},
        Ψ::Array{T},
        params::NamedTuple,
        efield::ReducedElectricField,
        SS, #::Snapshots,
        IP::IntegratorParameters{T},
        IC::ReducedIntegratorCache{T},
        p;
        save = true) where {T}

    # K needs to already be augmented for boundary conditions
    nsave = div(IP.nₜ, IP.nₛ-1)

    # initial conditions
    IC.zₓ .= Ψ' * vec(P₀.x)
    IC.zᵥ .= Ψ' * vec(P₀.v)
    IC.w .= vec(P₀.w)

    # effective timestep
    Δt = IP.dt * params.χ

    # save initial conditions
    ts = save_solution(SS, IC, Ψ, efield, IC.w, p, 0.0, 0, nsave, save)

    for it in 1:IP.nₜ
        # compute time
        t = it * IP.dt

        # half an advection step
        IC.zₓ .+= 0.5 .* Δt .* IC.zᵥ

        # evaluate electric field
        efield(IC.zₐ, IC.zₓ, IC.w, t)

        # acceleration step
        IC.zᵥ .+= Δt .* IC.zₐ

        # half an advection step
        IC.zₓ .+= 0.5 .* Δt .* IC.zᵥ

        ts = save_solution(SS, IC, Ψ, efield, IC.w, p, t, ts, nsave, save)
    end
end
