
using LinearAlgebra: dot


struct VPIntegratorParameters{T}
    dt::T          # time step
    nₜ::Int        # number of time steps
    nₛ::Int        # number of saved time steps
    nₕ::Int        # number of basis functions
    nₚ::Int        # number of particles

    t::Vector{T}

    function VPIntegratorParameters(dt::T, nₜ::Int, nₛ::Int, nₕ::Int, nₚ::Int) where {T}
        t = collect(range(0, stop=dt*nₜ, length=nₛ))
        new{T}(dt,nₜ,nₛ,nₕ,nₚ,t)
    end
end


mutable struct VPIntegratorCache{T}
    x::Vector{T}
    v::Vector{T}
    a::Vector{T}
    w::Vector{T}

    ρ::Vector{T}
    ϕ::Vector{T}
    rhs::Vector{T}

    X::Matrix{T}
    V::Matrix{T}
    A::Matrix{T}
    Φ::Matrix{T}

    W::Vector{T}
    K::Vector{T}
    M::Vector{T}
end

VPIntegratorCache(IP::VPIntegratorParameters{T}) where {T} = VPIntegratorCache(
                                                                zeros(T,IP.nₚ), # x
                                                                zeros(T,IP.nₚ), # v
                                                                zeros(T,IP.nₚ), # a
                                                                zeros(T,IP.nₚ), # w
                                                                zeros(T,IP.nₕ), # ρ
                                                                zeros(T,IP.nₕ), # ϕ
                                                                zeros(T,IP.nₕ), # rhs
                                                                zeros(T,IP.nₚ,IP.nₛ), # X
                                                                zeros(T,IP.nₚ,IP.nₛ), # V
                                                                zeros(T,IP.nₚ,IP.nₛ), # A
                                                                zeros(T,IP.nₕ,IP.nₛ), # Φ
                                                                zeros(T,IP.nₛ), # W
                                                                zeros(T,IP.nₛ), # K
                                                                zeros(T,IP.nₛ)  # M
                                                            )

function save_timestep!(IC::VPIntegratorCache, efield::ElectricField, ts)
    IC.X[:,ts] .= IC.x
    IC.V[:,ts] .= IC.v
    IC.A[:,ts] .= IC.a
    IC.Φ[:,ts] .= coefficients(efield)

    IC.W[ts] = energy(efield)
    IC.K[ts] = dot(IC.w .* IC.v, IC.v) / 2
    IC.M[ts] = dot(IC.w, IC.v)
end


function integrate_vp!(P::ParticleList{T},
                       efield::ElectricField,
                       parameters::NamedTuple,
                       IP::VPIntegratorParameters{T},
                       IC::VPIntegratorCache{T} = VPIntegratorCache(IP);
                       save = true) where {T}

    nsave = div(IP.nₜ, IP.nₛ-1)

    # effective timestep
    Δt = IP.dt * parameters.χ

    # initial conditions
    IC.x .= P.x[1,:]
    IC.v .= P.v[1,:]
    IC.w .= P.w[1,:]

    # save initial conditions
    if save
        update!(efield, IC.x, IC.w, 0.0)
        save_timestep!(IC, efield, 1)
    end

    ts = 1
    for it in 1:IP.nₜ
        # compute time
        t = (it+1) * IP.dt

        # half an advection step
        IC.x .+= 0.5 .* Δt .* IC.v

        # evaluate electric field
        efield(IC.a, IC.x, IC.w, t)

        # acceleration step
        IC.v .+= Δt .* IC.a

        # half an advection step
        IC.x .+= 0.5 .* Δt .* IC.v

        if save && it % nsave == 0
            update!(efield, IC.x, IC.w, t)
            save_timestep!(IC, efield, ts+1)
            ts += 1
        end
    end

    # P.x .= IC.x
    # P.v .= IC.v
end
