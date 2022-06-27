
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

function save_timestep!(IC::VPIntegratorCache, poisson::PoissonSolver, ts, χ)
    IC.X[:,ts] .= IC.x
    IC.V[:,ts] .= IC.v
    IC.A[:,ts] .= IC.a
    IC.Φ[:,ts] .= poisson.ϕ

    IC.W[ts] = dot(poisson.ϕ, poisson.S, poisson.ϕ) / 2 * χ^2
    IC.K[ts] = dot(IC.w .* IC.v, IC.v) / 2
    IC.M[ts] = dot(IC.w, IC.v)
end

function solve_potential!(IC::VPIntegratorCache, poisson::PoissonSolver, ts, χ, Φₑₓₜ, given_phi)
    if given_phi
        IC.ϕ .= Φₑₓₜ[:,ts]
    else
        solve!(poisson, IC.x, IC.w)
        poisson.ϕ ./= χ^2
        IC.ϕ .= poisson.ϕ
    end
end


function integrate_vp!(P::ParticleList{T},
                       poisson::PoissonSolver{T},
                       parameters::NamedTuple,
                       IP::VPIntegratorParameters{T},
                       IC::VPIntegratorCache{T} = VPIntegratorCache(IP);
                       given_phi = false,
                       Φₑₓₜ::Array{T} = zeros(T, IP.nₕ, IP.nₚ),
                       save = true) where {T}

    nsave = div(IP.nₜ, IP.nₛ-1)

    if given_phi
        @assert IP.nₛ == IP.nₜ + 1
    end

    # simulation parameter
    χ = parameters.χ

    # initial conditions
    IC.x .= P.x
    IC.v .= P.v
    IC.w .= P.w

    # save initial conditions
    if save
        solve_potential!(IC, poisson, 1, χ, Φₑₓₜ, given_phi)
        save_timestep!(IC, poisson, 1, χ)
    end

    tₛ = 1
    for t in 1:IP.nₜ
        # half an advection step
        IC.x .+= 0.5 .* IP.dt .* IC.v .* χ

        # solve for potential
        solve_potential!(IC, poisson, t+1, χ, Φₑₓₜ, given_phi)

        # evaluate electric field
        eval_field!(IC.a, poisson, IC.x)

        # acceleration step
        IC.v .+= IP.dt .* IC.a .* χ

        # half an advection step
        IC.x .+= 0.5 .* IP.dt .* IC.v .* χ

        if save && t % nsave == 0
            solve_potential!(IC, poisson, t+1, χ, Φₑₓₜ, given_phi)
            save_timestep!(IC, poisson, tₛ+1, χ)
            tₛ += 1
        end
    end

    # P.x .= IC.x
    # P.v .= IC.v
end
