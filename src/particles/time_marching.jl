# The reduced particle march of the analytic quadratic model: the leapfrog over the
# `ElectricField` interface of `src/electric_field.jl`, with the χ-scaled timestep.
#
# `IntegratorParameters` and its HDF5 pair, `ReducedIntegratorCache`, the leapfrog step, the save of
# one state, and the march. `leapfrog_step!` and `save_solution` are internal;
# `IntegratorParameters`, `ReducedIntegratorCache` and `reduced_integrate_vp` are not.

using ParticleMethods: ParticleList

"""
    IntegratorParameters(dt, nₜ, nₛ, nₕ, nₚ, nparam)

The parameters of a reduced sweep: the timestep `dt`, the number of time steps `nₜ`, the number of
saved times `nₛ` — the initial state included — the number of potential degrees of freedom `nₕ`,
the number of particles `nₚ` and the number of sampled equation parameters `nparam`. `t` holds
`nₛ` times, evenly spaced from `0` to `dt·nₜ`.

Every count is `Int`, so an `Int32` count is a `MethodError`. The element type `T` of `dt` is
shared by `P₀`, `Ψ`, `SS` and `IC` in [`reduced_integrate_vp`](@ref), so mixed storage types are
a `MethodError`. The timestep scale `params.χ` is converted to `T` by the march.

The HDF5 pair is `h5save(h5, IP; path = …)` and `IntegratorParameters(h5, path)`, whose round trip is
asserted `==`. `h5save` makes a group that no path holds yet with HDF5's own `create_group`.
"""
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

function IntegratorParameters(h5::HDF5.H5DataStore, path::AbstractString = "/")
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

function h5save(h5::HDF5.H5DataStore, IP::IntegratorParameters; path = "/")
    group = haskey(h5, path) ? h5[path] : create_group(h5, path)
    attributes(group)["dt"] = IP.dt
    attributes(group)["nt"] = IP.nₜ
    attributes(group)["ns"] = IP.nₛ
    attributes(group)["nh"] = IP.nₕ
    attributes(group)["np"] = IP.nₚ
    attributes(group)["nparam"] = IP.nparam
end

"""
    ReducedIntegratorCache(IP::IntegratorParameters, k)

The work buffers of the reduced march of basis size `k`: the reduced positions `zₓ`, velocities
`zᵥ` and accelerations `zₐ`, each `k` long, and the `nₚ` particle weights `w`. One cache serves a
whole march, which is why the march takes it as an argument.

The buffers follow the element type of `IP`, and their length is not checked against a basis: the
march's `k` is the column count of the basis it is handed.
"""
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

# `dot(w, sum(abs2, v; dims = 1)) / 2`, which is the kinetic energy for `nd == 1`: the `dims`
# reduction allocates its `(1, nₚ)` result, so the sum is a loop.
@inline function _kinetic_energy(w::AbstractVector{T}, v::AbstractMatrix{T}) where {T}
    s = zero(T)
    @inbounds for j in axes(v, 2)
        s += w[j] * abs2(v[1, j])
    end
    return s / 2
end

"""
    leapfrog_step!(IC::ReducedIntegratorCache, efield, Δt, t)

One step of the reduced leapfrog: half an advection at the scaled `Δt`, the acceleration of the
reduced field at the new position, and the second half advection. The symmetric composition is
second order; its leading phase-error contribution to check 1's relative `L²` error over the
horizon `T` is `β·T·(βΔt)²/(24√3)`.

It allocates nothing and is inferred: the cache, the field and both times are concrete. `Δt` is the
already scaled timestep — `Δt = IP.dt·T(params.χ)` — so half of it is `T(0.5)·Δt`, the scale being
outside the half.
"""
@inline function leapfrog_step!(
        IC::ReducedIntegratorCache{T}, efield::ReducedElectricField,
        Δt::T, t) where {T}
    half = T(0.5) * Δt

    IC.zₓ .+= half .* IC.zᵥ
    efield(IC.zₐ, IC.zₓ, IC.w, t)
    IC.zᵥ .+= Δt .* IC.zₐ
    IC.zₓ .+= half .* IC.zᵥ

    return IC
end

"""
    save_solution(SS, IC, Ψ, efield, w, p, t, it, nsave, save = true)

Write the state of the step `it` into the sample `p` of the snapshots `SS` when `it` is a multiple
of `nsave`: the positions and the velocities lifted by `Ψ`, the potential's coefficients, and the
three diagnostics. The saved index is `div(it, nsave) + 1`, so step `0` writes the initial state
and `nₜ/nsave + 1 = nₛ` states are written in all.

**The step index decides the save.** Step `it` is written when `nsave` divides `it`, so the save
index `div(it, nsave) + 1` runs over `1:nₛ` exactly once and no saved state is passed over.

The scalar momentum `M` is `dot(w, v)`, which is the momentum for `nd == 1` alone; the march refuses
another component count. It allocates nothing.
"""
function save_solution(
        SS::Snapshots{T}, IC::ReducedIntegratorCache{T}, Ψ::AbstractMatrix{T},
        efield::ReducedElectricField, w::AbstractVector{T}, p, t::T,
        it::Int, nsave::Int, save = true) where {T}
    if save && it % nsave == 0
        ts = div(it, nsave) + 1

        # views of the `(nd, nₚ)` positions and velocities
        x = @view SS.X[:, :, ts, p]
        v = @view SS.V[:, :, ts, p]

        # reconstruct high fidelity solution
        mul!(vec(x), Ψ, IC.zₓ)
        mul!(vec(v), Ψ, IC.zᵥ)

        # solve for potential and copy efield coefficients
        update!(efield, IC.zₓ, w, t)
        SS.Φ[:, ts, p] .= coefficients(efield)

        # diagnostics
        SS.W[ts, p] = energy(efield)
        SS.K[ts, p] = _kinetic_energy(w, v)
        SS.M[ts, p] = dot(w, v)
    end

    return nothing
end

# The guards of the march, each raised before the loop. The first three make `nsave` a positive
# integer that divides `nₜ` and reaches `nₛ`; the fourth keeps the saves inside the snapshots; the
# fifth is what `M = dot(w, v)` needs; the sixth refuses a zero-width basis, whose reduced state is
# empty; and the last is written so that `NaN` is refused too.
function _check_march(SS::Snapshots, Ψ::AbstractMatrix, IP::IntegratorParameters, params)
    IP.nₜ ≥ 1 || throw(ArgumentError(
        "the march needs at least one time step, got nₜ = $(IP.nₜ)"))
    IP.nₛ ≥ 2 || throw(ArgumentError(
        "the march saves the initial state and at least one more, got nₛ = $(IP.nₛ)"))
    IP.nₜ % (IP.nₛ - 1) == 0 || throw(ArgumentError(
        "nₜ = $(IP.nₜ) is not a multiple of nₛ - 1 = $(IP.nₛ - 1), so the saves would either " *
        "write past nₛ or drop the final state"))
    size(SS.X, 3) == IP.nₛ || throw(DimensionMismatch(
        "the snapshots hold $(size(SS.X, 3)) saved times; the parameters ask for $(IP.nₛ)"))
    size(SS.X, 1) == 1 || throw(ArgumentError(
        "the scalar momentum is defined for nd == 1; the snapshots hold $(size(SS.X, 1)) " *
        "components"))
    k = size(Ψ, 2)
    1 ≤ k ≤ size(Ψ, 1) || throw(ArgumentError(
        "the basis is $(size(Ψ, 1)) × $(size(Ψ, 2)); a reduced state is between 1 and " *
        "$(size(Ψ, 1)) long"))
    params.χ > 0 || throw(ArgumentError(
        "the χ scale is a positive timestep scale, got χ = $(params.χ)"))
    return nothing
end

"""
    reduced_integrate_vp(P₀, Ψ, params, efield, SS, IP, IC, p; save = true)

The reduced leapfrog march of the analytic quadratic model: the reduced state of `P₀` — the two
coordinates `Ψ'·vec(P₀.x)` and `Ψ'·vec(P₀.v)` — advanced `IP.nₜ` steps of `IP.dt` by
[`leapfrog_step!`](@ref), with the reduced acceleration of `efield` at every half step, and the
state of every `nsave`-th step written into the sample `p` of the snapshots `SS` by
[`save_solution`](@ref). `Ψ` is the position basis, and the same basis lifts the velocities: the
cotangent lift's `Φ`.

`params` is the caller's, and it must hold `χ`. **The χ scale is a timestep scale and nothing
else:** the step is `Δt = IP.dt·T(params.χ)`, converted to the march's element type before the
product, and the time handed to the field is `t = it·IP.dt`, unscaled. So a time-independent field
marched with `(dt, χ = 2)` equals the same field marched with `(2dt, χ = 1)` state for state.

It returns `nothing`; the states are in `SS`, and the integration parameters hold the times.

# Throws

- `ArgumentError` when `nₜ < 1`; when `nₛ < 2`; when `nₛ - 1` does not divide `nₜ`, so that the
  saves would write past `nₛ` or drop the final state; when the snapshots carry another number of
  components than one; when the basis has a zero column count or more columns than rows; and when
  `params.χ` is not `> 0`, which refuses `0`, a negative and `NaN` alike.
- `DimensionMismatch` when the snapshots do not hold `IP.nₛ` saved times.
- The `NamedTuple`'s own error when `params` holds no `χ`. There is no default.

# Cost

The setup is not claimed. The loop's two calls — [`leapfrog_step!`](@ref) and
[`save_solution`](@ref) — are inferred and allocate nothing, which the tests assert through
fixed-arity barriers.
"""
function reduced_integrate_vp(P₀::ParticleList{T},
        Ψ::AbstractMatrix{T},
        params::NamedTuple,
        efield::ReducedElectricField,
        SS::Snapshots{T},
        IP::IntegratorParameters{T},
        IC::ReducedIntegratorCache{T},
        p;
        save = true) where {T}
    _check_march(SS, Ψ, IP, params)

    nsave = div(IP.nₜ, IP.nₛ - 1)

    # initial conditions: the reduced coordinates of the particle state
    IC.zₓ .= Ψ' * vec(P₀.x)
    IC.zᵥ .= Ψ' * vec(P₀.v)
    IC.w .= vec(P₀.w)

    # effective timestep: the χ scale, converted to `T` before the product
    Δt = IP.dt * T(params.χ)

    # the initial state, at `zero(T)` and not `0.0`
    save_solution(SS, IC, Ψ, efield, IC.w, p, zero(T), 0, nsave, save)

    for it in 1:IP.nₜ
        # the time of the field stays unscaled
        t = it * IP.dt

        leapfrog_step!(IC, efield, Δt, t)

        save_solution(SS, IC, Ψ, efield, IC.w, p, t, it, nsave, save)
    end

    return nothing
end
