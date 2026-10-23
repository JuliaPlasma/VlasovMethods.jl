# The particle marches over the `ElectricField` interface of `src/electric_field.jl`: the
# full-order march of `integrate_vp`, on any field, and the reduced march of the analytic quadratic
# model, `reduced_integrate_vp`, both with the χ-scaled timestep.
#
# `IntegratorParameters` and its HDF5 pair, `ReducedIntegratorCache`, the leapfrog steps, the saves
# of one state, and the marches. The steps and the saves are internal; `IntegratorParameters`,
# `ReducedIntegratorCache`, `integrate_vp` and `reduced_integrate_vp` are not.

using ParticleMethods: ParticleList

"""
    IntegratorParameters(dt, nₜ, nₛ, nₕ, nₚ, nparam)

The parameters of a reduced sweep: the timestep `dt`, the number of time steps `nₜ`, the number of
saved times `nₛ` — the initial state included — the number of potential degrees of freedom `nₕ`,
the number of particles `nₚ` and the number of sampled equation parameters `nparam`. `t` holds
`nₛ` times, evenly spaced from `0` to `dt·nₜ`.

The constructor throws an `ArgumentError` when `nₜ < 1`, when `nₛ < 2`, and when `nₛ - 1` does not
divide `nₜ`, so that the saves of the march would write past `nₛ` or drop the final state. The
HDF5 constructor applies the same checks to the stored counts.

Every count is `Int`, so an `Int32` count is a `MethodError`. The element type `T` of `dt` is
shared by `P₀`, `Ψ`, `SS` and `IC` in [`reduced_integrate_vp`](@ref), so mixed storage types are
a `MethodError`. The timestep scale `params.χ` is converted to `T` by the march.

The HDF5 pair is `h5save(h5, IP; path = …)` and `IntegratorParameters(h5, path)`. `h5save` writes
into the group at `path`, and creates it when it does not exist.
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
        nₜ ≥ 1 || throw(ArgumentError("a sweep takes at least one time step, got nₜ = $nₜ"))
        nₛ ≥ 2 || throw(ArgumentError(
            "a sweep saves the initial state and at least one more, got nₛ = $nₛ"))
        nₜ % (nₛ - 1) == 0 || throw(ArgumentError(
            "nₜ = $nₜ is not a multiple of nₛ - 1 = $(nₛ - 1), so the saves would either " *
            "write past nₛ or drop the final state"))
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
    return group
end

"""
    ReducedIntegratorCache(IP::IntegratorParameters, k)

The work buffers of the reduced march of basis size `k`: the reduced positions `zₓ`, velocities
`zᵥ` and accelerations `zₐ`, each `k` long, and the `nₚ` particle weights `w`. One cache serves a
whole march, which is why the march takes it as an argument.

The buffers follow the element type of `IP`, and their length is not checked against a basis: the
march's `k` is the column count of the basis it is handed.
"""
struct ReducedIntegratorCache{T}
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

"""
    leapfrog_step!(IC::ReducedIntegratorCache, efield, Δt, t)

One step `Δt` of the reduced leapfrog: half a drift, the kick of the reduced field at the
half-drifted position and the time `t`, and the second half drift. The symmetric composition is
second order. It is inferred and allocates nothing.
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
    save_solution(SS, IC, Ψ, efield, p, t, ts)

Write the state of the cache `IC` at the time `t` into the save `ts` of the sample `p` of the
snapshots `SS`: the positions and the velocities lifted by `Ψ`, the potential's coefficients, and
the three diagnostics. The accelerations `SS.A` are not written.

The kinetic energy `K = ½∑ⱼwⱼvⱼ²` and the scalar momentum `M = ∑ⱼwⱼvⱼ` are those of one component,
`nd == 1`, which the march requires. It is inferred and allocates nothing.
"""
function save_solution(
        SS::Snapshots{T}, IC::ReducedIntegratorCache{T}, Ψ::AbstractMatrix{T},
        efield::ReducedElectricField, p, t::T, ts::Int) where {T}
    # the lifted positions and velocities of the one component
    x = @view SS.X[1, :, ts, p]
    v = @view SS.V[1, :, ts, p]
    mul!(x, Ψ, IC.zₓ)
    mul!(v, Ψ, IC.zᵥ)

    # solve for the potential and copy its coefficients
    update!(efield, IC.zₓ, IC.w, t)
    SS.Φ[:, ts, p] .= coefficients(efield)

    # diagnostics
    SS.W[ts, p] = energy(efield)
    SS.K[ts, p] = dot(v, Diagonal(IC.w), v) / 2
    SS.M[ts, p] = dot(IC.w, v)

    return nothing
end

# The guards of the march, each raised before the loop; the counts of `IP` are checked by its
# constructor. The first keeps the saves inside the snapshots; the second is what `M = dot(w, v)`
# needs; the third refuses a zero-width basis, whose reduced state is empty; and the last is written
# so that `NaN` is refused too.
function _check_march(SS::Snapshots, Ψ::AbstractMatrix, IP::IntegratorParameters, params)
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
[`leapfrog_step!`](@ref), with one kick of the reduced acceleration of `efield` per step. With
`nsave = nₜ/(nₛ - 1)`, the state of every `nsave`-th step `it`, the initial state included, is
written into the save `div(it, nsave) + 1` of the sample `p` of the snapshots `SS` by
[`save_solution`](@ref). `Ψ` is the position basis, and the same basis lifts the velocities: the
cotangent lift's `Φ`.

`params` is the caller's, and it must hold `χ`. **The χ scale is a timestep scale and nothing
else:** the step is `Δt = IP.dt·T(params.χ)`, converted to the march's element type before the
product, and the times handed to the field are unscaled: the kick of step `it` is at its midpoint
`(it - ½)·IP.dt`, and the save of step `it` at `it·IP.dt`. So a time-independent field marched
with `(dt, χ = 2)` equals the same field marched with `(2dt, χ = 1)` state for state.

**`params.χ` is the one source of the scale.** When `efield` is a `ReducedElectricField` over a
[`ScaledField`](@ref), the `χ` the `ScaledField` was built with must equal `params.χ`, and the
march takes the field as it is. Any other field is marched as it is.

It returns `nothing`; the states are in `SS`, and the integration parameters hold the times.

# Throws

- `ArgumentError` when the snapshots carry another number of components than one; when the basis
  has a zero column count or more columns than rows; when `params.χ` is not `> 0`, which
  refuses `0`, a negative and `NaN` alike; and when `efield` is over a `ScaledField` built with
  another `χ` than `params.χ`, naming both. The counts of `IP` are checked by
  [`IntegratorParameters`](@ref) itself.
- `DimensionMismatch` when the snapshots do not hold `IP.nₛ` saved times.
- The `NamedTuple`'s own error when `params` holds no `χ`. There is no default.

# Cost

The loop's two calls — [`leapfrog_step!`](@ref) and [`save_solution`](@ref) — are inferred and
allocate nothing.
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

    # the field at the scale of the parameters
    field = _with_scale(efield, params.χ)

    # initial conditions: the reduced coordinates of the particle state
    mul!(IC.zₓ, Ψ', vec(P₀.x))
    mul!(IC.zᵥ, Ψ', vec(P₀.v))
    IC.w .= vec(P₀.w)

    # effective timestep: the χ scale, converted to `T` before the product
    Δt = IP.dt * T(params.χ)

    save && save_solution(SS, IC, Ψ, field, p, zero(T), 1)

    for it in 1:IP.nₜ
        # the times of the field are unscaled
        leapfrog_step!(IC, field, Δt, (it - T(0.5)) * IP.dt)

        if save && it % nsave == 0
            save_solution(SS, IC, Ψ, field, p, it * IP.dt, div(it, nsave) + 1)
        end
    end

    return nothing
end

"""
    full_leapfrog_step!(x, v, a, w, efield, Δt, t)

One step `Δt` of the full-order leapfrog: half a drift, the deposit of the half-drifted positions
and the kick of the full-order field at the time `t`, and the second half drift. `x` and `v` are
the particle positions and velocities, `w` their weights, and `a` the acceleration buffer the
field writes. The symmetric composition is second order. It is inferred and allocates nothing.
"""
@inline function full_leapfrog_step!(
        x::AbstractVector{T}, v::AbstractVector{T}, a::AbstractVector{T},
        w::AbstractVector{T}, efield::ElectricField, Δt::T, t) where {T}
    half = T(0.5) * Δt

    x .+= half .* v
    update!(efield, x, w, t)
    efield!(efield, a, x)
    v .+= Δt .* a
    x .+= half .* v

    return nothing
end

"""
    save_full_solution(SS, x, v, w, efield, p, t, ts)

Write the state of a full-order march at the time `t` into the save `ts` of the sample `p` of the
snapshots `SS`: the positions `x` and the velocities `v`, the potential's coefficients, and the
three diagnostics. The accelerations `SS.A` are not written.

The kinetic energy `K = ½∑ⱼwⱼvⱼ²` and the scalar momentum `M = ∑ⱼwⱼvⱼ` are those of a single
component, `nd == 1`, which the march requires. It is inferred and allocates nothing.
"""
function save_full_solution(
        SS::Snapshots{T}, x::AbstractVector{T}, v::AbstractVector{T}, w::AbstractVector{T},
        efield::ElectricField, p, t::T, ts::Int) where {T}
    # the saved state of the one component
    copyto!(@view(SS.X[1, :, ts, p]), x)
    copyto!(@view(SS.V[1, :, ts, p]), v)

    # solve for the potential and copy its coefficients
    update!(efield, x, w, t)
    SS.Φ[:, ts, p] .= coefficients(efield)

    # diagnostics
    SS.W[ts, p] = energy(efield)
    SS.K[ts, p] = dot(v, Diagonal(w), v) / 2
    SS.M[ts, p] = dot(w, v)

    return nothing
end

# The guards of the full-order march: the three of `_check_march` that name no basis, with `χ` also
# required to be finite, as `ScaledField` requires, and the particle count of the list, which the
# snapshots must hold room for.
function _check_full_march(SS::Snapshots, P₀, IP::IntegratorParameters, params)
    size(SS.X, 3) == IP.nₛ || throw(DimensionMismatch(
        "the snapshots hold $(size(SS.X, 3)) saved times; the parameters ask for $(IP.nₛ)"))
    size(SS.X, 1) == 1 || throw(ArgumentError(
        "the scalar momentum is defined for nd == 1; the snapshots hold $(size(SS.X, 1)) " *
        "components"))
    length(P₀) == size(SS.X, 2) || throw(DimensionMismatch(
        "the particle list holds $(length(P₀)) particles; the snapshots hold " *
        "$(size(SS.X, 2))"))
    isfinite(params.χ) && params.χ > 0 || throw(ArgumentError(
        "the χ scale is a finite, positive timestep scale, got χ = $(params.χ)"))
    return nothing
end

"""
    integrate_vp(P₀, params, efield, SS, IP, p; save = true)

The full-order leapfrog march of an [`ElectricField`](@ref), the same march as
[`reduced_integrate_vp`](@ref) without a basis: the state of `P₀` advanced `IP.nₜ` steps of
`IP.dt` by [`full_leapfrog_step!`](@ref), with one kick of the full field of `efield` per step. A
reduced model marches the [`ScaledField`](@ref) of the same field, so the two marches of one
comparison take the same field and the same `χ`, which is `params.χ`.
With `nsave = nₜ/(nₛ - 1)`, the state of every `nsave`-th step `it`, the initial state included,
is written into the save `div(it, nsave) + 1` of the sample `p` of the snapshots `SS` by
[`save_full_solution`](@ref).

The march holds its own buffers, so `efield` is the only state it shares with the caller: `P₀` is
copied and never written, and the march needs no cache argument.

`params` is the caller's, and it must hold `χ`. **The χ scale is a timestep scale and nothing
else:** the step is `Δt = IP.dt·T(params.χ)`, converted to the march's element type before the
product, and the times handed to the field are unscaled: the kick of step `it` is at its midpoint
`(it - ½)·IP.dt`, and the save of step `it` at `it·IP.dt`. So a time-independent field marched
with `(dt, χ = 2)` equals the same field marched with `(2dt, χ = 1)` state for state.

**`params.χ` is the one source of the scale.** When `efield` is a [`ScaledField`](@ref), or a
`ReducedElectricField` over one, the `χ` the `ScaledField` was built with must equal `params.χ`,
and the field and its energy are divided by the `χ²` of the step. Any other field is marched as it
is.

It returns `nothing`; the states are in `SS`, and the integration parameters hold the times.

# Throws

- `ArgumentError` when the snapshots carry another number of components than one; when
  `params.χ` is not finite and `> 0`, which refuses `0`, a negative, `Inf` and `NaN` alike, as
  [`ScaledField`](@ref) does; and when `efield` is a `ScaledField`, or a `ReducedElectricField`
  over one, built with another `χ` than `params.χ`, naming both. The counts of `IP` are checked by
  [`IntegratorParameters`](@ref) itself.
- `DimensionMismatch` when the snapshots do not hold `IP.nₛ` saved times, and when the particle
  list holds another number of particles than the snapshots do.
- The `NamedTuple`'s own error when `params` holds no `χ`. There is no default.

# Cost

The loop's two calls — [`full_leapfrog_step!`](@ref) and [`save_full_solution`](@ref) — are
inferred and allocate nothing. The buffers of the state are allocated once, before the loop.
"""
function integrate_vp(P₀::ParticleList{T},
        params::NamedTuple,
        efield::ElectricField,
        SS::Snapshots{T},
        IP::IntegratorParameters{T},
        p;
        save = true) where {T}
    _check_full_march(SS, P₀, IP, params)

    nsave = div(IP.nₜ, IP.nₛ - 1)

    # the state of the march: a copy of the list, which the step raises in place
    x = copy(vec(P₀.x))
    v = copy(vec(P₀.v))
    w = copy(vec(P₀.w))
    a = zeros(T, length(P₀))

    # effective timestep: the χ scale, converted to `T` before the product. The times of the
    # field are unscaled, so a step of `IP.dt` in `IP.t` is one step of the march.
    Δt = IP.dt * T(params.χ)

    # the field at the scale of the parameters
    field = _with_scale(efield, params.χ)

    save && save_full_solution(SS, x, v, w, field, p, zero(T), 1)

    for it in 1:IP.nₜ
        full_leapfrog_step!(x, v, a, w, field, Δt, (it - T(0.5)) * IP.dt)

        if save && it % nsave == 0
            save_full_solution(SS, x, v, w, field, p, it * IP.dt, div(it, nsave) + 1)
        end
    end

    return nothing
end
