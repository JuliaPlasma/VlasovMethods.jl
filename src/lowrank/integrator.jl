"""
    lowrank_step!(Y, method, field, t, τ)

Advance the [`LowRankDistribution`](@ref) `Y` by one step of length `τ` from the time `t` of the
dynamical low-rank `method`, in place, and return `Y`. `field` is the right-hand side of the
method. The methods are [`BUG`](@ref), [`AugmentedBUG`](@ref) and [`MidpointBUG`](@ref), and
[`ProjectorSplitting`](@ref); the docstring of each method names its algorithm and the fields
that it takes.

`t` and `τ` are converted to the element type `T` of `Y` before their first use. A `τ` that is not
finite in `T` throws an `ArgumentError` before any write, whose message prints `τ` as it was
given: on a `Float32` state, the finite `τ = 1e300` reads `got τ = 1.0e300`. The step assigns new
factors to `Y`, possibly of another rank, and sets `nfixed = 0`. The entries of `Y` are not
checked.

[`LowRankIntegrator`](@ref) and [`run!`](@ref) take a run of steps.
"""
function lowrank_step! end

# τ in the element type T of the state, refused where it is not finite in T
function _finite_step(τ, ::Type{T}) where {T}
    τT = T(τ)
    isfinite(τT) || throw(ArgumentError(
        "the step τ is not finite in the element type $T of the state; got τ = $τ"))
    return τT
end

# One step of the classical RK4 tableau, b = (1, 2, 2, 1)/6, of ẏ = rhs(y, s) over [t, t + τ],
# with `rhs!(ẏ, y, s)` writing ẏ and returning it: every substep of the low-rank steps.
function _rk4(rhs!, y, t, τ)
    k1 = rhs!(similar(y), y, t)
    k2 = rhs!(similar(y), y .+ (τ / 2) .* k1, t + τ / 2)
    k3 = rhs!(similar(y), y .+ (τ / 2) .* k2, t + τ / 2)
    k4 = rhs!(similar(y), y .+ τ .* k3, t + τ)
    return y .+ (τ / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    LowRankIntegrator(Y, method, field, tspan, tstep)

The run of a dynamical low-rank `method` from the [`LowRankDistribution`](@ref) `Y` over the
span `tspan = (t₀, t₁)` with the step `tstep`, for a method with a [`lowrank_step!`](@ref) on
the workspace `field` that keeps the rank, or for `AugmentedBUG(tol; rmin, rmax)`, whose rank
stays within `rmax`; [`run!`](@ref) refuses a step that grows the rank beyond. The fields `Y`,
`method`, `field`, `tspan` and `tstep` hold them, the times promoted to one type. [`run!`](@ref)
advances `Y`.

Throws an `ArgumentError` for a step that is not finite and positive, for a span that is not
finite or has `t₁ < t₀`, and for a span whose grid `t₀:tstep:t₁` does not increase strictly, as
the `Float32` grid `16777216f0:1f0:16777218f0`, which repeats its first time.
"""
struct LowRankIntegrator{YT <: LowRankDistribution, MT, FT, TT <: Real} <: VlasovMethod
    Y::YT
    method::MT
    field::FT
    tspan::Tuple{TT, TT}
    tstep::TT

    function LowRankIntegrator(Y::YT, method::MT, field::FT, tspan::Tuple{Real, Real},
            tstep::Real) where {YT <: LowRankDistribution, MT, FT}
        t₀, t₁, Δt = promote(tspan..., tstep)
        isfinite(Δt) && Δt > 0 || throw(ArgumentError(
            "LowRankIntegrator takes a finite and positive step; got $tstep"))
        isfinite(t₀) && isfinite(t₁) && t₀ ≤ t₁ || throw(ArgumentError(
            "LowRankIntegrator takes a finite span (t₀, t₁) with t₀ ≤ t₁; got $tspan"))
        issorted(t₀:Δt:t₁; lt = <=) || throw(ArgumentError(
            "LowRankIntegrator takes a span whose times t₀:tstep:t₁ increase strictly; got " *
            "a repeated time for $tspan and the step $tstep"))
        new{YT, MT, FT, typeof(Δt)}(Y, method, field, (t₀, t₁), Δt)
    end
end

# the columns of the stored factors: the largest rank that the method can reach
_stored_rank(method, Y) = size(Y.S, 1)

# a method without rank bounds takes the rank of every state
_check_rank_bounds(method, r) = nothing

# the step of a run, which `lowrank_step!` converts to the element type of Y
function _check_run_step(Y, tstep)
    h = eltype(Y)(tstep)
    isfinite(h) && h > 0 || throw(ArgumentError(
        "the step $tstep is zero or not finite in the element type $(eltype(Y)) of the state"))
    return nothing
end

# the element type T, refused before the file opens where HDF5 has no datatype for it
function _check_hdf5_type(::Type{T}) where {T}
    storable = try
        HDF5.datatype(T)
        true
    catch
        false
    end
    storable || throw(ArgumentError(
        "run! writes the times and the state to HDF5, which has no datatype for $T"))
    return nothing
end

"""
    run!(int::LowRankIntegrator, h5file)

Advance the state `Y` of `int` over the times `t₀:tstep:t₁`, the initial time included, in place
with `lowrank_step!(Y, method, field, tₙ, tstep)`, write every state into the HDF5 file
`h5file`, and return `Y`. The file is created, or overwritten, and closed on return, also when a
step throws. Before the file is opened, `run!` throws an `ArgumentError` for a state whose rank
is outside `rmin:rmax` of an `AugmentedBUG(tol; rmin, rmax)`, for a `tstep` that is zero or not
finite in the element type of `Y`, and for an element type of the times or of `Y` that HDF5 does
not store, such as `Int128`, `Float16` or `BigFloat`; an existing file is then left unchanged.

Time `t₀` is the initial state, written before any step, and step `s` is the `s`-th call of
`lowrank_step!`, which writes the time `t₀ + s · tstep`. The file holds the datasets `t` of the
`nt` times, `rank` of the `nt` ranks, and `X` of size `(nx, R, nt)`, `S` of size `(R, R, nt)`
and `V` of size `(nv, R, nt)`, zero beyond the rank of each state, with `R` the rank of `Y` or,
for `AugmentedBUG(tol; rmin, rmax)`, `min(rmax, nx, nv)`.

A step that leaves `all(isfinite, Y.S)` false stops the run without an error: nothing is written
for it, its `rank` and factors stay zero, and the file attribute `status` is
`"failed at step <s>"`. A step that leaves a rank `r > R` stops the run the same way, with
`status` set to `"failed at step <s>: rank <r> beyond the <R> stored columns"`, and then throws
an `ArgumentError`. A step that throws sets `status` to `"failed at step <s>: <E>"`, with `E`
the name of the type of the exception, and the exception propagates from `run!`. A run that
completes sets `status` to `"ok"`.
"""
function run!(int::LowRankIntegrator, h5file)
    (; Y, method, field, tstep) = int
    t₀, t₁ = int.tspan
    times = t₀:tstep:t₁
    nt = length(times)
    nx, nv = size(Y)
    _check_rank_bounds(method, size(Y.S, 1))
    _check_run_step(Y, tstep)
    R = _stored_rank(method, Y)
    DT = eltype(Y)
    _check_hdf5_type(eltype(times))
    _check_hdf5_type(DT)

    h5open(h5file, "w") do h5
        h5t = create_dataset(h5, "t", eltype(times), (nt,); chunk = (1,))
        h5r = create_dataset(h5, "rank", Int, (nt,); chunk = (1,))
        h5X = create_dataset(h5, "X", DT, (nx, R, nt); chunk = (nx, R, 1))
        h5S = create_dataset(h5, "S", DT, (R, R, nt); chunk = (R, R, 1))
        h5V = create_dataset(h5, "V", DT, (nv, R, nt); chunk = (nv, R, 1))

        # the factors of Y in the first columns of R, zero beyond
        function store!(n)
            X, S, V = zeros(DT, nx, R), zeros(DT, R, R), zeros(DT, nv, R)
            X[:, axes(Y.X, 2)] .= Y.X
            S[axes(Y.S)...] .= Y.S
            V[:, axes(Y.V, 2)] .= Y.V
            h5t[n], h5r[n] = times[n], size(Y.S, 1)
            h5X[:, :, n], h5S[:, :, n], h5V[:, :, n] = X, S, V
            return nothing
        end

        store!(1)
        status, err = "ok", nothing
        for s in 1:(nt - 1)
            try
                lowrank_step!(Y, method, field, times[s], tstep)
            catch e
                attrs(h5)["status"] = "failed at step $s: $(nameof(typeof(e)))"
                rethrow()
            end
            if !all(isfinite, Y.S)
                status = "failed at step $s"
                break
            end
            r = size(Y.S, 1)
            if r > R
                status = "failed at step $s: rank $r beyond the $R stored columns"
                err = ArgumentError("step $s of the run changes the rank to $r, beyond the \
                    $R columns that run! stores; run! takes a method that keeps the rank, or \
                    AugmentedBUG(tol; rmin, rmax)")
                break
            end
            store!(s + 1)
        end
        attrs(h5)["status"] = status
        isnothing(err) || throw(err)
    end

    return Y
end
