# The parts that the three bump-on-tail drivers share: their options and run id, the domain, the
# draw of the particles, the field of a march, the sweeps, the accelerations of a training sweep,
# the reduced data and the writer of the two files that the HNN reads, and the growth-rate fit of
# the reference.
#
# The bump-on-tail problem at the wave number `κ` runs on the domain of the reference wave number
# `κ₀ = 0.3`, rescaled by `χ = κ/κ₀`, so that one draw of the particles serves every sample of a
# sweep. A march takes the scale from `params.χ`, for the timestep and for the field, so every march
# takes the field of `vp_field` with the same `params`.
#
# This file loads no plotting package, so `test/integration/bump_on_tail_scripts.jl` includes it,
# and every package it loads is in `test/Project.toml` as well as in `scripts/Project.toml`.

# `attributes` by name: Makie exports it too, and a driver that loads CairoMakie after this file
# would find the bare name ambiguous
using HDF5
using HDF5: attributes
using PoissonSolvers
using Random
using ReducedBasisMethods: ReducedBasis, h5save
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods
using VlasovMethods: efield!, update!

"The reference wave number `κ₀`: the domain is `[0, 2π/κ₀]`, and a sample at `κ` has `χ = κ/κ₀`."
const KAPPA_REFERENCE = 0.3

"""
The options that the three bump-on-tail drivers share, with the values of ReducedBasisMethods'
`scripts/bump_on_tail_1_training.jl`: the timestep `dt`, the final time, the number of particles,
the number of spline cells of the potential, the spline degree and the seed of the draw. Every
driver takes all of them, so that the three agree on the run id.
"""
const BUMP_ON_TAIL = (
    dt = 0.1, final_time = 25.0, particles = 10_000, cells = 16, degree = 3, seed = 1234)

"""
    run_id(o)

The name of the files of one bump-on-tail run, from the options `o` that set its training
snapshots.
"""
function run_id(o)
    return "bump_on_tail_np$(o.particles)_nh$(o.cells)_p$(o.degree)_dt$(o.dt)_" *
           "T$(o.final_time)_seed$(o.seed)"
end

"""
    step_count(final_time, dt)

The number of steps of `dt` to `final_time`, rounded to the nearest integer: `step_count(25.0,
0.1) == 250`, where `Int(div(25.0, 0.1))` gives `249`.
"""
step_count(final_time, dt) = round(Int, final_time / dt)

"""
    training_kappas()

The ten wave numbers of the training sweep, evenly spaced in `[0.1, 0.5]`: as wave numbers, the
`χ` range `0.1/κ₀ … 0.5/κ₀` of ten samples in a commented-out `Parameter` line of
ReducedBasisMethods' `scripts/bump_on_tail_1_training.jl`, whose active line runs the one sample
`χ = 1`. The sweep passes them to `Parameter` as a list of values, as it passes the test sweep.
"""
training_kappas() = collect(range(0.1, 0.5; length = 10))

"""
    test_sweep(seed; n = 10) -> (κ, χ)

The `n` wave numbers of the test sweep, drawn uniform in `[0.1, 0.5]` after `Random.seed!(seed)`
and sorted, and their scales `χ = κ/κ₀`.
"""
function test_sweep(seed; n = 10)
    Random.seed!(seed)
    κ = sort(0.1 .+ 0.4 .* rand(n))
    return κ, κ ./ KAPPA_REFERENCE
end

"""
    parameter_space(κ)

The `ParameterSpace` of a sweep over the wave numbers `κ`: one parameter, `χ = κ/κ₀`. The space
holds its samples in increasing order, so the sample `p` is `κ[p]` of a sorted `κ`, as both sweeps
are.
"""
parameter_space(κ) = ParameterSpace(Parameter(:χ, collect(κ ./ KAPPA_REFERENCE)))

"""
    potential(o)

The periodic spline potential of the options `o`, on `o.cells` cells of degree `o.degree` of the
domain `[0, 2π/κ₀]`. `PeriodicBasisSpline` takes the order, which is the degree plus one.
"""
function potential(o)
    return Potential(PeriodicBasisSpline((0.0, 2π / KAPPA_REFERENCE), o.degree + 1, o.cells))
end

"""
    initial_particles(o)

The bump-on-tail particles of the options `o`: `o.particles` particles drawn by accept–reject at
`κ₀` after `Random.seed!(o.seed)`. Every sample of a sweep starts from them.
"""
function initial_particles(o)
    Random.seed!(o.seed)
    dist = ParticleDistribution(1, 1, o.particles)
    initialize!(dist, BumpOnTail(κ = KAPPA_REFERENCE))
    return dist.particles
end

"""
    vp_field(potential, params)

The field of every bump-on-tail march: the Poisson field of `potential`, scaled by `params.χ`.
A march takes it with the same `params`, so the field and the timestep share one `χ`.
"""
vp_field(potential, params) = ScaledField(PoissonField(potential), params.χ)

"""
    accelerations!(SS, potential, w)

Write the accelerations of each saved state of the snapshots `SS` into `SS.A`: the field of the
march, `vp_field(potential, SS.paramspace(p))`, deposited from the saved positions with the
weights `w` and evaluated at them. A march saves the positions and the velocities, and not the
accelerations, which the field basis reduces.
"""
function accelerations!(SS::Snapshots, potential, w)
    a = zeros(eltype(SS.A), size(SS.A, 2))
    for p in axes(SS.A, 4)
        field = vp_field(potential, SS.paramspace(p))
        for ts in axes(SS.A, 3)
            x = @view SS.X[1, :, ts, p]
            update!(field, x, w, zero(eltype(SS.A)))
            SS.A[1, :, ts, p] .= efield!(field, a, x)
        end
    end
    return SS
end

"""
    reduced_data(basis, SS) -> (X, A)

The reduced positions and accelerations of every saved state of the snapshots `SS` in the
position basis `basis.Φ`, each of size `(k, nₛ, nparam)`: `X[:, ts, p] = Φᵀ SS.X[1, :, ts, p]`, and
`A` the same of `SS.A`.
"""
function reduced_data(basis::ReducedBasis, SS::Snapshots)
    Φ = basis.Φ
    X = zeros(size(Φ, 2), size(SS.X, 3), size(SS.X, 4))
    A = zero(X)
    for p in axes(SS.X, 4), ts in axes(SS.X, 3)

        X[:, ts, p] .= Φ' * SS.X[1, :, ts, p]
        A[:, ts, p] .= Φ' * SS.A[1, :, ts, p]
    end
    return X, A
end

"""
    hnn_data(X, A, Λ) -> (data, target)

The training data of the HNN from the reduced positions `X`, the reduced accelerations `A` and the
leading values `Λ` of the particle spectrum: each column of `X` divided by the square root of its
value of `Λ`, and `A` as it is. A network trained on the two maps a reduced position to its
acceleration, and the division gives the coefficients one scale.
"""
hnn_data(X::AbstractMatrix, A::AbstractMatrix, Λ::AbstractVector) = (X ./ sqrt.(Λ), A)

"""
    write_p39_files(projections, vectorfield, basis, X, A; commit, command, seed,
                    field_basis = nothing)

Write the two files that the HNN reads: the particle basis `basis` into `projections`, with the
field basis `field_basis` in its group `field` when one is given, and the reduced data `X` and `A`
of [`reduced_data`](@ref) into `vectorfield`. Each file carries the root attributes `commit`,
`command` and `seed`: the commit of VlasovMethods, the command of the driver and the seed of the
draw.
"""
function write_p39_files(projections, vectorfield, basis, X, A; commit, command, seed,
        field_basis = nothing)
    size(A) == size(X) || throw(DimensionMismatch(
        "the reduced accelerations are of size $(size(A)); the positions of size $(size(X))"))
    size(X, 1) == basis.k || throw(DimensionMismatch(
        "the reduced data hold $(size(X, 1)) coefficients; the basis holds k = $(basis.k)"))
    h5open(projections, "w") do io
        h5save(io, basis)
        field_basis === nothing || h5save(io, field_basis; path = "field")
        write_provenance(io, commit, command, seed)
    end
    h5open(vectorfield, "w") do io
        io["X"] = X
        io["A"] = A
        write_provenance(io, commit, command, seed)
    end
    return projections, vectorfield
end

"""
    write_provenance(io, commit, command, seed)

Write the root attributes `commit`, `command` and `seed` of the open HDF5 file `io`.
"""
function write_provenance(io, commit, command, seed)
    attributes(io)["commit"] = commit
    attributes(io)["command"] = command
    attributes(io)["seed"] = seed
    return io
end

"""
    git_commit(root = <the checkout that holds this file>)

The commit of the git checkout at `root`, `git rev-parse HEAD`. A directory that is no checkout
throws, so no file records a commit it cannot name.
"""
git_commit(root = dirname(@__DIR__)) = readchomp(`git -C $root rev-parse HEAD`)

"""
    command_line()

The command of the running driver: its file name and its arguments, so that the default run, which
has no arguments, still names its script.
"""
command_line() = join([basename(PROGRAM_FILE); ARGS], " ")

"""
    fit_growth_rate(t, W)

The growth rate of the field energy `W` at the times `t`, with the window and the count of the
reference fit, `VP-ROM-Plotting.ipynb` cell 36: the maxima at a window of 5 samples on each side,
and the first four of them, or all when there are fewer. A maximum is the one of `growth_rate`,
which is not the notebook's rule of five monotone samples on each side. It is
`growth_rate(t, W; window = 5, maxima = 1:m)`, half the slope of `log W`, and `NaN` where fewer
than two maxima leave no line to fit. A matrix `W` gives one rate per column.
"""
function fit_growth_rate(t, W::AbstractVector)
    # The maxima of `growth_rate` itself, so that the count and the fit find the same points. The
    # call of the internal `_find_maxima` is deliberate: no public function gives the count.
    m = min(4, length(VlasovMethods._find_maxima(W, 5)))
    m ≥ 2 || return NaN
    return growth_rate(t, W; window = 5, maxima = 1:m)
end

fit_growth_rate(t, W::AbstractMatrix) = [fit_growth_rate(t, w) for w in eachcol(W)]
