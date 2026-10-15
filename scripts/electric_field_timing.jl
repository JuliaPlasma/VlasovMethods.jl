# The A/B timing of the reduced electric field: the non-hyper-reduced path of
# `ReducedElectricField` against the DEIM path of `DEIMElectricField`.
#
# Run it in a fresh process, without code-coverage instrumentation:
#
#     julia --startup-file=no --project=scripts scripts/electric_field_timing.jl
#
# A path named as the argument, `reduced` or `deim`, times that path alone and prints its
# median. With no argument each path runs in a cold subprocess of its own, and the script exits
# with status 1 unless the DEIM median is below the reduced-path median. The sizes are the run
# of `bump_on_tail_3_testing.jl` of `ReducedBasisMethods`: 5·10⁴ particles, a reduced dimension
# of 20, and as many DEIM interpolation points. The timing runs in `Float64`.
#
# The timed call is one field evaluation of a march, `f(E, X, w, t)`: the deposit of the lifted
# positions, the solve of the potential, and the field evaluation, which the DEIM path performs
# at the interpolation points rather than at every particle. Both paths deposit every particle,
# so the deposit is not what separates them.

using LinearAlgebra
using PoissonSolvers
using Random
using ReducedBasisMethods: get_DEIM_indices
using Statistics: median
using VlasovMethods: DEIMElectricField, PoissonField, ReducedElectricField

const N = 50_000      # particles
const K = 20          # reduced dimension
const Kₑ = 20         # DEIM interpolation points
const RUNS = 5        # runs per path
const CALLS = 200     # field evaluations per run

const PATHS = ("reduced", "deim")

LinearAlgebra.BLAS.set_num_threads(1)

# The two paths over one full-order field, with the reduced coordinates of the same positions.
function fixtures(::Type{T}) where {T}
    Random.seed!(1234)
    q = T(0.2) .+ T(0.6) .* rand(T, N)
    w = fill(T(1) / N, N)
    field = PoissonField(Potential(PeriodicBasisSpline((zero(T), one(T)), 3, 16)))

    Random.seed!(42)
    Ψₓ = Matrix(qr(randn(T, N, K)).Q)
    Random.seed!(7)
    Ψₑ = Matrix(qr(randn(T, N, Kₑ)).Q)
    p = get_DEIM_indices(Ψₑ)

    X = Ψₓ' * q
    return (reduced = ReducedElectricField(field, Ψₓ, X),
        deim = DEIMElectricField(field, Ψₓ, Ψₑ, p, X), X = X, E = similar(X), w = w)
end

# The median of `RUNS` runs of `CALLS` field evaluations, in milliseconds per evaluation.
function measure(path::AbstractString)
    path in PATHS || throw(ArgumentError("the path is one of $PATHS, got $(repr(path))"))
    f = fixtures(Float64)
    field = getfield(f, Symbol(path))
    E, X, w = f.E, f.X, f.w
    t = 0.0

    # the first call compiles and warms; each run times what follows
    field(E, X, w, t)

    times = map(1:RUNS) do _
        t₀ = time_ns()
        for _ in 1:CALLS
            field(E, X, w, t)
        end
        return (time_ns() - t₀) / CALLS / 1e6
    end

    return median(times)
end

if isempty(ARGS)
    # every path in a cold process of its own, so that neither warms the other
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    project = @__DIR__
    script = @__FILE__

    medians = map(PATHS) do path
        out = read(`$julia --startup-file=no --project=$project $script $path`, String)
        return parse(Float64, strip(out))
    end

    for (path, m) in zip(PATHS, medians)
        println(
            rpad(path, 8), " median ", m, " ms per field evaluation of ", N, " particles")
    end
    exit(medians[2] < medians[1] ? 0 : 1)
else
    # one path: the value on stdout is all this process prints, so that the caller parses it
    println(measure(ARGS[1]))
end
