# The two Landau methods side by side at the Maxwell molecules, `γ = 0`: the error of each
# against the closed form of the rate of the second moments, and the wall-clock time of each, over
# the number of particles `N` and the number of basis functions `M`.
#
#   julia --project=scripts scripts/landau_compare.jl [--quick] [--lattice 12,20,40,60] \
#       [--ncell 4,8,15,16] [--repeats 3] [--runs-dir dir] [--results-dir dir]
#
# Each `n` of `--lattice` gives the weighted lattice of `scripts/landau_cases.jl`'s rate check with
# `N = n²` particles, on `[−3.4, 3.4]²`, weighted by the bi-Maxwellian of the variances `0.5` and
# `0.3` and normalised to one; each `ncell` of `--ncell` gives the quadratic spline on `ncell`
# cells per axis of `[−3.5, 3.5]²`, of `M = (ncell + 2)²` basis functions. For every pair
# `(N, M)`, `PairwiseLandau` and `Landau` evaluate `collisional_vectorfield!` at `γ = 0` and `ν = 1`
# on the same particles. The script measures the instantaneous `(d gap/dt) / gap`, with
# `gap = T₁₁ − T₂₂`, as `gap_rate` of `landau_cases.jl` does, its relative error against the closed
# form `−8` that `closed_gap_rate` reads from `rate_checks` of `scripts/reference_forms.jl`, and the
# smallest wall-clock time of `--repeats` evaluations after one evaluation that compiles. Each
# measured rate is printed with its sign: the spline method's vector field has the opposite sign of
# the appendix's (K96 of `KNOWN_ISSUES.md`), so its rate lies near `+8` where `N ≥ 400` and
# `M ≥ 289`, and elsewhere differs from `−8` by an amount that the table states.
#
# Every row states the rank condition of the appendix method, `dN ≥ M − 1` with `d = 2`, which
# the appendix needs for `K K⁺ = I` of the stacked matrix `K` of `M` rows and `dN` columns, and
# whether every support of a basis function holds a particle, which the appendix needs as well. A
# row where either fails is printed and drawn as it is measured.
#
# The script writes `landau_compare.md`, the table with its conditions, into `--results-dir`, the
# figures `landau_compare_pairwise.pdf` and `landau_compare_spline.pdf` there, the error of each
# row against its wall-clock time (`plot_cost`), and the record `landau_compare.hdf5` into
# `--runs-dir`.
#
# `--quick` takes the lattices `12` and `20` and the cells `8`, `15` and `16`. Its rows include
# `N = 144`, `M = 289`, where `dN = M − 1` and the condition holds at its bound, and `N = 144`,
# `M = 324`, where it fails.

include(joinpath(@__DIR__, "landau_cases.jl"))

# The half-width of the velocity domain of every spline of the sweep.
const COMPARE_DOMAIN = RATE_DOMAIN

# The dimension `d` of the velocity space in the rank condition `dN ≥ M − 1`.
const VELOCITY_DIMENSION = 2

"""
    rank_condition(npart, nbasis)

Whether `d N ≥ M − 1` holds with `d = VELOCITY_DIMENSION`, `N = npart` and `M = nbasis`: the
condition of the appendix method for the stacked matrix `K` of `M` rows and `dN` columns.
"""
rank_condition(npart, nbasis) = VELOCITY_DIMENSION * npart ≥ nbasis - 1

"""
    empty_supports(sdist, v)

The number of basis functions of `sdist` whose support holds none of the velocities `v`
(`2 × N`): those whose deposit of a unit weight from every particle is zero.
"""
function empty_supports(sdist, v)
    B = sdist.basis
    load = zeros(size(B)...)
    buffers = VlasovMethods._local_buffers(B, Float64)
    for p in axes(v, 2)
        VlasovMethods._deposit!(load, B, view(v, :, p), 1.0, buffers)
    end
    return count(iszero, load)
end

"""
    measure(method_name, n, ncell, repeats) -> NamedTuple

The row of `method_name` on the `n × n` lattice and `ncell` cells per axis: `N`, `M`, the
measured `(d gap/dt) / gap`, its closed form, the relative error, the smallest wall-clock time of
`repeats` evaluations of `collisional_vectorfield!`, the rank condition and the number of empty
supports.
"""
function measure(method_name, n, ncell, repeats)
    dist = weighted_lattice(n)
    sdist = landau_spline(ncell, COMPARE_DOMAIN)
    supports = empty_supports(sdist, dist.particles.v)
    model = METHODS[method_name].model(dist, CollisionEntropy(sdist); γ = 0, ν = 1.0)
    measured, T₁₁, T₂₂ = gap_rate(model, dist)
    closed = closed_gap_rate(T₁₁, T₂₂)
    v = copy(dist.particles.v)
    v̇ = similar(v)
    seconds = minimum(1:repeats) do _
        @elapsed VlasovMethods.collisional_vectorfield!(v̇, v, nothing, model)
    end
    npart, nbasis = n^2, length(sdist)
    return (method = method_name, npart, nbasis, measured, closed,
        error = abs(measured - closed) / abs(closed), seconds,
        rank = rank_condition(npart, nbasis), empty = supports)
end

const TABLE_HEADER = "| method | N | M | dN | dN ≥ M − 1 | particles in every support |" *
                     " (d gap/dt)/gap | closed form | relative error | wall-clock time [s] |"

function table_row(row)
    return @sprintf("| %s | %d | %d | %d | %s | %s | %+.4f | %+.4f | %.3e | %.3e |",
        row.method, row.npart, row.nbasis, VELOCITY_DIMENSION * row.npart,
        row.rank ? "holds" : "fails", row.empty == 0 ? "yes" : "no, $(row.empty) empty",
        row.measured, row.closed, row.error, row.seconds)
end

const CAPTION = "one marker per row of the table: the relative error of `(d gap/dt)/gap` " *
                "against the closed form `−8` and the wall-clock time of one evaluation of the " *
                "vector field; every row is drawn, those where the rank condition or a support " *
                "fails included."

function write_table(name, rows, quick, options)
    path = results_path(name)
    open(path, "w") do io
        println(io, "# The two Landau methods at the Maxwell molecules, `γ = 0`")
        println(io)
        println(io, "`", COMMAND, "`, commit `", COMMIT, "`, ", quick ? "quick" : "full",
            " sweep: the `n × n` lattices `n = ", join(options.lattice, ", "),
            "` on `[−", RATE_EXTENT, ", ", RATE_EXTENT, "]²`, weighted by the bi-Maxwellian ",
            "of the variances `", LATTICE_VARIANCES[1], "` and `", LATTICE_VARIANCES[2],
            "` and normalised to one; `ncell = ", join(options.ncell, ", "),
            "` quadratic cells per axis on `[−", COMPARE_DOMAIN, ", ", COMPARE_DOMAIN,
            "]²`; `ν = 1`; the smallest wall-clock time of ", options.repeats,
            " evaluations of `collisional_vectorfield!` after one that compiles.")
        println(io)
        println(
            io, "The rank condition of the appendix method is `dN ≥ M − 1`, `d = 2`, for ",
            "the stacked matrix `K` of `M` rows and `dN` columns; the appendix also needs a ",
            "particle in every support of a basis function. `Landau` solves each velocity ",
            "component alone (K96 of `KNOWN_ISSUES.md`). A row where either condition fails is ",
            "printed as it is measured.")
        println(io)
        println(io, TABLE_HEADER)
        println(io, "|:--|--:|--:|--:|:--|:--|--:|--:|--:|--:|")
        foreach(row -> println(io, table_row(row)), rows)
        println(io)
        for method_name in sort!(collect(keys(METHODS)))
            println(io, "`landau_compare_", method_name, ".pdf`, `", method_name, "`: ",
                CAPTION)
        end
    end
    @printf("  wrote %s\n", path)
    return path
end

function main(args)
    quick, options = pd_lb_options(args;
        paper = (lattice = [12, 20, 40, 60], ncell = [4, 8, 15, 16], repeats = 3),
        quick = (lattice = [12, 20], ncell = [8, 15, 16]))
    @printf("%s: %s sweep, lattices %s, cells %s, γ = 0\n",
        SCRIPT, quick ? "quick" : "full",
        join(options.lattice, ", "), join(options.ncell, ", "))

    rows = [measure(method_name, n, ncell, options.repeats)
            for method_name in sort!(collect(keys(METHODS)))
            for n in options.lattice for ncell in options.ncell]
    println(TABLE_HEADER)
    foreach(row -> println(table_row(row)), rows)

    write_table("landau_compare.md", rows, quick, options)
    column(method_name, field) = [getfield(row, field)
                                  for row in rows if row.method == method_name]
    write_record("landau_compare.hdf5", options;
        (Symbol(method_name, "_", field) => column(method_name, field)
        for method_name in keys(METHODS)
        for field in (:npart, :nbasis, :measured, :closed, :error, :seconds, :empty))...,
        (Symbol(method_name, "_rank") => Int.(column(method_name, :rank))
        for method_name in keys(METHODS))...)
    for method_name in sort!(collect(keys(METHODS)))
        save_figure("landau_compare_$(method_name).pdf",
            plot_cost(column(method_name, :seconds), column(method_name, :error)))
        println("  landau_compare_$(method_name).pdf: ", CAPTION)
    end

    report_time()
end

abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
