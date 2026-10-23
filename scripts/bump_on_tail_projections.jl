# Bump-on-tail, step 2 of 3: the reduced bases of the training sweep, and the files of the HNN.
#
# Reads the snapshots that `bump_on_tail_training.jl` wrote with the same options, and builds the
# cotangent-lift basis of the particles from the eigendecomposition of the Gram matrix of their
# positions and velocities, and the POD basis of the field from that of their accelerations. Each
# basis keeps the leading vectors that carry all but its tolerance of the energy of the spectrum.
# The defaults are the paper's tolerances, 1e-9 for the particles and 1e-5 for the field, which
# `VP-ROM-Projections.ipynb` of ReducedBasisMethods uses. The paper command is the default:
#
#     julia --project=scripts scripts/bump_on_tail_projections.jl --particle-tolerance 1e-9 --field-tolerance 1e-5
#
# Writes `runs_path("<run id>_projections.h5")`, the particle basis with the field basis in its
# group `field`, and `runs_path("<run id>_vectorfield.h5")`, the reduced positions `X` and
# accelerations `A` of every saved state, each of size `(k, nₛ, nparam)`. Both carry the root
# attributes `commit`, `command` and `seed`. Draws the two spectra into `results_path`.

include("options.jl")
include("bump_on_tail_common.jl")

using CairoMakie
using ReducedBasisMethods: CotangentLiftEVD, PODEVD, spectrumplot

const OPT = parse_arguments(ARGS; BUMP_ON_TAIL..., particle_tolerance = 1e-9,
    field_tolerance = 1e-5)

function main(o)
    id = run_id(o)
    SS = h5open(io -> Snapshots(io), runs_path(id * ".h5"), "r")

    particles = particle_basis(CotangentLiftEVD(), SS; tolerance = o.particle_tolerance)
    field = field_basis(PODEVD(), SS; tolerance = o.field_tolerance)
    println("k = $(particles.k), kₑ = $(field.k)")

    X, A = reduced_data(particles, SS)
    write_p39_files(runs_path(id * "_projections.h5"), runs_path(id * "_vectorfield.h5"),
        particles, X, A; commit = git_commit(), command = join(ARGS, " "), seed = o.seed,
        field_basis = field)

    save(results_path("item11_particle_spectrum.png"),
        spectrumplot(particles.Λ; title = "particles, k = $(particles.k)", ylabel = "λᵢ"))
    save(results_path("item11_field_spectrum.png"),
        spectrumplot(field.Λ; title = "field, kₑ = $(field.k)", ylabel = "λᵢ"))
    return nothing
end

main(OPT)
