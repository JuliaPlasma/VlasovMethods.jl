# What the three bump-on-tail scripts share: the parameters of the self-consistent set, the
# training and test κ sweeps, the initial condition, the field a march takes, the runid and the
# domain length.
#
# The field a march takes is built per sample, `vp_field(potential, params)`: the sample's scaled
# field, over the potential and at the χ of the same parameter tuple the march reads, so the scale
# reaches the march once, in its timestep and in its field.
#
# VlasovMethods and the standard library only, so that `test/integration/bump_on_tail_scripts.jl`
# includes this file without the environment of `scripts/`.

using Random
using VlasovMethods

# The parameters of the set, as the paper's scripts set them. `κ_REF` is the wave number of the
# initial condition and of the domain; a sample's χ is its κ over this, and it is the sample's time
# and field scale.
const κ_REF = 0.3
const ε_BOT = 0.03
const α_BOT = 0.1
const σ_BOT = 0.5
const v₀_BOT = 4.5
const DT_BOT = 0.1
const T_BOT = 25.0
const Nₜ_BOT = round(Int, T_BOT / DT_BOT)
const NCELLS_BOT = 16
const DEGREE_BOT = 3
const NPARTICLES_BOT = 10^4
const SEED_BOT = 1234
const NSAMPLES_BOT = 10
const RUNID_BOT = "BoT_np1e4_k_010_050_nchi_10_T25"

# The two ends of the κ range the paper sweeps; a sample's χ is its κ over `κ_REF`.
const κ_RANGE = (0.1, 0.5)

"The length of the periodic domain, ``2π/κ`` at the reference wave number."
domain_length() = 2π / κ_REF

"""
    training_chi() -> Vector{Float64}

The ten training χ of the paper's sweep: ten κ evenly spaced over `κ_RANGE`, each over `κ_REF`.
The list is explicit because `Parameter` of ReducedComplexityModeling rejects a one-sample range
whose bounds differ (K56), so the scripts pass samples rather than a bound pair everywhere.
"""
function training_chi()
    return collect(range(κ_RANGE[1], κ_RANGE[2]; length = NSAMPLES_BOT)) ./ κ_REF
end

"""
    test_chi() -> Vector{Float64}

The ten test χ: ten κ drawn uniformly from `κ_RANGE` under `Random.seed!(SEED_BOT)`, sorted, over
`κ_REF`, which is the paper's seeded test set.
"""
function test_chi()
    Random.seed!(SEED_BOT)
    κ = κ_RANGE[1] .+ (κ_RANGE[2] - κ_RANGE[1]) .* rand(NSAMPLES_BOT)
    return sort(κ) ./ κ_REF
end

"""
    initial_distribution(; n = NPARTICLES_BOT) -> ParticleDistribution

The initial condition of the set, drawn under `Random.seed!(SEED_BOT)`: the accept–reject sample
of `BumpOnTail` at the reference wave number, so one draw serves every sample of a sweep.
"""
function initial_distribution(; n = NPARTICLES_BOT)
    Random.seed!(SEED_BOT)
    dist = ParticleDistribution(1, 1, n)
    initialize!(dist, BumpOnTail(; ε = ε_BOT, κ = κ_REF, α = α_BOT, σ = σ_BOT, v₀ = v₀_BOT))
    return dist
end

"The parameters a march reads: the one χ of the sample, which the march applies itself."
march_parameters(χ) = (χ = χ,)

"""
    vp_field(potential, params) -> ScaledField

The field a march of this set takes, scaled at the sample's χ: `params` is the parameter tuple
handed to the same march, so the wrapper's χ and the march's are one value. A plain field would
leave the march a change of the time unit alone.
"""
vp_field(potential, params) = ScaledField(PoissonField(potential), params.χ)

"""
    commit_id(root = dirname(@__DIR__)) -> String

The commit of the VlasovMethods tree the scripts run in, from `git rev-parse HEAD` in `root`;
`"unknown"` where `root` is not a git checkout, so that the attribute every written file carries
exists either way.
"""
function commit_id(root = dirname(@__DIR__))
    try
        return readchomp(`git -C $root rev-parse HEAD`)
    catch error
        error isa Union{Base.IOError, ProcessFailedException} || rethrow()
        return "unknown"
    end
end
