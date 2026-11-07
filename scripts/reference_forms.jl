# The two closed forms of the 2D Landau reference cases, derived from the moment equations and
# checked against them: the rate at which the second moments of the bi-Maxwellian isotropise at the
# Maxwell molecules of `γ = 0`, and the BKW solution of Carrillo, Hu, Wang and Wu (2020), Appendix
# A, at this package's kernel.
#
# ## The moment equations
#
# For the homogeneous Landau equation
#
#     ∂_t f = ∇ · ∫ U(v - v_*) [ f_*(v) ∇f(v) - f(v) ∇_* f_*(v_*) ] dv_*
#
# with the kernel `U(u) = |u|² (I - û ûᵀ)` of the exponent `γ = 0` — the Maxwell molecules — a test
# function `φ` gives the weak form
#
#     d/dt ∫ φ f dv = - ½ ∫∫ f(v) f(v_*) U_ij(v - v_*) (∂_i φ(v) - ∂_i φ(v_*))
#                                                  (∂_j log f(v) - ∂_j log f(v_*)) dv dv_* ,
#
# which `moment_rate` evaluates. With `φ = |v|²` the φ-difference is `2u`, and `U(u) u = 0` for the
# kernel `|u|² (I - û ûᵀ)`, so the right side vanishes pointwise: the energy — the trace
# `T₁₁ + T₂₂` — is conserved for every `f`, as an identity of the kernel and not as a property of
# the profile. No check below rests on it; what is checked is the profile's energy, which is a
# quadrature of the profile and not of the weak form.
#
# ## The rate of the second moments
#
# Take `φ = v_k v_l` and `f = f_*` the bi-Maxwellian of the density one and the unit mass with
# `Σ = diag(T₁, T₂)`. Then `∂_j log f(v) = -v_j / T_j`, `∂_i φ(v) - ∂_i φ(v_*) = δ_ik u_l + δ_il
# u_k` for `u = v - v_*`, and `u` is `N(0, 2Σ)`, so
#
#     dT_kl/dt = ½ Σ_j (1/T_j) E[ (U_kj(u) u_l + U_lj(u) u_k) u_j ] .
#
# With `U₁₁ = u₂²`, `U₁₂ = U₂₁ = -u₁u₂` and `U₂₂ = u₁²` the only moment it needs is
# `E[u₁² u₂²] = (2T₁)(2T₂) = 4T₁T₂`, so
#
#     dT₁₁/dt = (1/T₁ - 1/T₂) E[u₁²u₂²] = 4 (T₂ - T₁) ,     dT₂₂/dt = 4 (T₁ - T₂) ,
#
# the correlations `T₁₂` stay zero, and the anisotropy decays at the rate
#
#     d(T₁₁ - T₂₂)/dt = - 8 (T₁₁ - T₂₂) .
#
# The constant is the kernel's and the dimension's alone: it does not depend on the temperature,
# which is what the Maxwell molecules' kernel `|u|²` gives. It is the rate at `ν = 1`, density one
# and unit mass, and a run at another `ν` scales it. The control of the check is the same closed
# form with the constant `4` in place of `8`, which the integral must refuse.
#
# ## The BKW solution
#
# Carrillo, Hu, Wang and Wu (2020, Appendix A) give, for their kernel `A(z) = (|z|² I - z ⊗ z)/16`,
# the profile
#
#     f(v) = (2πK)^(-d/2) exp(-|v|²/(2K)) [ ((d+2)K - d)/(2K) + (1 - K)|v|²/(2K²) ] ,
#     K(t) = 1 - (1/2) exp(-2 B (d - 1) t) ,     B = 1/16 ,  d = 2 ,
#
# that is `K(t) = 1 - exp(-t/8)/2`, of mass one and temperature one. This package's kernel is
# `U = |z|² (I - û ûᵀ) = 16 A`, and the operator is linear in the kernel, so every rate is sixteen
# times the paper's: at `d = 2` their profile with
#
#     K(t) = 1 - (1/2) exp(-2t) ,     dK/dt = exp(-2t) ,
#
# solves this package's equation at `ν = 1`. At `K = 1` it is the Maxwellian of temperature one; at
# `K = 1/2` it is `exp(-|v|²) |v|²/π`, which vanishes quadratically at the origin and is why the
# checks below take their times `t₀ > 0`. Its second moment is the identity at every `t`, by the
# moment equation above with `T₁ = T₂ = 1`, and its fourth moment is `M₄(K) = 8K(2 - K)`. The fourth
# moment's equation,
#
#     dM₄/dt = 16 (1 - K) dK/dt = 32 (1 - K)² ,
#
# is what fixes the rate `dK/dt` of the profile and not only its value: a profile whose `K(t)`
# decays at another rate breaks it, while the mass and the energy do not see `K` at all. The
# equation is autonomous, so no check here fixes the origin of `t`.
#
# ## What is checked
#
# `--check` evaluates both derivations and writes nothing. The rate: `rate_checks` integrates the
# weak form for `φ = v₁²` and `φ = v₂²` over the bi-Maxwellian at the four temperature ratios
# `1.2`, `1.5`, `2` and `3`, and compares each with `4 (T₂ - T₁)` and `4 (T₁ - T₂)`, and their
# difference with `-8 (T₁ - T₂)`, to the relative tolerance `TOLERANCE`; the control's constant `4`
# must miss by more than `1e-3` relative. The BKW: `bkw_rows` takes the profile's mass, energy,
# temperature and fourth moment by quadrature at four times inside `(0, 1)` and compares them with
# `1`, `2`, `1` and `8K(2 - K)`, and compares the fourth moment's equation — the quadrature of the
# weak form — with `16 (1 - K) dK/dt`, each to `TOLERANCE` relative. `sample_rows` samples the
# profile through `BKWDistribution` at `t₀ = 0.2` and requires each of the sample's first two moments
# within four standard errors, and its fourth moment within four of the sample's own standard
# errors.
#
# A run without `--check` prints the same numbers without failing, and writes them with the command
# line into `results_path("reference_forms.txt")`. The flags are those of `options.jl`, and every
# one of them but `--check` is read in either mode, so a mistyped flag is refused instead of run at
# its default; a bare `--check` reads none of them, so that it writes nothing at all and may run
# inside a checkout.
#
# ## Running it
#
#     julia --project=scripts scripts/reference_forms.jl --check
#     julia --project=scripts scripts/reference_forms.jl --results-dir /tmp/results

include(joinpath(@__DIR__, "options.jl"))

using LinearAlgebra: SymTridiagonal, eigen
using VlasovMethods: BKWDistribution, ParticleDistribution, bkw_parameter, initialize!

"The relative tolerance of every comparison of a closed form with the moment equations."
const TOLERANCE = 1e-10

"The number of Gauss–Hermite nodes per dimension of every quadrature of the checks."
const NODES = 16

"The temperature ratios of the rate check, `T₁/T₂` from `1.2` to `3`."
const TEMPERATURES = ((1.2, 1.0), (1.5, 1.0), (2.0, 1.0), (3.0, 1.0))

"The times of the BKW check, all inside `(0, 1)`."
const TIMES = (0.1, 0.3, 0.6, 0.9)

## the quadrature

"""
    gauss_hermite(n) -> (nodes, weights)

The `n`-point Gauss–Hermite rule of `∫ exp(-x²) g(x) dx`: its nodes are the eigenvalues of the
tridiagonal matrix of the Hermite recurrence, and its weights `√π` times the squared first component
of each eigenvector. It is exact for a polynomial of degree `2n - 1`.
"""
function gauss_hermite(n)
    rule = eigen(SymTridiagonal(zeros(n), [sqrt(k / 2) for k in 1:(n - 1)]))
    return rule.values, sqrt(π) .* abs2.(rule.vectors[1, :])
end

## the kernel and the moment equations

"The three independent entries `(U₁₁, U₁₂, U₂₂)` of the kernel `U(u) = |u|² (I - û ûᵀ)` at `u`."
kernel_entries(u₁, u₂) = (u₂^2, -u₁ * u₂, u₁^2)

"""
    moment_rate(dφ, dlogf, weight, T₁, T₂; n = NODES) -> Float64

The right side of the moment equation of the test function whose gradient is `dφ`, over a
distribution of the Gaussian envelope `N(0, diag(T₁, T₂))`: the weak form of the header, which is
minus one half of the double integral of `f f_* U_ij` against the two differences of gradients.
`dlogf(v)` is the gradient of `log f` and `weight(v)` the distribution's factor beyond the envelope,
one for a Maxwellian.

The integral is that of the product Gauss–Hermite rule of the envelope: with `v = √(2T₁) x₁` and
`√(2T₂) x₂` each Gaussian contributes the rule's weight, and the pair of dimensions the factor
`1/π²`. A polynomial integrand is integrated exactly, so three nodes per dimension already suffice
for a Maxwellian.
"""
function moment_rate(dφ, dlogf, weight, T₁, T₂; n = NODES)
    nodes, weights = gauss_hermite(n)
    σ₁, σ₂ = sqrt(2 * T₁), sqrt(2 * T₂)
    total = zero(Float64)
    for a in 1:n, b in 1:n, c in 1:n, d in 1:n
        v = (σ₁ * nodes[a], σ₂ * nodes[b])
        vstar = (σ₁ * nodes[c], σ₂ * nodes[d])
        U₁₁, U₁₂, U₂₂ = kernel_entries(v[1] - vstar[1], v[2] - vstar[2])
        dv, dvs = dφ(v), dφ(vstar)
        lv, lvs = dlogf(v), dlogf(vstar)
        δ₁, δ₂ = lv[1] - lvs[1], lv[2] - lvs[2]
        total += weights[a] * weights[b] * weights[c] * weights[d] * weight(v) *
                 weight(vstar) *
                 ((dv[1] - dvs[1]) * (U₁₁ * δ₁ + U₁₂ * δ₂) +
                  (dv[2] - dvs[2]) * (U₁₂ * δ₁ + U₂₂ * δ₂))
    end
    return -total / (2 * π^2)
end

"The gradient `(-v₁/T₁, -v₂/T₂)` of the logarithm of the bi-Maxwellian of the temperatures."
maxwellian_log_gradient(T₁, T₂) = v -> (-v[1] / T₁, -v[2] / T₂)

"The weight of a Maxwellian beyond its envelope, which `moment_rate` multiplies the integrand by."
maxwellian_weight(::Any) = 1.0

## the rate of the second moments

"""
    rate_checks(T₁, T₂; n = NODES) -> NamedTuple

The rate of the second moments of the bi-Maxwellian of the temperatures `T₁` and `T₂`, as the
moment equation's integral gives it: `rate₁₁` for the test function `v₁²`, `rate₂₂` for `v₂²`, and
`gap`, their difference, which is `d(T₁₁ - T₂₂)/dt`. The closed forms are `4 (T₂ - T₁)`,
`4 (T₁ - T₂)` and `-8 (T₁ - T₂)`.
"""
function rate_checks(T₁, T₂; n = NODES)
    logf = maxwellian_log_gradient(T₁, T₂)
    rate₁₁ = moment_rate(v -> (2v[1], 0.0), logf, maxwellian_weight, T₁, T₂; n)
    rate₂₂ = moment_rate(v -> (0.0, 2v[2]), logf, maxwellian_weight, T₁, T₂; n)
    return (rate₁₁ = rate₁₁, rate₂₂ = rate₂₂, gap = rate₁₁ - rate₂₂)
end

## the BKW profile

"""
    bkw_parameter_rate(t) -> Float64

`dK/dt = exp(-2t)`, the derivative in `t` of `bkw_parameter`. It stands as its own function because
the fourth moment's equation relates the profile to this rate: a `K(t)` that decays at another rate
breaks that equation, and the rate is what the equation checks.
"""
bkw_parameter_rate(t) = exp(-2 * t)

"""
    bkw_coefficients(K) -> (A, B)

The coefficients `A = (2K - 1)/K` and `B = (1 - K)/(2K²)` of the profile's quadratic factor
`A + B |v|²`, at `d = 2`. At `K ≤ 1/2` the factor is negative near the origin and the profile is no
distribution.
"""
bkw_coefficients(K) = ((2K - 1) / K, (1 - K) / (2K^2))

"The profile's quadratic factor `A + B |v|²` at `K`, which `moment_rate` takes as the weight."
function bkw_weight(K)
    A, B = bkw_coefficients(K)
    return v -> A + B * (v[1]^2 + v[2]^2)
end

"The gradient `(∂₁ log f, ∂₂ log f)` of the logarithm of the profile of parameter `K`."
function bkw_log_gradient(K)
    A, B = bkw_coefficients(K)
    return function (v)
        quadratic = A + B * (v[1]^2 + v[2]^2)
        return (-v[1] / K + 2 * B * v[1] / quadratic, -v[2] / K + 2 * B * v[2] / quadratic)
    end
end

"""
    bkw_fourth_moment(K) -> Float64

`M₄(K) = 8K(2 - K)`, the fourth moment `∫ |v|⁴ f dv` of the profile of parameter `K`: the
Gaussian's `8` at `K = 1`, and `6` at `K = 1/2`.
"""
bkw_fourth_moment(K) = 8 * K * (2 - K)

"""
    bkw_moment(K, φ; n = NODES) -> Float64

The moment `∫ φ f dv` of the profile of parameter `K`, by the Gauss–Hermite rule of its envelope
`N(0, K I)`: with `v = √(2K) x` each Gaussian contributes the rule's weight, and the pair of
dimensions the factor `1/π`.
"""
function bkw_moment(K, φ; n = NODES)
    nodes, weights = gauss_hermite(n)
    σ = sqrt(2 * K)
    weight = bkw_weight(K)
    total = zero(Float64)
    for a in 1:n, b in 1:n

        v = (σ * nodes[a], σ * nodes[b])
        total += weights[a] * weights[b] * φ(v) * weight(v)
    end
    return total / π
end

## the checks

"The relative difference `|value - closed| / |closed|` of a value and its closed form."
relative(value, closed) = abs(value - closed) / abs(closed)

"""
The row of a check: the quantity's name, the value, the reference it is compared with — its closed
form, or for a quantity that must vanish or stay small the bound — `miss`, how far the value is from
the reference in the units of the reference, and `tolerance`, which `miss` must not exceed.
"""
function row(name, value, closed)
    (; name, value, reference = closed, miss = relative(value, closed),
        tolerance = TOLERANCE)
end
function bounded(name, value, bound)
    (; name, value, reference = bound, miss = abs(value) / bound,
        tolerance = 1.0)
end

"""
    rate_rows(T₁, T₂; n = NODES) -> Vector{NamedTuple}

The three rows of the rate check at the temperatures `T₁` and `T₂`: `dT₁₁/dt`, `dT₂₂/dt` and
`d(T₁₁ - T₂₂)/dt` as the moment equation's integral gives them, each beside its closed form.
"""
function rate_rows(T₁, T₂; n = NODES)
    r = rate_checks(T₁, T₂; n)
    ratio = round(T₁ / T₂; digits = 2)
    return [row("dT₁₁/dt at T₁/T₂ = $ratio", r.rate₁₁, 4 * (T₂ - T₁)),
        row("dT₂₂/dt at T₁/T₂ = $ratio", r.rate₂₂, 4 * (T₁ - T₂)),
        row("d(T₁₁ - T₂₂)/dt at T₁/T₂ = $ratio", r.gap, -8 * (T₁ - T₂))]
end

"""
    control_rows(T₁, T₂; n = NODES) -> Vector{NamedTuple}

The control of the rate check: the closed form of `d(T₁₁ - T₂₂)/dt` with the constant `4` in place
of `8`, and how far the integral is from it, which must exceed `1e-3` relative.
"""
function control_rows(T₁, T₂; n = NODES)
    gap = rate_checks(T₁, T₂; n).gap
    ratio = round(T₁ / T₂; digits = 2)
    return [row("d(T₁₁ - T₂₂)/dt with the constant 4, T₁/T₂ = $ratio", gap, -4 * (T₁ - T₂))]
end

"""
    bkw_rows(t₀; n = NODES) -> Vector{NamedTuple}

The rows of the BKW check at the time `t₀`: the profile's mass, energy, temperature and fourth
moment, and the equation of its fourth moment. The mass and the energy of the profile do not depend
on `K`, and its fourth moment is `8K(2 - K)`; the fourth moment's equation is the integral of the
weak form against `16 (1 - K) dK/dt`.
"""
function bkw_rows(t₀; n = NODES)
    K = bkw_parameter(t₀)
    logf, weight = bkw_log_gradient(K), bkw_weight(K)
    squared = v -> v[1]^2 + v[2]^2
    fourth_rate = moment_rate(
        v -> (4 * squared(v) * v[1], 4 * squared(v) * v[2]), logf, weight,
        K, K; n)
    energy = bkw_moment(K, squared; n)
    name = "t₀ = $t₀"
    return [row("mass, $name", bkw_moment(K, v -> 1.0; n), 1.0),
        row("energy, $name", energy, 2.0),
        row("temperature, $name", energy / 2, 1.0),
        row("fourth moment, $name", bkw_moment(K, v -> squared(v)^2; n), bkw_fourth_moment(K)),
        row("fourth moment's equation, $name", fourth_rate,
            16 * (1 - K) * bkw_parameter_rate(t₀))]
end

"""
    sample_rows(t₀; npart = 100_000, seed = 20261107) -> Vector{NamedTuple}

The rows of the sample check at the time `t₀`: the deviations of the sample's two means and two
variances from the profile's `0` and `1`, in standard errors — `σ/√N` for a mean and `σ²√(2/N)` for
a variance — and the deviation of its fourth moment in the sample's own standard error, each bounded
by four. The sample is the one `initialize!` draws for a `BKWDistribution`, so that this check covers
the type the package exports.
"""
function sample_rows(t₀; npart = 100_000, seed = 20261107)
    dist = initialize!(ParticleDistribution(1, 2, npart), BKWDistribution(t₀; seed))
    v = dist.particles.v
    w = collect(vec(dist.particles.w))
    means = vec(sum(v .* w'; dims = 2))
    centered = v .- means
    variances = [sum(w .* centered[d, :] .^ 2) for d in 1:2]
    fourth = sum(w .* (v[1, :] .^ 2 .+ v[2, :] .^ 2) .^ 2)
    fourth_se = sqrt(sum((w .* ((v[1, :] .^ 2 .+ v[2, :] .^ 2) .^ 2 .- fourth)) .^ 2))
    name = "t₀ = $t₀"
    return [bounded("mean v₁, $name", means[1] / (1 / sqrt(npart)), 4.0),
        bounded("mean v₂, $name", means[2] / (1 / sqrt(npart)), 4.0),
        bounded("variance v₁, $name", (variances[1] - 1) / sqrt(2 / npart), 4.0),
        bounded("variance v₂, $name", (variances[2] - 1) / sqrt(2 / npart), 4.0),
        bounded("fourth moment, $name",
            (fourth - bkw_fourth_moment(bkw_parameter(t₀))) / fourth_se, 4.0)]
end

"""
    check_rows(; n = NODES) -> Vector{NamedTuple}

Every row of the two derivations' checks and of the sample check, each a `(name, value, reference,
miss, tolerance)` whose `miss` must not exceed its `tolerance`.
"""
function check_rows(; n = NODES)
    rows = NamedTuple[]
    for (T₁, T₂) in TEMPERATURES
        append!(rows, rate_rows(T₁, T₂; n))
    end
    for t₀ in TIMES
        append!(rows, bkw_rows(t₀; n))
    end
    append!(rows, sample_rows(0.2))
    return rows
end

"""
    control_check(; n = NODES) -> Vector{NamedTuple}

The rows of the control: one per temperature ratio, with the relative distance of the integral from
the closed form that has the constant `4` in place of `8`. Each must exceed `1e-3`.
"""
function control_check(; n = NODES)
    rows = NamedTuple[]
    for (T₁, T₂) in TEMPERATURES
        append!(rows, control_rows(T₁, T₂; n))
    end
    return rows
end

"""
    check(; n = NODES) -> String

The table of every row of the checks and of the control, and the `ArgumentError` that names each row
whose `miss` exceeds its `tolerance` and each control row whose distance does not.
"""
function check(; n = NODES)
    rows, controls = check_rows(; n), control_check(; n)
    io = IOBuffer()
    println(
        io, "# name", " "^30, lpad("value", 26), lpad("reference", 26), lpad("miss", 12))
    for r in rows
        println(io, rpad(r.name, 36), lpad(string(round(r.value; sigdigits = 12)), 26),
            lpad(string(round(r.reference; sigdigits = 12)), 26),
            lpad(string(round(r.miss; sigdigits = 3)), 12))
    end
    for r in controls
        println(io, rpad(r.name, 36), lpad(string(round(r.value; sigdigits = 12)), 26),
            lpad(string(round(r.reference; sigdigits = 12)), 26),
            lpad(string(round(r.miss; sigdigits = 3)), 12))
    end
    failed = filter(r -> !(r.miss ≤ r.tolerance), rows)
    isempty(failed) || throw(ArgumentError(
        "the closed forms and the moment equations disagree: " *
        join((r.name * " misses by " * string(r.miss) for r in failed), ", ")))
    missed = filter(r -> !(r.miss > 1e-3), controls)
    isempty(missed) || throw(ArgumentError("the control is not refused: " *
                        join((r.name for r in missed), ", ")))
    return String(take!(io))
end

"""
    main(args = ARGS)

The run: `--check` verifies the two derivations and writes nothing, and a run without it prints the
same table and writes it, with the command line, to `results_path("reference_forms.txt")`.
"""
function main(args = ARGS)
    check_only = "--check" in args
    rest = filter(!=("--check"), args)
    # every argument but `--check` is a flag of `options.jl`, in either mode, so that a mistyped one
    # is refused instead of run at its default. A bare `--check` reads none of them, and
    # `parse_arguments` would create the two output directories, which `--check` must not.
    if !(check_only && isempty(rest))
        parse_arguments(rest)
    end
    text = check()
    print(stdout, text)
    check_only && return nothing
    path = results_path("reference_forms.txt")
    write(path, "# command: " * join([basename(PROGRAM_FILE); args], " ") * "\n" * text)
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
