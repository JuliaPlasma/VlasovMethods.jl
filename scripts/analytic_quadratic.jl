# The analytic quadratic model: the full-order solution and its two reduced forms.
#
# The model is `H(q, p) = ½|p|² + ½β²|q|²` of 1000 particles, whose solution has the closed form
# `q(t) = q₀ cos(βt) + p₀/β sin(βt)`, `p(t) = p₀ cos(βt) - β q₀ sin(βt)`. The initial condition is
# the draw of ReducedBasisMethods' `reference/Vlasov_Analytic.jl`: `q₀ = 10·randn` after the seed
# 1, `p₀ = randn` after the seed 10, and then each particle whose `rand()` after the seed 100 is at
# least `1/(1 + 0.3)` moved into the bump, `p₀ ← 0.5·p₀ + 4`. The training set is the closed form
# of the six β = 5.95, 5.97, …, 6.05 at the saved times, and the test is β = 6.
#
# The full model is the closed form at β = 6 on the step dt = 0.001, of which every `nsave`-th
# state is saved. The two reduced forms of ReducedBasisMethods march the projected initial
# condition on the step 0.0001: the ODE form in the POD basis of k = 5 with RK4, and the PODE form
# in the cotangent-lift basis of k = 20 with Lobatto IIIA–IIIB(2). Each saves `--saves` states after
# t = 0, evenly spaced to the final time. The paper command runs to T = 1000 with 10 000 saved
# states; the default runs to T = 100 with 1000:
#
#     julia --project=scripts scripts/analytic_quadratic.jl --final-time 1000 --saves 10000
#
# Writes into `runs_path`: `analytic_full.h5`, `analytic_ode_k5.h5` and `analytic_pode_k20.h5`,
# each with the saved times `t`, the positions `q` and the velocities `p` of every saved state, and
# its Hamiltonian `H`; the reduced files also hold the reduced states `z`. Each has the root
# attributes `dt` and `nsave`. Draws the two spectra of the training set into `results_path`.

include("options.jl")

using CairoMakie
# an explicit only-list, and not `using GeometricIntegrators`: both it and VlasovMethods export
# `GeometricIntegrator`
using GeometricIntegrators: LobattoIIIAIIIB, RK4, integrate
# `attributes` by name: CairoMakie exports it too
using HDF5: attributes, h5open
using Random
using ReducedBasisMethods: ReducedBasis, PODSVD, CotangentLiftSVD, lift, project,
                           reduced_ode, reduced_pode, spectrumplot
using ReducedComplexityModeling: Parameter, ParameterSpace
using VlasovMethods: AnalyticHamiltonian

const OPT = parse_arguments(ARGS; final_time = 100.0, saves = 1000, particles = 1000)

const TEST_β = 6.0
const TRAINING_β = [5.95, 5.97, 5.99, 6.01, 6.03, 6.05]
const DT_FULL = 0.001
const DT_REDUCED = 0.0001

"The initial condition of the reference: positions, then velocities with a bump at `v₀ = 4`."
function initial_condition(n)
    Random.seed!(1)
    q₀ = 10 .* randn(n)
    Random.seed!(10)
    p₀ = randn(n)
    Random.seed!(100)
    for i in 1:n
        rand() ≥ 1 / (1 + 0.3) && (p₀[i] = 0.5 * p₀[i] + 4.0)
    end
    return q₀, p₀
end

"The closed form at `β`, as the stacked `[q; p]` of one column per time of `t`."
function solution(β, q₀, p₀, t)
    return [q₀ .* cos.(β .* t') .+ (p₀ ./ β) .* sin.(β .* t');
            p₀ .* cos.(β .* t') .- (β .* q₀) .* sin.(β .* t')]
end

"The number of steps of `dt` between two saves, which must be a whole number."
function steps_per_save(o, dt)
    nsave = round(Int, o.final_time / (o.saves * dt))
    isapprox(nsave * o.saves * dt, o.final_time; rtol = 1e-12) || throw(ArgumentError(
        "the final time $(o.final_time) is not $(o.saves) saves of a whole number of steps " *
        "of $dt"))
    return nsave
end

"""
The march of a reduced form in `saves` chunks of `nsave` steps, of which the last state of each is
saved: `step(z)` integrates one chunk from `z` and returns its last state. The model is autonomous,
so every chunk restarts the clock at zero.
"""
function march(step, z₀, saves)
    Z = zeros(length(z₀), saves + 1)
    Z[:, 1] .= z₀
    for j in 1:saves
        Z[:, j + 1] .= step(Z[:, j])
    end
    return Z
end

"Write one result: the times, the lifted states, the Hamiltonian, and any other dataset."
function write_result(file, t, x, h; dt, nsave, datasets...)
    n = size(x, 1) ÷ 2
    h5open(file, "w") do io
        io["t"] = collect(t)
        io["q"] = x[1:n, :]
        io["p"] = x[(n + 1):end, :]
        io["H"] = [h(0.0, c, NamedTuple()) for c in eachcol(x)]
        for (name, value) in pairs(datasets)
            io[string(name)] = value
        end
        attributes(io)["dt"] = dt
        attributes(io)["nsave"] = nsave
    end
    println("wrote ", file)
end

function main(o)
    t = range(0, o.final_time; length = o.saves + 1)
    q₀, p₀ = initial_condition(o.particles)
    h = AnalyticHamiltonian(TEST_β)

    # the full model, saved every `nsave` steps of `DT_FULL`
    nsave_full = steps_per_save(o, DT_FULL)
    nsave_full ≥ 10 || throw(ArgumentError(
        "the full model saves every $(nsave_full) steps; it saves every 10th step at most"))
    write_result(runs_path("analytic_full.h5"), t, solution(TEST_β, q₀, p₀, t), h;
        dt = DT_FULL, nsave = nsave_full)

    # the bases of the training set at the saved times
    S = reduce(hcat, [solution(β, q₀, p₀, t) for β in TRAINING_β])
    space = ParameterSpace(Parameter(:β, TRAINING_β))
    pod = ReducedBasis(PODSVD(), S, NamedTuple(), space; k = 5)
    cotangent = ReducedBasis(CotangentLiftSVD(), S, NamedTuple(), space; k = 20)
    save(results_path("item11_analytic_pod_spectrum.png"),
        spectrumplot(pod.Λ; title = "analytic, POD", ylabel = "σᵢ²"))
    save(results_path("item11_analytic_cotangent_spectrum.png"),
        spectrumplot(cotangent.Λ; title = "analytic, cotangent lift", ylabel = "σᵢ²"))

    # the two reduced forms, saved every `nsave` steps of `DT_REDUCED`
    nsave = steps_per_save(o, DT_REDUCED)
    span = (0.0, nsave * DT_REDUCED)
    x₀ = [q₀; p₀]

    ode(z) = integrate(reduced_ode(pod, h, span, DT_REDUCED, z), RK4()).q[end]
    Z = march(ode, project(pod, x₀), o.saves)
    write_result(runs_path("analytic_ode_k5.h5"), t, lift(pod, Z), h;
        dt = DT_REDUCED, nsave = nsave, z = Z)

    function pode(z)
        sol = integrate(
            reduced_pode(cotangent, h, span, DT_REDUCED, z), LobattoIIIAIIIB(2);
            f_abstol = 1e-14)
        return vcat(sol.q[end], sol.p[end])
    end
    Z = march(pode, project(cotangent, x₀), o.saves)
    write_result(runs_path("analytic_pode_k20.h5"), t, lift(cotangent, Z), h;
        dt = DT_REDUCED, nsave = nsave, z = Z)
    return nothing
end

main(OPT)
