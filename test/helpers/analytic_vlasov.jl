# The analytic quadratic-potential model of §4 and the closed form of
# `reference/Vlasov_Analytic.jl:61-62`. The check-1 file of the reduced march and the check-3–5
# file of the reduced forms both include it; it is listed in no `runtests.jl`.
#
# The draws are `Xoshiro`, as §5V.P37 states: `Xoshiro(seed)` and `Random.seed!(seed)` followed by
# `randn(n)` produce the same values, which is asserted in the check-1 file.

using LinearAlgebra: hcat, vcat
using Random
using ReducedComplexityModeling: ParameterSpace, Parameter

"The test `β` of §4: `6.00`, off-sample, where the training six run over `[5.95, 6.05]`."
const ANALYTIC_β = 6.0

"The six training `β` of §4, in the order the training snapshot is assembled in."
const TRAINING_β = (5.95, 5.97, 5.99, 6.01, 6.03, 6.05)

"The six-sample training `ParameterSpace`, whose `β` bounds are of the element type `T`."
function training_parameter_space(::Type{T}) where {T}
    ParameterSpace(Parameter(:β, T(5.95), T(6.05), length(TRAINING_β)))
end

"The one-sample `ParameterSpace` of the test run, whose `β` is the test `β`."
test_parameter_space(::Type{T}) where {T} = ParameterSpace(Parameter(:β, T(6), T(6), 1))

"""
    analytic_initial_condition(T; n = 1000) -> (q₀, p₀)

The initial condition of §4, as `reference/Vlasov_Analytic.jl:31-47` draws it: `q₀ = 10·randn(n)`
from `Xoshiro(1)`; `p₀ ~ N(0, 1)` from `Xoshiro(10)`, then, for a particle whose seeded draw of
`Xoshiro(100)` is `≥ 1/(1 + a)`, `p₀ ← σᵥ·p₀ + v₀` with the bump's `a = 0.3`, `σᵥ = 0.5` and
`v₀ = 4.0`. The training sweep and the test share it, as the reference does.
"""
function analytic_initial_condition(::Type{T}; n = 1000) where {T}
    q0 = T(10) .* randn(Xoshiro(1), T, n)
    p0 = randn(Xoshiro(10), T, n)

    a, σv, v0 = T(0.3), T(0.5), T(4)
    rng = Xoshiro(100)
    for i in 1:n
        rand(rng) >= 1 / (1 + a) && (p0[i] = σv * p0[i] + v0)
    end

    return q0, p0
end

"""
    analytic_solution(T, β, q₀, p₀, ts) -> (Q, P)

The closed form of `reference/Vlasov_Analytic.jl:61-62` at the times `ts`, as two `n × length(ts)`
matrices:

```math
q(t) = q_0 \\cos(βt) + \\frac{p_0}{β} \\sin(βt), \\qquad
p(t) = p_0 \\cos(βt) - β q_0 \\sin(βt).
```

It has no discretisation error, which is what makes it the reference of check 1.
"""
function analytic_solution(::Type{T}, β, q0, p0, ts) where {T}
    Q = hcat([(p0 ./ β) .* sin(β * t) .+ q0 .* cos(β * t) for t in ts]...)
    P = hcat([p0 .* cos(β * t) .- (β .* q0) .* sin(β * t) for t in ts]...)
    return Q, P
end

"""
    training_snapshots(T, q₀, p₀; step = 1) -> S

The stacked `2n × 6·m` training snapshot `[Q; P]` of the six training `β`, each trajectory saved
every `step` over `[0, 1000]` — `m = div(1000, step) + 1` times each, whose `t = 0` column is the
same for all six. The paper's spacing of `0.01` would be `9.6` GB (§5V.P29); §5V.P37's edge fixes
`step = 1`, which is 1001 times and 6006 columns.

`step` is a keyword because a core check's basis is a fixture: check 4's POD basis is built from the
training set of §4's `n = 200` at `step = 10`, which is `606` columns, where the `T = 1000` runs and
check 5's basis — the cotangent-lift SVD at §4's `n = 1000` — use the default `step = 1` and its
`6006` columns.
"""
function training_snapshots(::Type{T}, q0, p0; step = 1) where {T}
    n = length(q0)
    ts = T.(0:step:1000)
    S = Matrix{T}(undef, 2n, length(TRAINING_β) * length(ts))

    for (j, βi) in enumerate(TRAINING_β)
        Q, P = analytic_solution(T, T(βi), q0, p0, ts)
        range_j = ((j - 1) * length(ts) + 1):(j * length(ts))
        S[1:n, range_j] .= Q
        S[(n + 1):(2n), range_j] .= P
    end

    return S
end
