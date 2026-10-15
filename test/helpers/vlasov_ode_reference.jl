# The reference of check 4: the Galerkin form of the POD reduced equation, transcribed from
# `reference/Reduced_Vlasov_ODE.jl:64-84` (the plan-wide rule of 2026-10-10 puts a reference
# generator in `test/helpers/`). It is listed in no `runtests.jl`.
#
# The reference writes the RHS as the reference does: `x = U*z`, the second block of `x` is `ẋ`, and
# the first is negated by `E(x, β) = β²x`, so that `ż = UᵀJ_N∇H(Uz)`.

"""
    vlasov_ode_reference!(ż, t, z, β, U)

The reference's `RHS`: `x = U*z`, and then `ż = Uᵀ[x_{N+1:2N}; -β²x_{1:N}]` — the position block of
`U` multiplies the velocity half of the lifted state, and its velocity block multiplies the negated
position half — with `N = length(x)/2` and `U` the `2N × k` POD basis. `t` is ignored, as the
reference ignores it.
"""
function vlasov_ode_reference!(zdot, t, z, β, U)
    x = U * z
    N = div(length(x), 2)
    n = length(z)

    @inbounds for i in 1:n
        tmp = zero(eltype(zdot))
        for k in 1:N
            tmp += U[k, i] * x[N + k]
        end
        for k in 1:N
            tmp -= β^2 * U[N + k, i] * x[k]
        end
        zdot[i] = tmp
    end

    return zdot
end
