# Q: dispatch here on the entropy or entropy.dist? 
"""
Assumes final_dist contains basis for projection, and projects init_dist onto it.
Updates final_dist with projected Spline and coefficients.

"""
# function projection!(init_dist::ParticleDistribution, final_dist::SplineDistribution)

#     B = final_dist.basis
#     rhs = zero(final_dist.coefficients)
#     Mfact = cholesky!(galerkin_matrix(B))


#     # projection of delta functions to splines of @jipolanco 
#     # https://github.com/jipolanco/BSplineKit.jl/issues/48
#     for p in init_dist.particles
#         ilast, bs = B(p.v[1])  # same as `evaluate_all`
#         # Iterate over evaluated basis functions.
#         # The indices of the evaluated basis functions are ilast:-1:(ilast - k + 1),
#         # where k is the spline order.
#         for (δi, bi) ∈ pairs(bs)
#             i = ilast + 1 - δi
#             rhs[i] += bi * p.w[1]
#         end
#     end

#     # compute and store coefficients
#     ldiv!(final_dist.coefficients, Mfact, rhs)

#     # compute spline from projection coeffs
#     final_dist.spline = Spline(B, final_dist.coefficients)
# end


function projection(velocities::AbstractArray{VT}, dist::ParticleDistribution, final_dist::SplineDistribution) where {VT}
    rhs = zeros(VT, size(final_dist))
    coeffs = copy(rhs)
    Mfact = cholesky!(galerkin_matrix(final_dist.basis))
    # Mfact = cholesky!(final_dist.mass_matrix)

    # projection of delta functions to splines of @jipolanco 
    # https://github.com/jipolanco/BSplineKit.jl/issues/48
    for i in 1:length(velocities)
        ilast, bs = final_dist.basis(velocities[i])  # same as `evaluate_all`
        # Iterate over evaluated basis functions.
        # The indices of the evaluated basis functions are ilast:-1:(ilast - k + 1),
        # where k is the spline order.
        for (δi, bi) ∈ pairs(bs)
            i = ilast + 1 - δi
            rhs[i] += bi * dist.particles.w[1,i]
        end
    end

    # old implementation
    # for v in velocities
    #     for i in eachindex(basis)
    #         rhs[i] += basis[i,VT](v) * weight[1]
    #     end
    # end

    # compute and return coeffs
    ldiv!(coeffs, Mfact, rhs)

    return Spline(final_dist.basis, coeffs)
end

# TODO 
function projection!(init_dist::SplineDistribution, final_dist::ParticleDistribution)


end


function projection!(init_dist::ParticleDistribution, final_dist::ParticleDistribution)
    final_dist = init_dist
end


function projection!(init_dist::SplineDistribution, final_dist::SplineDistribution)
    final_dist = init_dist
end
