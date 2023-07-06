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

    # projection of delta functions to splines of @jipolanco 
    # https://github.com/jipolanco/BSplineKit.jl/issues/48
    for p in eachindex(velocities)
        ilast, bs = final_dist.basis(velocities[p])  # same as `evaluate_all`
        # Iterate over evaluated basis functions.
        # The indices of the evaluated basis functions are ilast:-1:(ilast - k + 1),
        # where k is the spline order.
        for (δi, bi) ∈ pairs(bs)
            i = ilast + 1 - δi
            rhs[i] += bi * dist.particles.w[1,p]
        end
    end

    # compute and return coeffs
    ldiv!(final_dist.coefficients, final_dist.mass_fact, rhs)

    return final_dist.spline
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
