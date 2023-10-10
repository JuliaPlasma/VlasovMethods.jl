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


function projection(velocities::AbstractVector{VT}, dist::ParticleDistribution{1,1}, final_dist::SplineDistribution{1,1}) where {VT}
    rhs = zeros(VT, size(final_dist))

    # projection of delta functions to splines of @jipolanco 
    # https://github.com/jipolanco/BSplineKit.jl/issues/48
    for p in eachindex(velocities)
        ilast, bs = final_dist.basis(velocities[p])  # same as `evaluate_all`
        # Iterate over evaluated basis functions.
        # The indices of the evaluated basis functions are ilast:-1:(ilast - k + 1),
        # where k is the spline order.
        if ilast > 0 && ilast <= size(final_dist)[1]
            for (δi, bi) ∈ pairs(bs)
                i = ilast + 1 - δi
                if i < 1 || i > 41
                    println("i=", i)
                    println("ilast = ", ilast)
                    println("bs = ", bs)
                else
                    rhs[i] += bi * dist.particles.w[1,p]
                end
            end
        end
    end

    # compute and return coeffs
    ldiv!(final_dist.coefficients, final_dist.mass_fact, rhs)
    # final_dist.coefficients ./= 4.
    # final_dist.coefficients ./= BSplineKit.knots(final_dist.basis)[end] - BSplineKit.knots(final_dist.basis)[1]

    return final_dist.spline
end


function projection(velocities::AbstractMatrix{VT}, dist::ParticleDistribution{1,2}, final_dist::SplineDistribution{1,2}) where {VT}
    rhs = zeros(VT, size(final_dist))

    # projection of delta functions to splines of @jipolanco 
    # https://github.com/jipolanco/BSplineKit.jl/issues/48
    for p in eachindex(velocities[1,:])
        ilast1, bs1 = final_dist.basis(velocities[1,p])  # same as `evaluate_all`, first component
        ilast2, bs2 = final_dist.basis(velocities[2,p])  # same as `evaluate_all`, second component
        # Iterate over evaluated basis functions.
        # The indices of the evaluated basis functions are ilast:-1:(ilast - k + 1),
        # where k is the spline order.
        for (δi, bi) ∈ pairs(bs1)
            for (δi2, b2) ∈ pairs(bs2)
                i = ilast1 + 1 - δi
                j = ilast2 + 1 - δi2
                if i > 0 && i <= size(final_dist)[1] && j > 0 && j <= size(final_dist)[2]
                    rhs[i,j] += bi * b2 * dist.particles.w[1,p]
                elseif ilast1 <= 0 || ilast1 > size(final_dist)[1] || ilast2 <= 0 || ilast2 > size(final_dist)[2]
                    println("WARNING: particle outside domain")
                end
            end
        end
    end

    # compute and return coeffs
    # ldiv!(final_dist.coefficients, final_dist.mass_fact, rhs)
    invm = inv(final_dist.mass_fact)
    final_dist.coefficients .= invm * rhs * invm

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
