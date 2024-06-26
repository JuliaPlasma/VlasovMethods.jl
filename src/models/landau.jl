struct Landau{XD, VD, DT <: DistributionFunction{XD,VD}, ET <: Entropy, T} <: VlasovModel
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 
    
    function Landau(dist::DistributionFunction{XD,VD}, ent::Entropy; ν::T=1.) where {XD, VD, T}
        new{XD, VD, typeof(dist), typeof(ent), T}(dist, ent, ν)
    end
end

function compute_J(sdist::SplineDistribution{1,2})
    int = zeros(size(sdist))
    B = sdist.basis
    d_start = BSplineKit.knots(B)[1] 
    d_end = BSplineKit.knots(B)[end]

    for k in 1:size(sdist)
        i, j = ij_from_k(k, length(B))
        integrand(v) = B[i](v[1]) * B[j](v[2]) * (1. + log(abs(sdist.spline(v)))) 
        int[k], _ = hcubature(integrand,[d_start, d_start], [d_end, d_end], atol = 1e-4, rtol=1e-4)
    end

    ldiv!(sdist.mass_fact, int)

    return int
end

function compute_U(v_α, v_β)
    n = length(v_α)
    U = zeros(eltype(v_α),(n,n))

    norm_diff = norm(v_α - v_β)

    if v_α != v_β
        for i in CartesianIndices(U)
            if i[1] == i[2]
                U[i] += 1/norm_diff
            end
            U[i] -= (v_α[i[1]] - v_β[i[1]])*(v_α[i[2]] - v_β[i[2]])./norm_diff^3
        end
    end
    return U
end

function Landau_rhs!(v̇, v, v_array, L, B, sdist)
    # computes rhs for a single particle, assuming that the projection and other particle velocities are taken from the previous timestep 
    # params.L is the vector L_k, which depends on the projection
    # params.idist.particles.v is the particle velocities
    # params.fdist.basis is the spline basis

    K = size(sdist) # length of tensor product basis (this is the square of the 1d basis length)

    for α in axes(v_array, 2)
        U = compute_U(v, v_array[:,α])
        for k in 1:K # could make this more precise by using find_knot_interval function 
            v̇ .+=  L[k] * U * (eval_bfd(B, k, v) - eval_bfd(B, k, v_array[:, α]))
        end
    end

    return v̇
end

# particle-to-particle version
function Landau_rhs(v, params)
    # computes rhs for a single particle, assuming that the projection and other particle velocities are taken from the previous timestep 
    # params.L is the vector L_k, which depends on the projection
    # params.idist.particles.v is the particle velocities
    # params.fdist.basis is the spline basis
    v̇ = zero(v)
    K = size(params.sdist) # length of tensor product basis (this is the square of the 1d basis length)

    ind1, res1 = evaluate_der_2d(params.B, v)

    for α in axes(params.v_array, 2)
        U = compute_U(v, params.v_array[:,α])
        ind_α, res_α = evaluate_der_2d(params.B, params.v_array[:, α])

        for (i, k) in pairs(ind1)
            if k > 0 && k ≤ K
                v̇ .+=  params.dist.particles.w[1,α] * params.L[k] * U * res1[:, i]
            end
        end

        for (i2, k2) in pairs(ind_α)
            if k2 > 0 && k2 ≤ K
                v̇ .-=  params.dist.particles.w[1,α] * params.L[k2] * U * res_α[:, i2]
            end
        end
    end

    return v̇
end


function compute_K_plus(v_array, dist, sdist)
    M = size(sdist)
    K1 = zeros(M, size(v_array,2))
    K2 = zeros(M, size(v_array,2))

    for α in axes(v_array, 2)
        klist, der_array = evaluate_der_2d(sdist.basis, v_array[:,α])
        for (i, k) in pairs(klist)
            if k > 0 && k <= M
                K1[k, α] = dist.particles.w[1,α] * der_array[1,i]
                K2[k, α] = dist.particles.w[1,α] * der_array[2,i]
            end
        end
    end

    if rank(K1) < M  || rank(K2) < M
        println("K1 or K2 not full rank")
        @show size(K1,1) - rank(K1)
        @show size(K2,1) - rank(K2)
    end

    return pinv(K1), pinv(K2)
end

function f_Maxwellian(v)
    return 1/(2π) * exp(- norm(v)^2 / 2)
end

function L_integrand(v, k, sdist)
    v1 = [v[1], v[2]]
    v2 = [v[3], v[4]]
    
    id_list_1, _ = evaluate_der_2d(sdist.basis, v1)
    id_list_2, _ = evaluate_der_2d(sdist.basis, v2)
    
    if (k[1] in id_list_1 || k[1] in id_list_2) && (k[2] in id_list_1 || k[2] in id_list_2)

        integrand = transpose(eval_bfd(sdist.basis, k[1],v1) .- eval_bfd(sdist.basis, k[1], v2)) * sdist.spline(v1) * compute_U(v1, v2)
        return integrand * sdist.spline(v2) * (eval_bfd(sdist.basis, k[2],v1) .- eval_bfd(sdist.basis, k[2], v2))

    else
        return zero(eltype(v))
    end
end

function compute_L_ij(sdist)
    L = zeros(size(sdist), size(sdist))
    d_start = BSplineKit.knots(sdist.basis)[1] 
    d_end = BSplineKit.knots(sdist.basis)[end]
    
    Threads.@threads for k in CartesianIndices(L)
        integrand(v) = L_integrand(v, k, sdist)
        L[k], _ = hcubature(integrand, [d_start, d_start, d_start, d_start], [d_end, d_end, d_end, d_end], atol = 1e-4, rtol = 1e-4)
    end
    
    return L .* 0.5
end

# spline-to-spline? version 
function Landau_rhs_2!(v̇, t, v, params)
    # v̇ = zero(v)
    # project v onto params.sdist
    S = projection(v, params.dist, params.sdist)

    # compute K matrices 
    K1_plus, K2_plus = compute_K_plus(v, params.dist, params.sdist)

    # compute L_ij matrix
    Lij = compute_L_ij(params.sdist)

    # compute J vector
    J = compute_J(params.sdist)

    # solve for vector field
    v̇[1,:] .= K1_plus * Lij * J  
    v̇[2,:] .= K2_plus * Lij * J  

    return v̇

end
