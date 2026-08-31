struct LandauCache{T, PT <: ParticleDistribution, ST <: SplineDistribution{T}} <: Cache{T}
    pdist::PT
    sdist::ST

    v::Matrix{T}
    v̇::Matrix{T}

    J::Vector{T}
    L::Matrix{T}
    LJ::Vector{T}
    K1::Matrix{T}
    K2::Matrix{T}

    function LandauCache{T}(pdist, sdist) where {T}
        M = length(sdist)
        N = size(pdist.particles.v, 2)

        v = zero(pdist.particles.v)
        v̇ = zero(pdist.particles.v)

        J = zeros(T, M)
        L = zeros(T, (M, M))
        LJ = zeros(T, M)
        K1 = zeros(T, (M, N))
        K2 = zeros(T, (M, N))

        new{T, typeof(pdist), typeof(sdist)}(pdist, sdist, v, v̇, J, L, LJ, K1, K2)
    end
end

function LandauCache(pdist::ParticleDistribution{T}, sdist::SplineDistribution{T}) where {T}
    LandauCache{T}(pdist, sdist)
end

function Cache(AT, c::LandauCache{DT, PT, ST}) where {DT, PT, ST}
    LandauCache{AT}(c.pdist, similar(AT, c.sdist))
end
function CacheType(AT, c::LandauCache{DT, PT, ST}) where {DT, PT, ST}
    LandauCache{AT, PT, similar_type(AT, c.sdist)}
end

struct Landau{
    D, XD, VD, DT <: DistributionFunction{XD, VD}, ET <: Entropy, CT <: CacheDict} <:
       VlasovModel
    dist::DT    # distribution function
    entropy::ET # entropy
    ν::D        # collision frequency

    cache::CT

    function Landau(dist::DistributionFunction{XD, VD}, ent::Entropy; ν::D = 1.0) where {
            D, XD, VD}
        cache = CacheDict(LandauCache(dist, ent.dist))
        new{D, XD, VD, typeof(dist), typeof(ent), typeof(cache)}(dist, ent, ν, cache)
    end
end

Cache(AT, l::Landau) = LandauCache{AT}(l.pdist, similar(AT, l.sdist))
CacheType(AT, l::Landau) = LandauCache{AT, typeof(l.pdist), similar_type(AT, l.sdist)}

function integrand_J(v::AbstractArray{T}, B, i, j, sdist::SplineDistribution) where {T}
    B[i, T](v[1]) * B[j, T](v[2]) * (one(T) + log(sdist.spline(v)))
end

function integrand_J(v::AbstractArray{T}, params) where {T}
    integrand_J(v, params.B, params.i, params.j, params.sdist)
end

# function compute_J(sdist::SplineDistribution{1,2})
#     T = eltype(sdist)
#     int = zeros(T, length(sdist))
#     B = sdist.basis
#     d_start = T(BSplineKit.knots(B)[1])
#     d_end = T(BSplineKit.knots(B)[end])
#     domain = ([d_start, d_start], [d_end, d_end])

#     for k in 1:length(sdist)
#         i, j = ij_from_k(k, length(B))
#         integrand(v) = integrand_J(v, B, i, j, sdist)
#         # integand(v) = B[i, T](v[1]) * B[j, T](v[2]) * (1. + log(abs(sdist.spline(v))))
#         int[k], _ = HCubature.hcubature(integrand,[d_start, d_start], [d_end, d_end], atol = 1e-6, rtol=1e-6)
#     end

#     ldiv!(sdist.mass_fact, int)

#     return int
# end

function compute_J!(J, sdist::SplineDistribution{T, 1, 2}, n, ::Landau) where {T}
    for k in 1:length(sdist)
        i, j = ij_from_k(k, length(sdist.basis))
        params = (sdist = sdist, B = sdist.basis, i = i, j = j)
        J[k] = gauss_quad_2d(integrand_J, sdist.basis, n, params)
    end

    ldiv!(sdist.mass_fact, J)

    return J
end

function kernel(v_α::AbstractVector{T}, v_β::AbstractVector{T}, ::Landau) where {T}
    norm_diff = euclidean(v_α, v_β)

    if v_α != v_β
        U11 = - (v_α[1] - v_β[1]) * (v_α[1] - v_β[1]) / norm_diff^3 + inv(norm_diff)
        U12 = - (v_α[1] - v_β[1]) * (v_α[2] - v_β[2]) / norm_diff^3
        U21 = - (v_α[2] - v_β[2]) * (v_α[1] - v_β[1]) / norm_diff^3
        U22 = - (v_α[2] - v_β[2]) * (v_α[2] - v_β[2]) / norm_diff^3 + inv(norm_diff)

        return @SMatrix [U11 U12;
                         U21 U22]
    else
        return @SMatrix [zero(T) zero(T);
                         zero(T) zero(T)]
    end
end

# particle-to-particle version
# function Landau_rhs(v, params)
#     # computes rhs for a single particle, assuming that the projection and other particle velocities are taken from the previous timestep 
#     # params.L is the vector L_k, which depends on the projection
#     # params.idist.particles.v is the particle velocities
#     # params.fdist.basis is the spline basis
#     v̇ = zero(v)
#     K = size(params.sdist) # length of tensor product basis (this is the square of the 1d basis length)
#     n = length(v)
#     U = zeros(n,n)

#     ind_s, res_s = evaluate_der_2d(params.B, v)

#     for α in axes(params.v_array, 2)
#         vα = params.v_array[:,α]
#         compute_U!(U, v, vα)
#         ind_α, res_α = evaluate_der_2d(params.B, vα)

#         for (i, k) in pairs(ind_s)
#             if k > 0 && k ≤ K
#                 v̇ .+= params.dist.particles.w[1,α] * params.L[k] * (U * res_s[:, i])
#             end
#         end

#         for (i, k) in pairs(ind_α)
#             if k > 0 && k ≤ K
#                 v̇ .-= params.dist.particles.w[1,α] * params.L[k] * (U * res_α[:, i])
#             end
#         end
#     end

#     return v̇
# end

function compute_K!(K1, K2, v_array::AbstractArray{T}, sdist, landau::Landau) where {T}
    for α in axes(v_array, 2)
        klist, der_array = evaluate_der_2d(sdist.basis, v_array[:, α])
        for (i, k) in pairs(klist)
            if k > 0 && k <= length(sdist)
                K1[k, α] = landau.dist.particles.w[1, α] * der_array[1, i]
                K2[k, α] = landau.dist.particles.w[1, α] * der_array[2, i]
            end
        end
    end

    return K1, K2
end

# function compute_K_plus(v_array::AbstractArray{T}, dist, sdist) where {T}
#     K1, K2 = compute_K(v_array, dist, sdist)

#     if rank(K1) < length(sdist) || rank(K2) < length(sdist)
#         println("K1 or K2 not full rank")
#         @show size(K1,1) - rank(K1)
#         @show size(K2,1) - rank(K2)
#     end

#     return pinv(K1), pinv(K2)
# end

# function L_integrand_vec(v::AbstractVector{T}, params) where T
#     v1 = [v[1], v[2]]
#     v2 = [v[3], v[4]]

#     id_list_1 = evaluate_der_2d_indices(params.sdist.basis, v1)
#     id_list_2 = evaluate_der_2d_indices(params.sdist.basis, v2)

#     if (params.k[1] in id_list_1 || params.k[1] in id_list_2) && (params.k[2] in id_list_1 || params.k[2] in id_list_2)
#         # U = compute_U(v1, v2)
#         basis_derivative = zeros(T, 2)
#         eval_bfd!(basis_derivative, params.sdist.basis, params.k[1], v1, 0, 1)
#         eval_bfd!(basis_derivative, params.sdist.basis, params.k[1], v2, 1, -1)

#         integrand = transpose(basis_derivative) * params.sdist.spline(v1) * compute_U(v1, v2)
#         # integrand = transpose(basis_derivative) * params.sdist.spline(v1) * compute_U(v1, v2)

#         eval_bfd!(basis_derivative, params.sdist.basis, params.k[2], v1, 0, 1)
#         eval_bfd!(basis_derivative, params.sdist.basis, params.k[2], v2, 1, -1)

#         return (integrand * params.sdist.spline(v2) * basis_derivative)::T

#     else
#         return zero(T)
#     end
# end

# function compute_L_ij(sdist)
#     T = eltype(sdist)
#     L = zeros(T, (length(sdist), length(sdist)))
#     B = sdist.basis
#     M = length(B)
#     knots = T.(BSplineKit.knots(sdist.basis))
#     # d_start = T(BSplineKit.knots(sdist.basis)[1])
#     # d_end = T(BSplineKit.knots(sdist.basis)[end])
#     # domain = ([d_start, d_start, d_start, d_start], [d_end, d_end, d_end, d_end])

#     Threads.@threads for k in CartesianIndices(L)
#         i1, j1 = ij_from_k(k[1], M)
#         i2, j2 = ij_from_k(k[2], M)
#         if k[1] ≤ k[2] && length(BSplines.common_support(B[i1], B[i2])) > 1 && length(BSplines.common_support(B[j1], B[j2])) > 1

#             irange = knots[BSplines.common_support(B[i1], B[i2])]
#             jrange = knots[BSplines.common_support(B[j1], B[j2])]

#             domain = ([irange[1], jrange[1], irange[1], jrange[1]], [irange[end], jrange[end], irange[end], jrange[end]])

#             params = (k = k, sdist = sdist)
#             prob = IntegralProblem(L_integrand_vec, domain, params)
#             sol = Integrals.solve(prob, HCubatureJL(); abstol = 1e-5, reltol = 1e-5)
#             L[k] = sol.u
#         end
#     end

#     Threads.@threads for k in CartesianIndices(L)
#         if k[1] > k[2] 
#             L[k] = L[k[2], k[1]]
#         end
#     end

#     return L .* 0.5
# end

function L_integrand(
        v1::AbstractVector{T}, v2::AbstractVector{T}, sdist, i, j, landau::Landau) where {T}
    basis_derivative1 = eval_bfd(sdist.basis, i, v1) - eval_bfd(sdist.basis, i, v2)
    basis_derivative2 = eval_bfd(sdist.basis, j, v1) - eval_bfd(sdist.basis, j, v2)

    sdist.spline(v1) * dot(basis_derivative1, kernel(v1, v2, landau) * basis_derivative2) *
    sdist.spline(v2)
end

function L_integrand(v1::AbstractVector{T}, v2::AbstractVector{T}, params, landau::Landau) where {T}
    id_list_1 = evaluate_der_2d_indices(params.sdist.basis, v1)
    id_list_2 = evaluate_der_2d_indices(params.sdist.basis, v2)

    for i in eachindex(params.k)
        params.k[i] in id_list_1 || params.k[i] in id_list_2 || return zero(T)
    end

    L_integrand(v1, v2, params.sdist, params.k[1], params.k[2], landau)
end

function compute_L!(L, sdist::SplineDistribution{T, 1, 2}, n::Int, landau::Landau) where {T}
    integrand = (v1, v2, params) -> L_integrand(v1, v2, params, landau)

    for i in axes(L, 1)
        for j in axes(L, 2)[i:end]
            # i1, j1 = ij_from_k(i, M)
            # i2, j2 = ij_from_k(j, M)

            # iknots = BSplines.common_support(B[i1], B[i2])
            # jknots = BSplines.common_support(B[j1], B[j2])

            params = (k = (i, j), sdist = sdist)
            L[i, j] = gauss_quad(integrand, sdist.basis, n, params) / 2
        end
    end

    for i in axes(L, 1)
        for j in axes(L, 2)[begin:(i - 1)]
            L[i, j] = L[j, i]
        end
    end

    return L
end

# spline-to-spline? version 
function collisional_vectorfield!(v̇::AbstractArray{ST}, v::AbstractArray{ST}, params, landau::Landau) where {ST}
    cache = landau.cache[ST]

    # project v onto params.sdist
    # println("sdist")
    sdist = cache.sdist

    # println("projection")
    projection(v, landau.dist, sdist)

    # compute K matrices 
    # println("compute K")
    # K1_plus, K2_plus = compute_K_plus(v, params.dist, sdist)
    compute_K!(cache.K1, cache.K2, v, sdist, landau)

    if rank(cache.K1) < length(sdist) || rank(cache.K2) < length(sdist)
        println("K1 or K2 not full rank")
        @show size(cache.K1, 1) - rank(cache.K1)
        @show size(cache.K2, 1) - rank(cache.K2)
    end

    # compute L_ij matrix
    # println("computing L")
    compute_L!(cache.L, sdist, params.n, landau)

    # compute J vector
    # println("computing J")
    compute_J!(cache.J, sdist, params.n, landau)

    # solve for vector field
    mul!(cache.LJ, cache.L, cache.J)

    v̇[1, :] .= cache.K1 \ cache.LJ
    v̇[2, :] .= cache.K2 \ cache.LJ

    # v̇[1,:] .*= -1
    # v̇[2,:] .*= -1

    # return v̇
    return nothing
end
