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

Base.eltype(::LandauCache{T}) where {T} = T

# `l.dist` and `l.entropy.dist`: the fields are named `dist` and `entropy`, so `l.pdist` and
# `l.sdist` would have thrown. Unreachable, `CacheDict`'s parent being a `LandauCache`.
Cache(AT, l::Landau) = LandauCache{AT}(l.dist, similar(AT, l.entropy.dist))
function CacheType(AT, l::Landau)
    LandauCache{AT, typeof(l.dist), similar_type(AT, l.entropy.dist)}
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

@doc raw"""
    compute_J!(J, sdist, ::Landau)

The vector

```math
\mathbb{L}_k = \sum_i \mathbb{M}^{-1}_{ik} \int_\Omega \varphi_i(v)
               \left[ 1 + \log \Big( \sum_j f_j \varphi_j(v) \Big) \right] dv
```

of the Landau manuscript's `eq:defn_Lk` — the ``L^2`` projection of ``1 + \log f_s`` onto the
spline basis.

Written as a projection rather than as one quadrature per coefficient. The integrand is
sampled once on the tensor-product Gauß-Legendre grid the basis already carries and contracted
against the one-dimensional tabulations, so the cost is the grid size rather than the grid size
times the number of coefficients. The earlier version called a `gauss_quad_2d` that looped over
**every** cell pair for **each** of the ``M`` coefficients, and took the quadrature domain from
a hard-coded `knots[3:(end-1)]` slice — asymmetric, independent of the spline order and of the
boundary condition, and a different subrange from the one the mass matrix was assembled on, so
the two halves of `eq:defn_Lk` were not computed over the same domain.
"""
function compute_J!(J, sdist::SplineDistribution{T, XD, VD}, ::Landau) where {T, XD, VD}
    q = sdist.quadrature
    fs = sdist.spline
    X = quadrature_nodes(q)

    g = Array{T, VD}(undef, quadrature_grid_size(q)...)
    for (r, pt) in enumerate(Iterators.product(X...))
        f = fs(pt)
        # The manuscript's own appendix opens with a commented-out section titled
        # "Positivity-preserving projections", stating that f_s "has to be non-negative
        # everywhere", and does not solve it. Writing this as 0.5*log(f^2) — which an earlier
        # version of the metriplectic operator did — returns log|f| and hides the violation
        # behind a finite number, while the H-theorem it is used to prove reverses wherever
        # f_s < 0.
        f > 0 || throw(ErrorException(
            "the projected distribution is non-positive, f_s = $(f) at v = $(pt), so " *
            "log f_s in eq:defn_Lk is undefined there. An L² projection of a particle " *
            "distribution undershoots where the sampling is thin; use more particles or a " *
            "coarser spline basis."))
        g[r] = 1 + log(f)
    end

    l2_projection!(reshape(J, size(sdist)), q, g)
    return J
end

@doc raw"""
    kernel(v_α, v_β, ::Landau)

The Landau collision kernel

```math
U(u) = \frac{1}{|u|} \left( \mathbb{I} - \frac{u \otimes u}{|u|^2} \right) ,
\qquad u = v_\alpha - v_\beta ,
```

the orthogonal projector onto the complement of ``u``, scaled by ``1/|u|``.

This is the standard Landau tensor and it is **not** what the manuscript's `eq:landau_kernel`
as typeset says: that carries a prefactor ``1/|u|^3`` in front of the bracketed projector, i.e.
an extra ``|u|^{-2}``. The intended reading is
``|u|^{-3} ( |u|^2 \mathbb{I} - u \otimes u )``, which is this. The properties the
conservation proofs rest on hold for the expression implemented here: ``U`` is symmetric,
``U(v_\alpha, v_\beta) = U(v_\beta, v_\alpha)``, and ``U u = 0`` exactly.

!!! warning "The coincident-point value is a regularisation, not a limit"
    ``U`` diverges as ``|u| \to 0``; there is no value to return there. Zero is returned so
    that a product quadrature whose two factors share nodes does not produce `Inf`, but the
    error this introduces does not vanish under refinement — and because both factors of the
    double integral in [`compute_L!`](@ref) use the *same* Gauß-Legendre nodes, it fires on
    every diagonal cell pair rather than on a set of measure zero. The kernel is integrable in
    two dimensions, so the integral itself is finite; what is needed is a quadrature that
    knows about the singularity, or offset grids for the two factors. Neither is implemented.
"""
function kernel(v_α::AbstractVector{T}, v_β::AbstractVector{T}, ::Landau) where {T}
    u1 = v_α[1] - v_β[1]
    u2 = v_α[2] - v_β[2]
    n2 = u1 * u1 + u2 * u2

    iszero(n2) && return @SMatrix [zero(T) zero(T); zero(T) zero(T)]

    nd = sqrt(n2)
    inv_n = inv(nd)
    inv_n3 = inv_n / n2

    U11 = inv_n - u1 * u1 * inv_n3
    U12 = -u1 * u2 * inv_n3
    U22 = inv_n - u2 * u2 * inv_n3

    return @SMatrix [U11 U12;
                     U12 U22]
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

@doc raw"""
    compute_K!(K1, K2, v_array, sdist, landau)

The matrix ``\mathbb{K}_{k,\alpha} = w_\alpha \, \nabla \varphi_k(v_\alpha)`` of the Landau
manuscript's appendix, split into its two velocity components.

!!! note "Both matrices are cleared first"
    They are cache arrays, allocated once and written only at the entries whose basis functions
    the current particle positions overlap. Without the `fill!` the entries written on earlier
    calls survive: as particles move between cells, ``\mathbb{K}`` accumulates spurious
    nonzeros from previous Picard iterations and previous time steps, and
    ``\mathbb{K}^{+}`` is then the pseudo-inverse of a matrix that is not
    ``w_\alpha \nabla \varphi_k(v_\alpha)``.

Every index is bounded **per component** before being flattened. Testing only the flat index
does not bound the components: with `i = 0`, `j = 3`, `M = 10` the flat index
`(j-1)M + i = 20` passes a `1 ≤ k ≤ M²` test and decodes to `(10, 2)`, so an out-of-support
contribution is aliased onto an unrelated basis function instead of being discarded. Wrapping
goes through `basis_index`, which is periodic-aware, so the gradient is now evaluated with the
same wrap as the spline it differentiates — the two disagreed before.
"""
function compute_K!(K1, K2, v_array::AbstractArray{T}, sdist, landau::Landau) where {T}
    fill!(K1, zero(eltype(K1)))
    fill!(K2, zero(eltype(K2)))

    B = sdist.basis
    bs = bases(B)
    lin = LinearIndices(B)
    w = landau.dist.particles.w

    vals = ntuple(k -> zeros(T, local_width(bs[k])), 2)
    ders = ntuple(k -> zeros(T, local_width(bs[k])), 2)

    for α in axes(v_array, 2)
        v = view(v_array, :, α)
        wα = w[1, α]

        j₀ = ntuple(k -> evaluate_all!(vals[k], bs[k], v[k], 0), 2)
        ntuple(k -> evaluate_all!(ders[k], bs[k], v[k], 1), 2)

        for t1 in eachindex(vals[1]), t2 in eachindex(vals[2])

            i = basis_index(bs[1], j₀[1] + t1 - 1)
            j = basis_index(bs[2], j₀[2] + t2 - 1)
            (1 ≤ i ≤ size(B, 1) && 1 ≤ j ≤ size(B, 2)) || continue

            k = lin[i, j]
            K1[k, α] += wα * ders[1][t1] * vals[2][t2]
            K2[k, α] += wα * vals[1][t1] * ders[2][t2]
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

@doc raw"""
    gradient_tabulations(sdist)

The sparse tables ``D^c_{k,a} = \partial_c \varphi_k(v_a)`` of the two velocity components of
the basis gradient at every point of the tensor-product quadrature grid, as an `M × Q` pair.

Built as Kronecker products of the one-dimensional tabulations the quadrature already holds:
``\partial_1 (\varphi_i \otimes \varphi_j) = \varphi_i' \otimes \varphi_j``, so
`D¹ = kron(Φ₂, Φ₁')` and `D² = kron(Φ₂', Φ₁)` — with the row and column orderings that
`LinearIndices` and the grid enumeration already use. Nothing is evaluated that the quadrature
had not already tabulated.
"""
function gradient_tabulations(sdist::SplineDistribution{T, XD, 2}) where {T, XD}
    q1, q2 = quadratures(sdist.quadrature)
    Φ1, dΦ1 = basis_values(q1, 0), basis_values(q1, 1)
    Φ2, dΦ2 = basis_values(q2, 0), basis_values(q2, 1)
    return kron(Φ2, dΦ1), kron(dΦ2, Φ1)
end

@doc raw"""
    compute_L!(L, sdist, landau; chunk = 256)

The symmetric matrix

```math
\mathbb{L}_{ij} = \frac{1}{2} \int_\Omega \! \int_\Omega
  \big( \nabla \varphi_i(v) - \nabla \varphi_i(v') \big) \cdot
  f_s(v) \, U(v,v') \, f_s(v') \cdot
  \big( \nabla \varphi_j(v) - \nabla \varphi_j(v') \big) \, dv \, dv'
```

of the Landau manuscript's `eq:discrete-landau-matrix`.

# Method

Expanding the two differences gives four terms. With ``U`` symmetric in its arguments and in
its indices, relabelling ``v \leftrightarrow v'`` shows the first and fourth are equal and the
second and third are equal, so the factor of one half cancels and

```math
\mathbb{L} = \sum_{c,d} D^c \operatorname{diag}(s \odot A^{cd}) (D^d)^{\mathsf T}
           - \sum_{c,d} \big( D^c \operatorname{diag}(s) \big) \, U^{cd} \,
             \big( D^d \operatorname{diag}(s) \big)^{\mathsf T} ,
\qquad A^{cd}_a = \sum_b s_b \, U^{cd}_{ab} ,
```

with ``s_a = w_a f_s(v_a)`` and ``D^c`` the sparse gradient tables above. This is an identity,
not an approximation: the same quadrature sum, rearranged.

Both terms are sparse-times-dense products. The first is a single pass over the grid; the
second is the genuinely coupled one, and its ``Q \times Q`` kernel is built and consumed in
row blocks of `chunk` so that nothing of that size is ever held.

# What this replaces

The earlier implementation filled ``M^2/2`` entries, and for **each** of them called a
`gauss_quad` that looped over every cell **quadruple** and then every node quadruple, applying
the support test *inside* the integrand so the iterations were executed and discarded. That is
``O(M^2 \, n^4 \, n_q^4)`` against the ``O(Q^2)`` here — at ``M = 100``, ``n = 10``,
``n_q = 3`` roughly ``4 \times 10^9`` integrand calls through closures, against ``4 \times
10^7`` flops in BLAS. It is also why the driver scripts ran with `n = 1`, which is what made
the coincident-node problem in [`kernel`](@ref) fire on every diagonal cell.
"""
function compute_L!(L, sdist::SplineDistribution{T, XD, 2}, landau::Landau;
        chunk::Int = 256) where {T, XD}
    q = sdist.quadrature
    X = quadrature_nodes(q)
    W = quadrature_weights(q)
    Q1, Q2 = length(X[1]), length(X[2])
    Q = Q1 * Q2
    M = length(sdist)
    fs = sdist.spline

    # The grid, enumerated with the first axis fastest -- the same order the Kronecker
    # products in `gradient_tabulations` produce for their columns.
    pts = Vector{SVector{2, T}}(undef, Q)
    s = Vector{T}(undef, Q)
    r = 0
    for a2 in 1:Q2, a1 in 1:Q1

        r += 1
        p = SVector(X[1][a1], X[2][a2])
        pts[r] = p
        s[r] = W[1][a1] * W[2][a2] * fs(p)
    end

    D = gradient_tabulations(sdist)

    # A^{cd}_a = Σ_b s_b U^{cd}(v_a, v_b). Symmetric in (c,d), so three components.
    A11 = zeros(T, Q)
    A12 = zeros(T, Q)
    A22 = zeros(T, Q)
    for a in 1:Q
        t11 = zero(T)
        t12 = zero(T)
        t22 = zero(T)
        for b in 1:Q
            U = kernel(pts[a], pts[b], landau)
            sb = s[b]
            t11 += sb * U[1, 1]
            t12 += sb * U[1, 2]
            t22 += sb * U[2, 2]
        end
        A11[a] = t11
        A12[a] = t12
        A22[a] = t22
    end

    Acd = ((A11, A12), (A12, A22))

    fill!(L, zero(T))

    # First term: a diagonal weighting of the grid, one sparse contraction per component pair.
    for c in 1:2, d in 1:2

        L .+= Matrix(D[c] * Diagonal(s .* Acd[c][d]) * D[d]')
    end

    # Second term: the coupled one. Xc = D^c diag(s), and the Q×Q kernel is built in row
    # blocks so that only `chunk × Q` of it exists at a time.
    Xs = (D[1] * Diagonal(s), D[2] * Diagonal(s))
    Ublock = Matrix{T}(undef, min(chunk, Q), Q)

    for lo in 1:chunk:Q
        hi = min(lo + chunk - 1, Q)
        nb = hi - lo + 1

        for c in 1:2, d in 1:2

            Ub = view(Ublock, 1:nb, :)
            for (ii, a) in enumerate(lo:hi)
                @inbounds for b in 1:Q
                    Ub[ii, b] = kernel(pts[a], pts[b], landau)[c, d]
                end
            end
            # (M × nb) * (nb × Q) * (Q × M)
            L .-= Matrix(view(Xs[c], :, lo:hi) * Ub * Xs[d]')
        end
    end

    # Symmetric by construction; enforce it exactly, as the manuscript's L_ij is symmetric and
    # the positive semi-definiteness of the bracket is stated for the symmetric form.
    for i in 1:M, j in 1:(i - 1)

        m = (L[i, j] + L[j, i]) / 2
        L[i, j] = m
        L[j, i] = m
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
    compute_K!(cache.K1, cache.K2, v, sdist, landau)

    # compute L_ij matrix
    compute_L!(cache.L, sdist, landau)

    # compute J vector
    compute_J!(cache.J, sdist, landau)

    # solve for vector field
    mul!(cache.LJ, cache.L, cache.J)

    # `\` on a rank-deficient least-squares problem is the minimum-norm solution — Julia
    # routes it through `qr(A, ColumnNorm())` and LAPACK's complete-orthogonal path, which is
    # xGELSY — so this is `pinv(K) * LJ` without forming the pseudo-inverse.
    #
    # The step from eq(147) to eq:particle_ode_final in the appendix uses K K⁺ = I, which needs
    # K to have full row rank. The earlier version computed `rank(K1)` and `rank(K2)` — two
    # SVDs of an M×N matrix, on *every* vector-field evaluation — printed a line when the
    # hypothesis failed, and then proceeded with the simplified formula anyway. Detecting the
    # failure of a hypothesis and continuing regardless is the thing to avoid; the rank check
    # is therefore gone from the hot path, and belongs in a diagnostic script instead.
    v̇[1, :] .= cache.K1 \ cache.LJ
    v̇[2, :] .= cache.K2 \ cache.LJ

    # The collision frequency was declared and stored but never read anywhere in this file, so
    # `Landau(dist, ent; ν = 0.1)` silently ran at ν = 1. Every other collision model here
    # multiplies by it.
    v̇ .*= landau.ν

    return nothing
end
