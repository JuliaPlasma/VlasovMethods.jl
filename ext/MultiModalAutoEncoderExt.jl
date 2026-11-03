module MultiModalAutoEncoderExt

using GeometricMachineLearning: BilinearStiefelLayer, LowRankStiefelLayer
using LinearAlgebra: Diagonal, Symmetric, cholesky, cond, dot, opnorm, svd, svdvals
using SimpleSplines: AbstractBSplineBasis, BSplineBasis, PeriodicBSplineBasis, Spline,
                     SplineQuadrature,
                     basis, basis_index, coefficients, degree, derivative_matrix,
                     domainlength,
                     evaluate_all!, knotvector, mass_matrix, nbasis, stiffness_matrix,
                     weighted_matrix

import VlasovMethods: build_spline_parity_factors, build_grid_parity_factors,
                      tie_continuum_frontends!, dual_parity_modes, klimontovich_deposition,
                      klimontovich_bridge_loss, kruskal_oberman_free_energy

# The parity split of `BilinearStiefelLayer` as an orthogonal matrix: row `k ≤ h` is the even
# combination of index `k` with its mirror `n + 1 - k`, row `h + k` the odd one, both scaled by
# 1/√2. `S * c` gives the even half on top of the odd half.
function _parity_matrix(::Type{T}, n::Integer) where {T}
    h = n ÷ 2
    S = zeros(T, n, n)
    s = inv(sqrt(T(2)))
    for k in 1:h
        S[k, k] = s
        S[k, n + 1 - k] = s
        S[h + k, k] = s
        S[h + k, n + 1 - k] = -s
    end
    return S
end

_dense_mass(b::AbstractBSplineBasis) = Matrix(mass_matrix(SplineQuadrature(b)))

# The velocity axis must be mirror-symmetric about v = 0, with an even number of points or
# basis functions, for the parity split to pair index `k` with the reflection of its point.
function _check_symmetric(points::AbstractVector{T}, what::AbstractString) where {T}
    scale = max(abs(first(points)), abs(last(points)))
    tol = 10 * eps(T) * scale
    abs(first(points) + last(points)) ≤ tol ||
        throw(ArgumentError("the velocity domain [$(first(points)), $(last(points))] is not " *
                            "symmetric: v_min must equal -v_max"))
    all(k -> abs(points[k] + points[end + 1 - k]) ≤ tol, eachindex(points)) ||
        throw(ArgumentError("the velocity $(what) is not symmetric under v ↦ -v"))
    return nothing
end

function build_spline_parity_factors(basis_x::AbstractBSplineBasis{T},
        basis_v::BSplineBasis{T}, mx::Integer, mv::Integer) where {T}
    _check_symmetric(knotvector(basis_v), "knot vector")
    nv = nbasis(basis_v)
    iseven(nv) ||
        throw(ArgumentError("the velocity basis has $(nv) functions, an odd number: the " *
                            "parity split needs an even count, so pick N cells with N + p even"))
    Rx = cholesky(Symmetric(_dense_mass(basis_x))).U
    S = _parity_matrix(T, nv)
    P = S * _dense_mass(basis_v) * S'
    h = nv ÷ 2
    Rp = cholesky(Symmetric(P[1:h, 1:h])).U
    Rm = cholesky(Symmetric(P[(h + 1):end, (h + 1):end])).U
    Rpm = zeros(T, nv, nv)
    Rpm[1:h, 1:h] .= Rp
    Rpm[(h + 1):end, (h + 1):end] .= Rm
    Rv = S' * Rpm * S
    return BilinearStiefelLayer(nbasis(basis_x), nv, mx, mv;
        metric_factors = (Matrix(Rx), Rv))
end

function build_grid_parity_factors(grid_x::AbstractRange{T}, grid_v::AbstractRange{T},
        mx::Integer, mv::Integer) where {T}
    _check_symmetric(grid_v, "grid")
    nv = length(grid_v)
    iseven(nv) ||
        throw(ArgumentError("the velocity grid has $(nv) points, an odd number: the parity " *
                            "split needs an even count"))
    Wx = Diagonal(fill(sqrt(T(step(grid_x))), length(grid_x)))
    Wv = Diagonal(fill(sqrt(T(step(grid_v))), nv))
    return BilinearStiefelLayer(length(grid_x), nv, mx, mv; metric_factors = (Wx, Wv))
end

# The orthonormal matrix nearest to `A`, its polar factor.
function _polar(A::AbstractMatrix)
    F = svd(A)
    return F.U * F.Vt
end

# A `StiefelManifold` keeps its matrix in the field `A` and defines no `setindex!`.
function _copy_basis!(dest::AbstractMatrix, src::AbstractMatrix)
    return hasproperty(dest, :A) ? copyto!(dest.A, src) : copyto!(dest, src)
end

function tie_continuum_frontends!(spline::Tuple{BilinearStiefelLayer, NamedTuple},
        grid::Tuple{BilinearStiefelLayer, NamedTuple},
        dlra::Tuple{LowRankStiefelLayer, NamedTuple},
        E_x::AbstractMatrix, E_v::AbstractMatrix)
    spline_layer, ps_spline = spline
    grid_layer, ps_grid = grid
    dlra_layer, ps_dlra = dlra
    Rx, Rv = spline_layer.metric_factors
    Wx, Wv = grid_layer.metric_factors
    dlra_layer.metric_factors == grid_layer.metric_factors ||
        throw(ArgumentError("the DLRA layer must carry the metric factors of the grid layer"))
    size(E_x) == (size(Wx, 1), size(Rx, 1)) ||
        throw(DimensionMismatch("E_x must have size $((size(Wx, 1), size(Rx, 1))); got $(size(E_x))"))
    size(E_v) == (size(Wv, 1), size(Rv, 1)) ||
        throw(DimensionMismatch("E_v must have size $((size(Wv, 1), size(Rv, 1))); got $(size(E_v))"))
    T = eltype(Rv)
    nv, ng = size(Rv, 1), size(Wv, 1)
    Sv = _parity_matrix(T, nv)[1:(nv ÷ 2), :]
    Sg = _parity_matrix(T, ng)[1:(ng ÷ 2), :]
    # A spline mode in L²-orthonormal coordinates `a` is the function with coefficients
    # `R \ a`; its grid coordinates are its nodal values times √w. The velocity modes are the
    # even ones, which the layer's parity split reaches through the top half of `S`.
    Φx = _polar(Wx * (E_x * (Rx \ Matrix(ps_spline.Φx))))
    Φv = _polar(Sg * (Wv * (E_v * (Rv \ (Sv' * Matrix(ps_spline.Φv))))))
    for (name, ps) in (("ps_grid", ps_grid), ("ps_dlra", ps_dlra)),
        (f, Φ) in ((:Φx, Φx), (:Φv, Φv))

        size(ps[f]) == size(Φ) ||
            throw(DimensionMismatch("$(name).$(f) must have size $(size(Φ)); got $(size(ps[f]))"))
    end
    for ps in (ps_grid, ps_dlra)
        _copy_basis!(ps.Φx, Φx)
        _copy_basis!(ps.Φv, Φv)
    end
    return ps_grid, ps_dlra
end

# One basis for a vector of splines, and their coefficients as the columns of a matrix.
function _modes(splines::AbstractVector{<:Spline}, what::AbstractString)
    b = basis(first(splines))
    all(s -> basis(s) == b, splines) ||
        throw(ArgumentError("the $(what) modes must share one basis"))
    return b, reduce(hcat, coefficients.(splines))
end

function dual_parity_modes(x_modes::AbstractVector{<:Spline},
        even_v::AbstractVector{<:Spline}, odd_v::AbstractVector{<:Spline}, df0;
        nq::Integer = 24)
    length(even_v) == length(odd_v) ||
        throw(ArgumentError("the even and odd velocity spaces have dimensions " *
                            "$(length(even_v)) and $(length(odd_v)): they must be equal"))
    bx, Φx = _modes(x_modes, "x")
    bv, Ψp = _modes(even_v, "even velocity")
    bm, Ψm = _modes(odd_v, "odd velocity")
    bm == bv ||
        throw(ArgumentError("the odd velocity modes must share the basis of the even ones"))
    T = eltype(Φx)
    qx = SplineQuadrature(bx)
    Mx = Φx' * mass_matrix(qx) * Φx
    Kx = Φx' * stiffness_matrix(qx) * Φx
    L = domainlength(bx)
    # ∫ φ'² / ∫ φ² is k² for a mode of wavenumber k, at least (2π/L)² for one that is not
    # constant; the round-off of a projected constant leaves it near (eps/h)².
    all(i -> Kx[i, i] * L^2 > sqrt(eps(T)) * Mx[i, i], axes(Kx, 1)) ||
        throw(ArgumentError("an x mode is constant: ∂x is not invertible on the x space"))
    # (D_x)_{ik} = ∫ φᵢ' φₖ, and `derivative_matrix` is ∫ Λₐ Λ_b'. Its singular values are the
    # wavenumbers of the modes, so a singular `Dx` is one whose smallest singular value is at the
    # round-off of ‖∂x‖ on the modes. The condition number does not serve: a `Dx` whose entries
    # are all round-off can have a small one.
    Dx = Φx' * derivative_matrix(qx)' * Φx
    minimum(svdvals(Dx)) > 100 * eps(T) * sqrt(opnorm(Kx) * opnorm(Mx)) ||
        throw(ArgumentError("∂x is singular on the x space: each cos mx needs its sin mx"))
    qv = SplineQuadrature(bv; nq = nq)
    Kv = Ψp' * weighted_matrix(qv, v -> T(df0(v)), 0, 0) * Ψm
    cond(Kv) ≤ inv(100 * eps(T)) ||
        throw(ArgumentError("K_v = ∫ f₀' ψ⁺ ψ⁻ is numerically singular: " *
                            "cond(K_v) > 1/(100 eps)"))
    B = -kron(Kv, Dx)
    return kron(Ψm, Φx) / B
end

# The deposition: the dense Cholesky factors of the two mass matrices, the weights, and the
# two bases the particles are evaluated in.
struct KlimontovichDepositionLayer{T, BX <: AbstractBSplineBasis{T},
    BV <: AbstractBSplineBasis{T}, FT, WT <: AbstractVector{T}}
    basis_x::BX
    basis_v::BV
    Mx::FT
    Mv::FT
    weights::WT
end

function klimontovich_deposition(
        basis_x::PeriodicBSplineBasis{T}, basis_v::BSplineBasis{T},
        weights::AbstractVector{T}) where {T}
    Mx = cholesky(Symmetric(_dense_mass(basis_x)))
    Mv = cholesky(Symmetric(_dense_mass(basis_v)))
    return KlimontovichDepositionLayer(basis_x, basis_v, Mx, Mv, weights)
end

function (l::KlimontovichDepositionLayer{T})(x::AbstractVector{T},
        v::AbstractVector{T}) where {T}
    length(x) == length(v) == length(l.weights) ||
        throw(DimensionMismatch("x, v and the weights must have one entry per particle; " *
                                "got $(length(x)), $(length(v)) and $(length(l.weights))"))
    bx, bv = l.basis_x, l.basis_v
    c = zeros(T, nbasis(bx), nbasis(bv))
    ex = zeros(T, degree(bx) + 1)
    ev = zeros(T, degree(bv) + 1)
    for i in eachindex(x, v, l.weights)
        jx = evaluate_all!(ex, bx, x[i])
        jv = evaluate_all!(ev, bv, v[i])
        # The block of a periodic or a clamped basis has its p + 1 indices in range, a point
        # outside the clamped domain included: its block is that of the end cell, with zeros.
        for t in eachindex(ev)
            b = basis_index(bv, jv + t - 1)
            wv = l.weights[i] * ev[t]
            for s in eachindex(ex)
                c[basis_index(bx, jx + s - 1), b] += ex[s] * wv
            end
        end
    end
    return (l.Mx \ c) / l.Mv
end

# The bridge loss: the deposition, the spline frontend it feeds and the ε of the denominator.
struct KlimontovichBridgeLoss{
    DT <: KlimontovichDepositionLayer, LT <: BilinearStiefelLayer, T}
    deposition::DT
    spline_layer::LT
    ε::T
end

function klimontovich_bridge_loss(deposition::KlimontovichDepositionLayer{T},
        spline_layer::BilinearStiefelLayer; ε::Real = eps(T)) where {T}
    return KlimontovichBridgeLoss(deposition, spline_layer, T(ε))
end

function (L::KlimontovichBridgeLoss)(z::NamedTuple{(:q, :p)}, x::AbstractVector,
        v::AbstractVector, ps::NamedTuple)
    zs = L.spline_layer(L.deposition(x, v), ps)
    num = sum(abs2, z.q - zs.q) + sum(abs2, z.p - zs.p)
    return num / (sum(abs2, zs.q) + sum(abs2, zs.p) + L.ε)
end

# The Kruskal–Oberman free energy on Fourier modes in x and nodes in v: the wavenumber of each
# row of a state, the nodes and their weight, and the weight w₀ = -v / f₀' at the nodes.
struct KruskalObermanFreeEnergy{T, KT <: AbstractVector{T}, VT <: AbstractVector{T}}
    k::KT
    v::VT
    Δv::T
    df0::VT
    w0::VT
end

function kruskal_oberman_free_energy(df0, grid_v::AbstractRange{T},
        wavenumbers::AbstractVector{<:Integer}) where {T}
    any(iszero, wavenumbers) &&
        throw(ArgumentError("the x space holds the constant mode k = 0, which has no " *
                            "conjugate: pass the wavenumbers k ≠ 0 only"))
    allunique(abs.(wavenumbers)) ||
        throw(ArgumentError("the wavenumbers repeat a mode: k and -k give the same modes " *
                            "cos kx, sin kx, so each |k| may occur once"))
    _check_symmetric(grid_v, "grid")
    iseven(length(grid_v)) ||
        throw(ArgumentError("the velocity grid has $(length(grid_v)) points, an odd " *
                            "number: the parity split needs an even count"))
    v = collect(grid_v)
    f = T.(df0.(v))
    all(q -> v[q] ≤ 0 || f[q] < 0, eachindex(v)) ||
        throw(ArgumentError("f₀' ≥ 0 at a node v > 0, where the weight -v/f₀' of the free " *
                            "energy must be positive; for a non-monotone f₀ use dual_parity_modes"))
    k = T.(repeat(wavenumbers; inner = 2))
    return KruskalObermanFreeEnergy(k, v, T(step(grid_v)), f, -v ./ f)
end

function (H::KruskalObermanFreeEnergy)(C::AbstractMatrix)
    size(C) == (length(H.k), length(H.v)) ||
        throw(DimensionMismatch("the state must have size $((length(H.k), length(H.v))); " *
                                "got $(size(C))"))
    ρ = vec(sum(C; dims = 2)) .* H.Δv
    Φ = ρ ./ H.k .^ 2
    kinetic = H.Δv * sum(q -> H.w0[q] * sum(abs2, view(C, :, q)), axes(C, 2))
    value = (kinetic + dot(ρ, Φ)) / 2
    # δH/δf = w₀ δf + Φ, the gradient with respect to `C` divided by the node weight
    U = C .* H.w0' .+ Φ
    gradient = H.Δv .* U
    # The skew operator [0 𝒜; 𝒜 0] on the parity parts of δH/δf: 𝒜 = f₀' ∂ₓ maps the odd part
    # to the rate of δf⁺ and the even part to the rate of δf⁻, since f₀' is odd. 𝒜 is linear, so
    # the rate of δf = δf⁺ + δf⁻ is 𝒜 applied to the whole of δH/δf.
    rhs = _advect(H, U)
    return (value = value, gradient = gradient, rhs = rhs)
end

# 𝒜 U = f₀'(v) ∂ₓ U on the rows cos kx, sin kx: ∂ₓ takes the pair (a, b) to (k b, -k a).
function _advect(H::KruskalObermanFreeEnergy, U::AbstractMatrix)
    dU = similar(U)
    @views for i in 1:2:size(U, 1)
        k = H.k[i]
        dU[i, :] .= k .* U[i + 1, :] .* H.df0
        dU[i + 1, :] .= -k .* U[i, :] .* H.df0
    end
    return dU
end

end
