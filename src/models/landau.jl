struct Landau{XD, VD, DT <: DistributionFunction{XD,VD}, ET <: Entropy, T} <: VlasovModel
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 
    
    function Landau(dist::DistributionFunction{XD,VD}, ent::Entropy; ν::T=1.) where {XD, VD, T}
        new{XD, VD, typeof(dist), typeof(ent), T}(dist, ent, ν)
    end
end

function stuff(S::SplineDistribution{1,2}, i, j, v::AbstractVector{T}) where T
    B = S.basis
    # @show v
    a = B[i,T](v[1])
    b = B[j,T](v[2])

    Sv = S.spline(v)
    if a == 0. || b == 0.
        return zero(T)
    # elseif Sv < 0.0 
    #     println("WARNING: NEGATIVE f values")
    #     return zero(T)
    else 
        # vgrid = -6:0.1:6
        # f_plot(x,y) = S.spline(x,y)
        # plot(vgrid, vgrid, f_plot, st=:surface, xlabel="x", ylabel="y")
        # savefig("spline-test")
        return a * b * (one(T) .+ log(Sv.^2) / 2.)
    end
end

function compute_matrices(distribution::SplineDistribution{1,2})
    # B = distribution.basis
    n = length(distribution.basis)
    k = BSplineKit.order(distribution.basis) 
    d_start = BSplineKit.knots(distribution.basis)[1] + 2.
    d_end = BSplineKit.knots(distribution.basis)[end] - 2.

    # int1 = zeros(eltype(distribution),n)
    int2 = zeros(eltype(distribution),(n,n))

    invm = inv(distribution.mass_fact)

    for i in eachindex(distribution.basis)
        # int1[i], _ = quadgk(v -> distribution.basis[i](v), d_start, d_end, rtol=1e-8)
        for j in eachindex(distribution.basis)
            if j ∈ i-k+1:i+k-1
                # println("i=", i, ", j =", j)
                # f(v::AbstractVector) = B[i](v[1]) * B[j](v[2]) * log(distribution.spline(v))
                integrand(v::AbstractVector) = stuff(distribution, i, j, v)
                int2[i,j], _ = hcubature(integrand, [d_start, d_start], [d_end, d_end], rtol=1e-6)
            end
        end
    end

    L_ij = invm * int2 * invm
    # K_ij = invm * int1 * int1' * invm
 
    return L_ij
end

function compute_U(v_α::Vector{T}, v_β::Vector{T}) where T
    n = length(v_α)
    U = zeros(T,(n,n))

    if v_α != v_β
        for i in CartesianIndices(U)
            if i[1] == i[2]
                U[i] += 1/norm(v_α - v_β)
            end
            U[i] -= (v_α[i[1]] - v_β[i[1]])*(v_α[i[2]] - v_β[i[2]])./norm(v_α - v_β)^3
        end
    end
    return U
end

function Landau_rhs!(v̇, v::AbstractArray{ST}, params, t) where {ST}
    @show t
    # @show v̇
    # @show v
    dist = params.model.ent.cache[ST]
    S = projection(v, params.idist, dist)
    B = basis(S)
    k = BSplineKit.order(B)
    B_t = BSplineKit.knots(B)
    L = compute_matrices(dist)
    # L, K = compute_matrices(dist)

    for γ in axes(v̇,2)
        ilast1, _ = find_knot_interval(B_t, v[1, γ])
        ilast2, _ = find_knot_interval(B_t, v[2, γ])
        v̇[:,γ] .= zeros(ST, 2)
        for α in axes(v,2)
            A = zeros(ST, 2)
            ilast3, _ = find_knot_interval(B_t, v[1, α])
            ilast4, _ = find_knot_interval(B_t, v[2, α])

            for i in ilast1:-1:ilast1-k+1
                # A = sum((L_ij[i,j] + K_ij[i,j]) * (eval_bfd(B, i, j, v[:,γ]) -  eval_bfd(B, i, j, v[:,α])) for j in eachindex(B), i in eachindex(B))
                for j in ilast2:-1:ilast2-k+1
                    if i > 0 && j > 0 && i <= size(dist)[1] && j <= size(dist)[2]
                        A .+= (L[i,j]) * eval_bfd(B, i, j, v[:, γ])
                        # A .+= (L[i,j] + K[i, j]) * eval_bfd(B, i, j, v[:, γ])
                    end
                end
            end

            for i in ilast3:-1:ilast1-k+1
                # A = sum((L_ij[i,j] + K_ij[i,j]) * (eval_bfd(B, i, j, v[:,γ]) -  eval_bfd(B, i, j, v[:,α])) for j in eachindex(B), i in eachindex(B))
                for j in ilast4:-1:ilast2-k+1
                    if i > 0 && j > 0 && i <= size(dist)[1] && j <= size(dist)[2]
                        A .-= (L[i,j]) * eval_bfd(B, i, j, v[:, α])
                        # A .-= (L[i,j] + K[i, j]) * eval_bfd(B, i, j, v[:, α])
                    end
                end
            end
            # A = sum((L_ij[i,j] + K_ij[i,j]) * (eval_bfd(S, i, j, v[:,γ]) - eval_bfd(S, i, j, v[:,α])) for i in eachindex(B), j in eachindex(B))
            v̇[:,γ] -= params.idist.particles.w[1,α]*compute_U(v[:,γ], v[:,α])*A
        end
    end
end


function DiffEqIntegrator(model::Landau{1,2}, tspan::Tuple, tstep::Real)
    # parameters for computing vector field
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # u0 = copy(model.dist.particles.v[1,:])
    # construct DifferentialEquations ODEProblem
    equ = DifferentialEquations.ODEProblem(
        Landau_rhs!,
        copy(model.dist.particles.v),
        tspan,
        params
    )

    # choose integrator
    # int = DifferentialEquations.TRBDF2()
    # int = DifferentialEquations.Trapezoid()
    int = DifferentialEquations.Tsit5()
    # int = DifferentialEquations.ImplicitMidpoint(autodiff=false)

    DiffEqIntegrator(model, equ, int, tstep)
 end
