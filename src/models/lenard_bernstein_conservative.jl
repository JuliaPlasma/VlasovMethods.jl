struct CLBCache{T, PT <: ParticleDistribution, ST <: SplineDistribution{T}} <: Cache{T}
    pdist::PT
    sdist::ST

    function CLBCache{T}(pdist, sdist) where {T}
        new{T, typeof(pdist), typeof(sdist)}(pdist, sdist)
    end
end

function CLBCache(pdist::ParticleDistribution{T}, sdist::SplineDistribution{T}) where {T}
    CLBCache{T}(pdist, sdist)
end

function Cache(AT, c::CLBCache{DT, PT, ST}) where {DT, PT, ST}
    CLBCache{AT}(c.pdist, similar(AT, c.sdist))
end
function CacheType(AT, c::CLBCache{DT, PT, ST}) where {DT, PT, ST}
    CLBCache{AT, PT, similar_type(AT, c.sdist)}
end

struct ConservativeLenardBernstein{
    XD, VD, DT <: DistributionFunction{XD, VD}, ET <: Entropy, T, CT <: CacheDict} <:
       CollisionOperator
    dist::DT    # distribution function
    ent::ET     # entropy 
    ν::T        # collision frequency 

    cache::CT

    function ConservativeLenardBernstein(
            dist::DistributionFunction{XD, VD}, ent::Entropy; ν::T = 1.0) where {XD, VD, T}
        # `scripts/verify_conservation.jl` measures what this warning is about, and the
        # distinction it draws is worth stating exactly, because the obvious reading is wrong.
        #
        # The particle sums Σ w v̇ and Σ w v v̇ vanish at round-off on *any* basis: A₁ and A₂ are
        # solved from precisely those two constraints, which are algebraic conditions on the
        # particle sums and involve the basis only through f_s'/f_s. Momentum and energy of the
        # *particles* are therefore conserved whatever boundary condition is used.
        #
        # What the span requirement controls is the agreement between the two representations:
        # ∫ vᵏ f_s dv equals Σ_α w_α v_α^k only if vᵏ lies in the span, and the manuscript's
        # sentence at main.tex:226 is about that. Measured on a cubic basis over [-8,8]: the
        # clamped basis reproduces 1, v and v² to 1e-14, a periodic one gets v wrong by 1.2 and
        # v² by 0.12, and a Dirichlet-recombined one gets even the constant wrong by 4e-2. So on
        # a bad basis the diagnostics computed from f_s disagree with the particles', and the
        # H-theorem — which is a statement about f_s — rests on a distribution whose moments are
        # not the ones being conserved.
        m = polynomial_reproduction(ent.dist.basis)
        if m < 2
            @warn """
            ConservativeLenardBernstein built on a basis that reproduces polynomials only up \
            to degree $(m). The particle momentum and energy are still conserved — A₁ and A₂ \
            are solved from those two constraints directly — but the moments of the spline \
            representation f_s no longer agree with the particles': ∫ vᵏ f_s dv = Σ w v^k \
            requires vᵏ in the span. Any diagnostic read off f_s, the entropy included, is \
            therefore measuring a distribution with the wrong moments. Free() gives degree p, \
            Periodic() gives 0 because v is not periodic, Dirichlet() gives -1 because every \
            basis function vanishes at the ends.""" maxlog=1
        end
        cache = CacheDict(CLBCache(dist, ent.dist))
        new{XD, VD, typeof(dist), typeof(ent), T, typeof(cache)}(dist, ent, ν, cache)
    end
end

function Cache(AT, clb::ConservativeLenardBernstein)
    CLBCache{AT}(clb.dist, similar(AT, clb.ent.dist))
end
function CacheType(AT, clb::ConservativeLenardBernstein)
    CLBCache{AT, typeof(clb.dist), similar_type(AT, clb.ent.dist)}
end

@doc raw"""
    compute_coefficients(distribution, particle_dist, vp)

The coefficients ``A_1``, ``A_2`` of the conservative Lenard-Bernstein drift, from the
requirement that the operator conserve momentum and energy:

```math
\begin{pmatrix} n_h & n_h u_h \\ n_h u_h & n_h \varepsilon_h \end{pmatrix}
\begin{pmatrix} A_1 \\ A_2 \end{pmatrix}
= - \sum_\alpha w_\alpha \begin{pmatrix} 1 \\ v_\alpha \end{pmatrix}
  \frac{f_s'(v_\alpha)}{f_s(v_\alpha)} ,
\qquad
\begin{aligned}
n_h &= \textstyle\sum_\alpha w_\alpha , \\
n_h u_h &= \textstyle\sum_\alpha w_\alpha v_\alpha , \\
n_h \varepsilon_h &= \textstyle\sum_\alpha w_\alpha v_\alpha^2 ,
\end{aligned}
```

solved by Cramer's rule. This is `eq:coefficient_lin_system` of the Lenard-Bernstein
manuscript.

!!! note "The particle weights belong in every sum"
    They were absent from all five of them. Both the matrix and the right-hand side are linear
    in ``w``, so for **uniform** weights the factor cancels and ``A_1, A_2`` come out
    unchanged — which is why this never showed up: every initialiser in `examples/` sets
    `w = 1/npart`. What the earlier version actually annihilated was
    ``\sum_\alpha \dot{v}_\alpha`` rather than ``\sum_\alpha w_\alpha \dot{v}_\alpha``, so with
    the non-uniform weights the importance-sampling initialiser produces, momentum and energy
    conservation were lost.
"""
function compute_coefficients(
        distribution::SplineDistribution, particle_dist::ParticleDistribution,
        vp::AbstractArray{VT}) where {VT}
    w = particle_dist.particles.w
    fs = distribution.spline
    dfs = derivative(fs)

    n = zero(VT)
    nu = zero(VT)
    neps = zero(VT)
    B1 = zero(VT)
    B2 = zero(VT)

    for α in eachindex(vp)
        wα = w[1, α]
        vα = vp[α]
        fα = fs(vα)

        # f_s is an L² projection of a particle distribution, so it undershoots into negative
        # values where the sampling is thin and is exactly zero in cells with no particles.
        # f'/f is then either a division by zero or a term of the wrong sign, and the entropy
        # of the metriplectic form is a log of a negative number. Both manuscripts know this:
        # the Landau appendix opens with a commented-out section titled "Positivity-preserving
        # projections" stating that f_s "has to be non-negative everywhere", and the
        # Lenard-Bernstein H-theorem is stated "provided that f is non-negative". Neither
        # solves it. Failing here names the open problem instead of continuing with a value
        # whose sign is wrong.
        fα > 0 || throw(ErrorException(
            "the projected distribution is non-positive, f_s(v) = $(fα) at v = $(vα), so " *
            "f_s'/f_s is not meaningful there. This is the positivity problem both " *
            "manuscripts name and neither solves: an L² projection of a particle " *
            "distribution undershoots where the sampling is thin. Use more particles, a " *
            "coarser spline basis, or a positivity-preserving projection."))

        n += wα
        nu += wα * vα
        neps += wα * vα^2

        g = dfs(vα) / fα
        B1 -= wα * g
        B2 -= wα * vα * g
    end

    det = n * neps - nu^2
    A1 = (neps * B1 - nu * B2) / det
    A2 = -(nu * B1 - n * B2) / det

    if isnan(A1) || isnan(A2)
        # The velocity spread has collapsed: n εₕ = (n uₕ)² makes the 2×2 system singular.
        # The earlier message blamed the timestep, which is not the cause of any of the ways
        # this can happen.
        throw(ErrorException(
            "the coefficient system is singular: n = $(n), nu = $(nu), neps = $(neps), " *
            "determinant $(det). The particle velocities have no spread in energy about " *
            "their mean, so momentum and energy conservation are not independent conditions."))
    end

    return A1, A2
end

# function compute_coefficients(distribution::SplineDistribution{1,2}, particle_dist::ParticleDistribution, vp::AbstractArray{VT}) where {VT}
#     n = zeros(T, 2)
#     nu = similar(n)
#     neps = similar(n)
#     n, nu, neps = compute_f_densities(distribution, vp)
#     B1, B2 = compute_df_densities(distribution, vp)
#     B1 *= -1 
#     B2 *= -1 
#     A1 = (neps * B1 - nu * B2) / (n * neps - (nu)^2)
#     A2 = -  (nu * B1 - n * B2) / (n * neps - (nu)^2)

#     return A1, A2
# end

@doc raw"""
The conservative Lenard-Bernstein right-hand side,

```math
\dot{v}_\alpha = - \nu \left( \frac{f_s'(v_\alpha)}{f_s(v_\alpha)} + A_1 + A_2 v_\alpha \right) ,
```

which is `eq:velocity_ode` of the Lenard-Bernstein manuscript **with the opposite sign**.

!!! note "The sign here is the correct one and the manuscript's is not"
    The manuscript writes ``\partial_t f = \partial_v (a f)`` with
    ``a = \nu (f'/f + A_1 + A_2 v)`` and then ``\dot{v}_\alpha = a(v_\alpha)``. Continuity is
    ``\partial_t f + \partial_v (U f) = 0``, so ``\partial_t f = +\partial_v(a f)`` forces
    ``U = -a`` and hence ``\dot{v} = -a``. The check: for a Maxwellian of width
    ``\sigma > 1`` — too broad — ``a = \nu v (1 - \sigma^{-2}) > 0`` for ``v > 0``, so
    ``\dot v = +a`` would push particles outward and *widen* an already too-wide distribution
    instead of relaxing it. Conservation is unaffected either way, the constraints being
    homogeneous, so no published result changes; but the ODE as typeset relaxes backwards in
    time.
"""
function CLB_rhs!(v̇, v::AbstractVector{ST}, params, t) where {ST}
    # `.sdist`: the cache is a CLBCache, not a SplineDistribution, so the projection below
    # would have been a MethodError.
    dist = params.model.cache[ST].sdist

    fs = projection(v, params.idist, dist)
    dfdv = derivative(fs)

    A = compute_coefficients(dist, params.idist, v)
    v̇ .= -params.ν .* (dfdv.(v) ./ fs.(v) .+ (A[1] .+ A[2] .* v))
end

function CLB_rhs_GI!(v, t, q::AbstractArray{ST}, params) where {ST}
    dist = params.model.cache[ST].sdist

    fs = projection(q, params.idist, dist)
    dfdv = derivative(fs)

    A = compute_coefficients(dist, params.idist, q)

    v .= -params.ν .* (dfdv.(q) ./ fs.(q) .+ (A[1] .+ A[2] .* q))
end

# used for plotting
function CLB_rhs(v::AbstractVector{ST}, params, fs::Spline) where {ST}
    dist = params.model.cache[ST].sdist

    dfdv = derivative(fs)

    A = compute_coefficients(dist, params.idist, v)

    return -params.ν .* (dfdv.(v) ./ fs.(v) .+ (A[1] .+ A[2] .* v))
end

function DiffEqIntegrator(model::ConservativeLenardBernstein{1, 1}, tspan::Tuple, tstep::Real)
    # parameters for computing vector field
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # u0 = copy(model.dist.particles.v[1,:])
    # construct DifferentialEquations ODEProblem
    equ = DifferentialEquations.ODEProblem(
        CLB_rhs!,
        copy(model.dist.particles.v[1, :]),
        tspan,
        params
    )

    # choose integrator
    # int = DifferentialEquations.TRBDF2()
    int = DifferentialEquations.Trapezoid()

    DiffEqIntegrator(model, equ, int, tstep)
end

function GeometricIntegrator(model::ConservativeLenardBernstein, tspan::Tuple, tstep::Real)
    # collect parameters
    # params = (ϕ = model.potential, model = model)
    params = (ν = model.ν, idist = model.dist, fdist = model.ent.dist, model = model)
    # create geometric problem
    equ = GeometricEquations.ODEProblem(
        CLB_rhs_GI!,
        tspan, tstep, copy(model.dist.particles.v[1, :]);
        parameters = params)

    # create integrator
    int = Integrators.GeometricIntegrator(equ, Integrators.ImplicitMidpoint())
    # int = Integrators.Integrator(equ, Integrators.ImplicitMidpoint())
    # int = Integrators.Integrator(equ, Integrators.RK438())
    # int = Integrators.Integrator(equ, Integrators.CrankNicolson())

    # put together splitting method
    GeometricIntegrator(model, equ, int)
    # return int
end
