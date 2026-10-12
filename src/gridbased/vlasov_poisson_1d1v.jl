@doc raw"""
    _green_matrix(dist::GridDistribution, potential::Potential)

The ``n_x \times n_x`` Green's matrix ``G`` of `potential`'s backend on the grid of `dist`, which
gives the field of a density in one product:

```math
\rho_i = h_v \sum_j f_{ij} , \qquad \phi = G \rho .
```

It solves the backend for the ``n_x`` unit vectors and symmetrises the result. The symmetry
makes ``\partial H / \partial f = h_x h_v h`` exact, so the discrete Hamiltonian is quadratic in
``f`` and an implicit midpoint step conserves it. A spline basis is refused: its density at the
nodes is not a Galerkin load vector. `PoissonSolverMatrixFree` bounds its residual, not its
error, so in `Float32` on a fine grid `G` is less accurate than ``eps(Float32)``. The errors are
listed under [`VlasovPoisson`](@ref).
"""
function _green_matrix(dist::GridDistribution, potential::Potential)
    b = basis(potential)
    solver = potential.solver
    (solver isa PoissonSolverFFT || solver isa PoissonSolverMatrixFree) ||
        throw(ArgumentError(
            "the grid 1D1V Vlasov–Poisson scheme solves the field through PoissonSolvers'" *
            " PoissonSolverFFT or PoissonSolverMatrixFree; got a $(typeof(solver))"))

    nx, nv = size(dist.values)
    length(b) == nx || throw(DimensionMismatch(
        "the potential has $(length(b)) degrees of freedom on a distribution of " *
        "$nx x-nodes"))

    b.xgrid[1:nx] == dist.x || throw(ArgumentError(
        "the potential's x-grid $(collect(b.xgrid[1:nx])) is not the distribution's " *
        "$(collect(dist.x)); build the distribution from the same domain, in the potential's " *
        "element type, with `GridDistribution(values, xdomain, vdomain)`"))

    # Fewer than three nodes per direction is Arakawa's own check, which the bracket the run
    # uses performs too; this one surfaces it at construction. It is discarded.
    Arakawa(nx, nv, step(dist.x), step(dist.v))

    # the unit vectors are solved in the potential's element type and `G` is stored in the
    # distribution's; the two differ only where the exact grid comparison above holds in both
    T, S = eltype(dist.values), eltype(b)
    G = stack(T.(PoissonSolvers.solve(solver, S.(1:nx .== c))) for c in 1:nx)

    return (G + G') / 2
end

@doc raw"""
    _grid_1d1v_rhs!(ḟ, t, f, params)

Arakawa's discretisation of the 1D1V Vlasov–Poisson right-hand side, as the vector field of a
`GeometricEquations.ODEProblem`:

```math
\dot{f} = - [f, h] , \qquad h_{ij} = \frac{v_j^2}{2} + \phi_i , \qquad \phi = G \rho ,
\qquad \rho_i = h_v \sum_j f_{ij} ,
```

for a charge of ``+1``, with `G` from [`_green_matrix`](@ref). `params` holds the bracket
`arakawa`, the velocity nodes `v` and `G`; the state `f` is `vec(dist.values)`. The bracket is
doubly periodic, so a distribution that reaches the ends of the velocity range wraps mass from
``v_b`` to ``v_a`` (see [`VlasovPoisson`](@ref)).
"""
function _grid_1d1v_rhs!(ḟ, t, f, params)
    (; arakawa, v, G) = params
    ϕ = G * (arakawa.hv .* vec(sum(reshape(f, arakawa.nx, arakawa.nv); dims = 2)))
    h = vec(ϕ .+ v' .^ 2 ./ 2)
    ḟ .= poisson_apply(arakawa, f, h) ./ (arakawa.hx * arakawa.hv)

    return ḟ
end

@doc raw"""
    GeometricIntegrator(model::VlasovPoisson{1, 1, <:GridDistribution}, tspan, tstep,
                        method = Integrators.ImplicitMidpoint())

The `GeometricIntegrator` that steps the grid 1D1V Vlasov–Poisson model `model` over `tspan`
with step `tstep`, in the state `vec(model.distribution.values)`.

`method` is any `GeometricIntegrators` method for an `ODEProblem`. `ImplicitMidpoint()` keeps
``H``; `RK4()` is the one `scripts/landau_damping_1d1v.jl` uses, because the implicit midpoint's
dense Newton Jacobian is too large on that script's ``n_v``. The scheme runs in `Float32` and
`Float64`; FFTW has no `Float16` plan. The vector field is [`_grid_1d1v_rhs!`](@ref).
"""
function GeometricIntegrator(model::VlasovPoisson{1, 1, <:GridDistribution}, tspan::Tuple,
        tstep::Real, method = Integrators.ImplicitMidpoint())
    dist = model.distribution
    params = (arakawa = Arakawa(size(dist.values)..., step(dist.x), step(dist.v)),
        v = dist.v, G = model.work)

    equ = GeometricEquations.ODEProblem(
        _grid_1d1v_rhs!, tspan, tstep, copy(vec(dist.values));
        parameters = params)

    int = Integrators.GeometricIntegrator(equ, method)

    return GeometricIntegrator(model, equ, int)
end

# The grid model's state is the whole node array rather than a column of positions, so its
# write-back reshapes instead of touching `dist.particles.v`. The particle models' method is in
# `src/methods/geometric_integrator.jl`.
function _write_back!(model::VlasovPoisson{1, 1, <:GridDistribution}, q)
    copyto!(model.distribution.values, q)

    return model.distribution
end
