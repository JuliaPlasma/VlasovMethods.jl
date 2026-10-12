@doc raw"""
    _green_matrix(dist::GridDistribution, potential::Potential)

The ``n_x \times n_x`` Green's matrix ``G`` of `potential`'s backend on the grid of `dist`, which
gives the field of a density in one product:

```math
\rho_i = h_v \sum_j f_{ij} , \qquad \phi = G \rho .
```

It is assembled once, at model construction, by solving the backend for the ``n_x`` unit vectors,
and symmetrised. ``\phi`` is linear in ``\rho``, so ``G \rho`` is the backend's
``\phi`` to round-off, and `G`'s symmetry makes ``\partial H / \partial f = h_x h_v h`` exact:
the discrete Hamiltonian ``H = h_x h_v \sum_{ij} f_{ij} v_j^2 / 2 + (h_x/2) \sum_i \rho_i \phi_i``
is quadratic in ``f``, and this is what lets an implicit midpoint step conserve it. `G` is
mean-free, the background being neutralising. `PoissonSolverMatrixFree` bounds its residual and
not its error, which may exceed that bound by the condition number of its operator — growing as
``n_x^2`` — so on a fine grid in `Float32` its cell solves are accurate to far fewer digits than
``eps(Float32)``, and `G` is then only as accurate as those solves.

The constructor of a grid `VlasovPoisson` calls this, so its three errors are raised there:

  - `ArgumentError` unless the potential's solver is a `PoissonSolverFFT` or a
    `PoissonSolverMatrixFree`. On a spline basis the density sampled at the nodes is not a
    Galerkin load vector, and the field would not be the one the scheme conserves with.
  - `DimensionMismatch` unless the potential has `size(dist.values, 1)` degrees of freedom.
  - `ArgumentError` unless the potential's own x-nodes are the distribution's. Both are laid out
    by `range(a, b, nx + 1)` over the same domain, so a matching grid compares exactly, and a
    grid built in another way is refused rather than run.

`dist.x` and `dist.v` also decide the bracket: `Arakawa(size(dist.values)..., step(dist.x),
step(dist.v))` is the one the run uses, and its own check refuses fewer than three nodes per
direction.
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

for a charge of ``+1`` and ``-\phi'' = \rho - \bar{\rho}``, with `G` the Green's matrix that
[`_green_matrix`](@ref) assembles. `params.arakawa` is the bracket, `params.v` the velocity
nodes and `params.G` that matrix; the state `f` is the node values of a `GridDistribution` in
the order of `vec(dist.values)`.

The bracket is doubly periodic, so the ``v``-index of the last node wraps onto the first with
the spacing ``h_v = (v_b - v_a)/(n_v - 1)``. **The scheme is periodic in `v`; it is valid while
`f` vanishes to `eps(T)·maximum(f)` at `v_a` and `v_b`**, where `T` is the element type of the
state: a distribution that reaches the ends of the velocity range of a `GridDistribution` is
then solved as the periodic problem, the wrap transporting mass from ``v_b`` to ``v_a``, while
the torus identities it is checked by still hold.
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

`method` is any `GeometricIntegrators` method for an `ODEProblem`, and it is the one thing a
caller chooses: `ImplicitMidpoint()` keeps ``H``, `ExplicitMidpoint()` is of the same second
order and does not, and `RK4()` is what `scripts/landau_damping_1d1v.jl` integrates with — the
implicit midpoint's dense Newton Jacobian is out of reach on that script's ``n_v``.

The scheme runs in `Float32` and `Float64`; `Float16` is not supported, since FFTW has no
`Float16` plan.

The vector field is [`_grid_1d1v_rhs!`](@ref), whose parameters carry a fresh `Arakawa` bracket
over `model.distribution`'s grid, its velocity nodes, and the Green's matrix `model.work`. The
scheme is **periodic in `v`**, and valid only while `f` vanishes to `eps(T)·maximum(f)` at both
ends of the velocity range: see [`_grid_1d1v_rhs!`](@ref).
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
