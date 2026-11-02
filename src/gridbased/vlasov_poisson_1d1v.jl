using GeometricBrackets: poisson_apply!

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

# The buffers of the grid right-hand side in the element type `T` of a state, on an `nx × nv`
# grid: the density `ρ`, the field `ϕ` and the Hamiltonian `h`. A cache belongs to one task.
struct Grid1D1VCache{T} <: Cache{T}
    ρ::Vector{T}
    ϕ::Vector{T}
    h::Vector{T}
end

function Grid1D1VCache{T}(nx::Integer, nv::Integer) where {T}
    Grid1D1VCache{T}(zeros(T, nx), zeros(T, nx), zeros(T, nx * nv))
end

Base.eltype(::Type{<:Grid1D1VCache{T}}) where {T} = T
Cache(AT, c::Grid1D1VCache) = Grid1D1VCache{AT}(length(c.ρ), length(c.h) ÷ length(c.ρ))
CacheType(AT, ::Grid1D1VCache) = Grid1D1VCache{AT}

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

`params` also holds buffers in the element type `T` of the distribution: `kinetic`, the
``v_j^2 / 2``, and `ρ`, `ϕ` and `h`. A state and a derivative of element type `T` take a method
that writes into these buffers and allocates nothing. A state of another element type, such as
the `ForwardDiff.Dual` numbers of the implicit midpoint's Jacobian, takes the buffers `ρ`, `ϕ`
and `h` of its own element type from `params.cache`, a `CacheDict` whose parent holds the
buffers of `T`: the first call of a type builds them, and a later call allocates nothing.
Either method throws a `DimensionMismatch` before it writes any buffer when `f` or `ḟ` does not
hold ``n_x n_v`` values, and an `ArgumentError` when `ḟ` shares memory with `f`. A `params`
without these buffers takes a method that allocates its intermediate arrays. The buffers belong
to one task: two concurrent calls on one `params` are unsafe, because they share them.
"""
function _grid_1d1v_rhs!(ḟ::AbstractVector{T}, t, f::AbstractVector{T},
        params::NamedTuple{(:arakawa, :v, :G, :kinetic, :ρ, :ϕ, :h, :cache),
            <:Tuple{Arakawa{T}, Vararg{Any}}}) where {T}
    (; arakawa, G, kinetic, ρ, ϕ, h) = params
    return _grid_1d1v_buffered!(ḟ, f, arakawa, G, kinetic, ρ, ϕ, h)
end

function _grid_1d1v_rhs!(ḟ::AbstractVector, t, f::AbstractVector{ST},
        params::NamedTuple{(:arakawa, :v, :G, :kinetic, :ρ, :ϕ, :h, :cache)}) where {ST}
    (; arakawa, G, kinetic) = params
    (; ρ, ϕ, h) = params.cache[ST]
    return _grid_1d1v_buffered!(ḟ, f, arakawa, G, kinetic, ρ, ϕ, h)
end

# The right-hand side of `_grid_1d1v_rhs!` on the buffers `ρ`, `ϕ` and `h`, of the element type
# of `f`.
function _grid_1d1v_buffered!(ḟ, f, arakawa, G, kinetic, ρ, ϕ, h)
    T = eltype(ρ)
    nx, nv = arakawa.nx, arakawa.nv
    length(f) == length(ḟ) == nx * nv || throw(DimensionMismatch(
        "_grid_1d1v_rhs! on a $nx × $nv grid takes a state f and a derivative ḟ of nx * nv " *
        "values each; got $(length(f)) and $(length(ḟ))"))

    fill!(ρ, zero(T))
    for j in 1:nv
        ρ .+= view(f, ((j - 1) * nx + 1):(j * nx))
    end
    ρ .*= arakawa.hv
    mul!(ϕ, G, ρ)
    for j in 1:nv
        view(h, ((j - 1) * nx + 1):(j * nx)) .= ϕ .+ kinetic[j]
    end
    poisson_apply!(ḟ, arakawa, f, h)
    ḟ ./= arakawa.hx * arakawa.hv

    return ḟ
end

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
    T = eltype(dist.values)
    nx, nv = size(dist.values)
    buffers = Grid1D1VCache{T}(nx, nv)
    params = (arakawa = Arakawa(nx, nv, step(dist.x), step(dist.v)),
        v = dist.v, G = model.work, kinetic = collect(dist.v) .^ 2 ./ 2,
        ρ = buffers.ρ, ϕ = buffers.ϕ, h = buffers.h, cache = CacheDict(buffers))

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
