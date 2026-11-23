using GeometricBrackets: poisson_apply!

@doc raw"""
    VelocityClosure

The closure of the velocity grid of the grid 1D1V Vlasov–Poisson model: what the Arakawa bracket
reads at and beyond the end nodes ``v_1 = v_a`` and ``v_{n_v} = v_b``. The two closures are
[`PeriodicVelocityClosure`](@ref) and [`BoundedVelocityClosure`](@ref). A closure is the
`closure` keyword of the model's [`GeometricIntegrator`](@ref) and of
[`VlasovMethods.LowRankFields1D1V`](@ref).
"""
abstract type VelocityClosure end

@doc raw"""
    PeriodicVelocityClosure()

The periodic closure of the velocity grid, the default: the Arakawa bracket wraps the velocity
direction, so the node beyond ``v_b`` is ``v_a``. The unknowns are all ``n_v`` nodes. A run is valid
only while ``f`` vanishes at both ends of the velocity range.
"""
struct PeriodicVelocityClosure <: VelocityClosure end

@doc raw"""
    BoundedVelocityClosure()

The bounded closure of the velocity grid, with zero inflow. The unknowns are the nodes
``j = 2, …, n_v - 1``, and ``f`` and the Hamiltonian are read as zero on and beyond the end nodes
``j = 1`` and ``j = n_v``, with ``\chi`` the indicator of the unknowns. The field is the periodic
bracket of the cut data, read at the unknowns,

```math
\dot{f} = \chi \circ [P h, \chi \circ f] , \qquad
P h = \chi \circ (\phi \otimes 1 + 1 \otimes \tfrac{v^2}{2}) , \qquad
\phi = G \rho , \qquad \rho_i = h_v \sum_{j=2}^{n_v-1} f_{ij} ,
```

and zero at the end nodes. On grid functions that vanish on and beyond the end nodes, Arakawa's
form ``T(a, b, c) = \langle a, [b, c] \rangle`` is totally antisymmetric, so ``\dot{f}`` is the
field of the cyclic average of ``T``, and ``\langle P h, \dot{f} \rangle = 0``: the discrete
energy ``H = h_x h_v \sum_i \sum_{j=2}^{n_v-1} f_{ij} v_j^2 / 2 + (h_x / 2) \rho^\mathsf{T} G \rho``
is kept by the field, and by the implicit midpoint rule. The mass is not kept: its rate is the
flux ``-h_x h_v \sum_i \sum_{j \in \{1, n_v\}} [P h, \chi \circ f]_{ij}`` through the end nodes.
Where ``f`` vanishes on the two outermost nodes at each end, ``j \in \{1, 2, n_v - 1, n_v\}``,
the two closures give the same field.
"""
struct BoundedVelocityClosure <: VelocityClosure end

# The guard of the `closure` keyword of `fname`: a `TypeError` for a subtype of `VelocityClosure`
# other than the two closures above, which would otherwise run as the periodic closure.
function _check_closure(fname, closure)
    closures = Union{PeriodicVelocityClosure, BoundedVelocityClosure}
    closure isa closures || throw(TypeError(fname, "keyword closure", closures, closure))
    return nothing
end

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
# grid: the density `ρ`, the field `ϕ`, the Hamiltonian `h`, and `g`, the state with its end
# columns set to zero, which only the bounded closure uses: it has `ng` values, `nx * nv` for the
# bounded closure and 0 for the periodic one. A cache belongs to one task.
struct Grid1D1VCache{T} <: Cache{T}
    ρ::Vector{T}
    ϕ::Vector{T}
    h::Vector{T}
    g::Vector{T}
end

function Grid1D1VCache{T}(nx::Integer, nv::Integer, ng::Integer = 0) where {T}
    Grid1D1VCache{T}(zeros(T, nx), zeros(T, nx), zeros(T, nx * nv), zeros(T, ng))
end

Base.eltype(::Type{<:Grid1D1VCache{T}}) where {T} = T
function Cache(AT, c::Grid1D1VCache)
    Grid1D1VCache{AT}(length(c.ρ), length(c.h) ÷ length(c.ρ), length(c.g))
end
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

# the guard on the lengths of the state and the derivative of both closures
function _check_grid_1d1v(ḟ, f, nx, nv)
    length(f) == length(ḟ) == nx * nv || throw(DimensionMismatch(
        "_grid_1d1v_rhs! on a $nx × $nv grid takes a state f and a derivative ḟ of nx * nv " *
        "values each; got $(length(f)) and $(length(ḟ))"))
    return nothing
end

# The right-hand side of `_grid_1d1v_rhs!` on the buffers `ρ`, `ϕ` and `h`, of the element type
# of `f`.
function _grid_1d1v_buffered!(ḟ, f, arakawa, G, kinetic, ρ, ϕ, h)
    T = eltype(ρ)
    nx, nv = arakawa.nx, arakawa.nv
    _check_grid_1d1v(ḟ, f, nx, nv)

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
    _grid_1d1v_bounded_rhs!(ḟ, t, f, params)

The right-hand side of the grid 1D1V Vlasov–Poisson model with the
[`BoundedVelocityClosure`](@ref), as the vector field of a `GeometricEquations.ODEProblem`:

```math
\dot{f} = \chi \circ [P h, \chi \circ f] , \qquad
P h = \chi \circ (\phi \otimes 1 + 1 \otimes \tfrac{v^2}{2}) , \qquad
\phi = G \rho , \qquad \rho_i = h_v \sum_{j=2}^{n_v-1} f_{ij} ,
```

with ``\chi`` the indicator of the unknowns ``j = 2, …, n_v - 1`` and ``[\cdot, \cdot]`` the doubly
periodic bracket of [`_grid_1d1v_rhs!`](@ref). It writes zero in columns ``1`` and ``n_v`` of
`ḟ`, and reads no value of `f` in those columns: it copies the unknowns of `f` into the buffer
`g` of the cache, whose end columns are zero, so a non-finite value there does not reach `ḟ`.

`params` is that of [`_grid_1d1v_rhs!`](@ref), and its cache holds `g` too. A state
of the element type `T` of the buffers allocates nothing; a state of another element type, such
as the `ForwardDiff.Dual` numbers of the implicit midpoint's Jacobian, takes the buffers of its
own element type from `params.cache`. It throws the errors of [`_grid_1d1v_rhs!`](@ref), with
its messages, before it writes: a `DimensionMismatch` when `f` or `ḟ` does not hold ``n_x n_v``
values, and an `ArgumentError` when `ḟ` shares memory with `f`. At ``n_v = 3`` the one unknown
row has no velocity neighbour among the unknowns, and the field is zero.
"""
function _grid_1d1v_bounded_rhs!(ḟ::AbstractVector{T}, t, f::AbstractVector{T},
        params::NamedTuple{(:arakawa, :v, :G, :kinetic, :ρ, :ϕ, :h, :cache),
            <:Tuple{Arakawa{T}, Vararg{Any}}}) where {T}
    (; arakawa, G, kinetic, ρ, ϕ, h) = params
    g = parent(params.cache).g
    return _grid_1d1v_bounded_buffered!(ḟ, f, arakawa, G, kinetic, ρ, ϕ, h, g)
end

function _grid_1d1v_bounded_rhs!(ḟ::AbstractVector, t, f::AbstractVector{ST},
        params::NamedTuple{(:arakawa, :v, :G, :kinetic, :ρ, :ϕ, :h, :cache)}) where {ST}
    (; arakawa, G, kinetic) = params
    (; ρ, ϕ, h, g) = params.cache[ST]
    return _grid_1d1v_bounded_buffered!(ḟ, f, arakawa, G, kinetic, ρ, ϕ, h, g)
end

# The right-hand side of `_grid_1d1v_bounded_rhs!` on the buffers `ρ`, `ϕ`, `h` and `g`, of the
# element type of `f`. The unknowns are copied into `g` and its end columns are set to zero, not
# multiplied by zero, so that a NaN there is not carried into the bracket.
function _grid_1d1v_bounded_buffered!(ḟ, f, arakawa, G, kinetic, ρ, ϕ, h, g)
    T = eltype(ρ)
    nx, nv = arakawa.nx, arakawa.nv
    _check_grid_1d1v(ḟ, f, nx, nv)
    Base.mightalias(ḟ, f) && throw(ArgumentError(
        "Pc shares memory with û or c, which the stencil reads while it writes Pc"))
    column(j) = ((j - 1) * nx + 1):(j * nx)
    unknowns = (nx + 1):((nv - 1) * nx)

    fill!(g, zero(T))
    copyto!(g, first(unknowns), f, first(unknowns), length(unknowns))

    fill!(ρ, zero(T))
    for j in 2:(nv - 1)
        ρ .+= view(g, column(j))
    end
    ρ .*= arakawa.hv
    mul!(ϕ, G, ρ)
    fill!(h, zero(T))
    for j in 2:(nv - 1)
        view(h, column(j)) .= ϕ .+ kinetic[j]
    end
    poisson_apply!(ḟ, arakawa, g, h)
    ḟ ./= arakawa.hx * arakawa.hv
    fill!(view(ḟ, column(1)), zero(eltype(ḟ)))
    fill!(view(ḟ, column(nv)), zero(eltype(ḟ)))

    return ḟ
end

@doc raw"""
    GeometricIntegrator(model::VlasovPoisson{1, 1, <:GridDistribution}, tspan, tstep,
                        method = Integrators.ImplicitMidpoint();
                        closure::VelocityClosure = PeriodicVelocityClosure())

The `GeometricIntegrator` that steps the grid 1D1V Vlasov–Poisson model `model` over `tspan`
with step `tstep`, in the state `vec(model.distribution.values)`.

`method` is any `GeometricIntegrators` method for an `ODEProblem`. `ImplicitMidpoint()` keeps
``H``; `RK4()` is the one `scripts/landau_damping_1d1v.jl` uses, because the implicit midpoint's
dense Newton Jacobian is too large on that script's ``n_v``. The scheme runs in `Float32` and
`Float64`; FFTW has no `Float16` plan.

`closure` is the closure of the velocity grid. With [`PeriodicVelocityClosure`](@ref), the
default, the vector field is [`_grid_1d1v_rhs!`](@ref). With
[`BoundedVelocityClosure`](@ref), the vector field is [`_grid_1d1v_bounded_rhs!`](@ref), and the
initial state is `vec(model.distribution.values)` with its columns ``1`` and ``n_v`` set to zero;
`model.distribution` is not changed until `run!` writes the states back. A `closure` of another
type throws a `TypeError`.
"""
function GeometricIntegrator(model::VlasovPoisson{1, 1, <:GridDistribution}, tspan::Tuple,
        tstep::Real, method = Integrators.ImplicitMidpoint();
        closure::VelocityClosure = PeriodicVelocityClosure())
    _check_closure(:GeometricIntegrator, closure)
    dist = model.distribution
    T = eltype(dist.values)
    nx, nv = size(dist.values)
    bounded = closure isa BoundedVelocityClosure
    buffers = Grid1D1VCache{T}(nx, nv, bounded ? nx * nv : 0)
    params = (arakawa = Arakawa(nx, nv, step(dist.x), step(dist.v)),
        v = dist.v, G = model.work, kinetic = collect(dist.v) .^ 2 ./ 2,
        ρ = buffers.ρ, ϕ = buffers.ϕ, h = buffers.h, cache = CacheDict(buffers))

    q₀ = copy(vec(dist.values))
    if bounded
        q₀[1:nx] .= zero(T)
        q₀[((nv - 1) * nx + 1):end] .= zero(T)
    end
    rhs! = bounded ? _grid_1d1v_bounded_rhs! : _grid_1d1v_rhs!

    equ = GeometricEquations.ODEProblem(rhs!, tspan, tstep, q₀; parameters = params)

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
