
struct VlasovPoisson{XD, VD, DT <: DistributionFunction{<:Any, XD, VD}, PT <: Potential,
    WT <: AbstractVector} <: VlasovModel
    distribution::DT
    potential::PT

    # the buffer of the `local_width(basis)` basis values that the deposit and the field
    # evaluation fill, so that neither allocates
    work::WT

    function VlasovPoisson(dist::DistributionFunction{<:Any, XD, VD}, potential) where {
            XD, VD}
        work = _deposit_buffer(potential)
        new{XD, VD, typeof(dist), typeof(potential), typeof(work)}(dist, potential, work)
    end
end

function update_potential!(model::VlasovPoisson)
    projection!(model.potential, model.distribution)
    PoissonSolvers.update!(model.potential)
end

# The field of a step is deposited from the positions in the first row of `z`, with the weights
# `w`, and not from `model.distribution`: the integrator steps its own copy of the state, and each
# member of a training-data ensemble carries its own weights.
function update_potential!(model::VlasovPoisson, z::AbstractMatrix, w::AbstractMatrix)
    _deposit!(model.potential, model.work, z, w)
    PoissonSolvers.update!(model.potential)
end

# The electric field `-∂ₓϕ` at `x`, summed over the basis functions that do not vanish there,
# with `work` as the buffer for their derivatives. `evaluate_all!` reduces a periodic argument
# into the domain, so an unwrapped position needs no wrap here.
function _electric_field(potential::Potential, work::AbstractVector, x::Number)
    b = basis(potential)
    c = PoissonSolvers.coefficients(potential)
    j₀ = evaluate_all!(work, b, x, 1)
    ∂ϕ = zero(eltype(c))
    for (t, value) in pairs(work)
        # the padding of a recombined basis, as in `projection!`
        iszero(value) && continue
        ∂ϕ += c[basis_index(b, j₀ + t - 1)] * value
    end
    return -∂ϕ
end

####################################################
# Define Splitting Method for Vlasov-Poisson Model #
####################################################

###########################################################
# Vlasov-Poisson 1D1V splitting fields for particles      #
###########################################################

# Vector field for advection
function v_advection!(ż, t, z, params)
    for i in axes(ż, 2)
        ż[1, i] = z[2, i]
        ż[2, i] = 0
    end
end

# Vector field for acceleration. It holds the model, so that the parameters are the weights
# alone, which HDF5 can store.
struct VAcceleration{MT}
    model::MT
end

function (v::VAcceleration)(ż, t, z, params)
    update_potential!(v.model, z, params.w)
    work = v.model.work
    for i in axes(ż, 2)
        ż[1, i] = 0
        ż[2, i] = _electric_field(v.model.potential, work, z[1, i])
    end
end

# Solution for advection
function s_advection!(z, t, z̄, t̄, params)
    for i in axes(z, 2)
        z[1, i] = z̄[1, i] + (t-t̄) * z̄[2, i]
        z[2, i] = z̄[2, i]
    end
end

# Solution for acceleration
struct SAcceleration{MT}
    model::MT
end

function (s::SAcceleration)(z, t, z̄, t̄, params)
    update_potential!(s.model, z̄, params.w)
    work = s.model.work
    for i in axes(z, 2)
        z[1, i] = z̄[1, i]
        z[2, i] = z̄[2, i] + (t - t̄) * _electric_field(s.model.potential, work, z̄[1, i])
    end
end

@doc raw"""
    VlasovHamiltonian(model)

The energy of the Vlasov–Poisson model `model`, as an invariant `h(t, q, params)`,

```math
h(t, q, w) = \frac{1}{2}\sum_i w_i v_i^2 + \frac{1}{2}\int (\partial_x\phi)^2 \,\mathrm{d}x ,
```

where ``\phi`` solves the Poisson problem for a fresh deposit of the positions in `q` with the
weights `params.w`. The call leaves that deposit in `model.potential`. Each acceleration substep
deposits its own state before it reads the field, so this does not change a trajectory.

The field term is the Gauss quadrature of ``(\partial_x\phi)^2`` on the basis of the potential,
tabulated once by the constructor, so a call allocates nothing. Its terms are non-negative, so the
sum keeps their relative accuracy also for a small field.
"""
struct VlasovHamiltonian{MT, XT <: AbstractVector, WT <: AbstractVector}
    model::MT
    nodes::XT
    weights::WT
end

function VlasovHamiltonian(model::VlasovPoisson)
    quadrature = SplineQuadrature(basis(model.potential))
    return VlasovHamiltonian(
        model, quadrature_nodes(quadrature), quadrature_weights(quadrature))
end

function (h::VlasovHamiltonian)(t, q, params)
    model = h.model
    update_potential!(model, q, params.w)

    e = zero(eltype(q))
    w = params.w
    for i in axes(q, 2)
        e += w[1, i] * q[2, i]^2
    end

    # `_electric_field` is `-∂ₓϕ`; its square is the field energy's integrand.
    f = zero(eltype(q))
    work = model.work
    for (x, ω) in zip(h.nodes, h.weights)
        ∂ϕ = _electric_field(model.potential, work, x)
        f += ω * ∂ϕ^2
    end

    return e / 2 + f / 2
end

function SplittingMethod(
        model::VlasovPoisson{1, 1, <:ParticleDistribution}, tspan::Tuple, tstep::Real)
    # the live view of the model's own weights, so that the deposit sees a change to them
    params = (w = model.distribution.particles.w,)
    _splitting_method(model, (v_advection!, VAcceleration(model)),
        (s_advection!, SAcceleration(model)), params, tspan, tstep, GeometricIntegrators.Strang())
end

@doc raw"""
    sodeensemble(model, members, tspan, tstep)

The training data of `model` as a `GeometricEquations.SODEEnsemble`, one sample per
`ParticleDistribution` in `members`.

The initial condition of a member `m` is `(q = copy(m.particles.z),)`, and its parameters are its
own weights, `(w = copy(m.particles.w),)`, because an importance-sampled draw gives each member
different weights. Both are copies, so a later change to `m` does not reach the problem. The
invariant `h` is the `VlasovHamiltonian` of `model`.

`GeometricIntegrators.integrate(ens, Strang())` solves every member into one `EnsembleSolution`,
which `h5save` and `h5load` round-trip through the HDF5 extension of `GeometricSolutions`. The
invariant is evaluated along a stored solution once per saved state, the initial state included,
with

```julia
compute_invariant(sol.t, sol[i][:q], parameters(ens)[i], invariants(equation(ens)).h)
```

`integrate` itself evaluates no invariant.

All members share the potential and the work buffer of `model`. Solve them one after another, as
`integrate` does: members solved concurrently read each other's fields, with no error.

An empty member list, and a member whose weights and state have different particle counts, throw
an `ArgumentError`. The model and the members have one element type `T`; a member of another
element type raises a `MethodError`.
"""
function sodeensemble(model::VlasovPoisson{1, 1, <:ParticleDistribution{T}},
        members::AbstractVector{<:ParticleDistribution{T, 1, 1}},
        tspan::Tuple, tstep::Real) where {T}
    isempty(members) && throw(ArgumentError("sodeensemble needs at least one member"))

    ics = [(q = copy(m.particles.z),) for m in members]

    params = map(members) do m
        nw = size(m.particles.w, 2)
        nz = size(m.particles.z, 2)
        nw == nz || throw(ArgumentError(
            "the weights of a member have $nw particles, its state has $nz"))
        (w = copy(m.particles.w),)
    end

    return GeometricEquations.SODEEnsemble(
        (v_advection!, VAcceleration(model)),
        (s_advection!, SAcceleration(model)),
        tspan, tstep, ics;
        parameters = params,
        invariants = (h = VlasovHamiltonian(model),))
end
