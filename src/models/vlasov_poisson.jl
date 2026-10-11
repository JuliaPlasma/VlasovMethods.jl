
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

# The integrator steps its own copy of the particle state, so the field of a step is deposited
# from the positions in the first row of `z`, and not from `model.distribution`, which the
# integrator never writes. The weights come in as an argument for the same reason: the members of
# a training-data ensemble draw their own, so a deposit that read them off the model would be
# wrong for every member but one.
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

# Vector field for acceleration. It carries the model instead of taking it out of the parameters:
# a training-data ensemble's parameters are its members' weights alone, so that an
# `EnsembleSolution` round-trips through HDF5.
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

# Solution for Lorentz force
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

The energy of the Vlasov–Poisson model of `model`, as a callable invariant of signature
`h(t, q, params)`,

```math
h(t, q, w) = \frac{1}{2}\sum_i w_i v_i^2 + \frac{1}{2}\int (\partial_x\phi)^2 \,\mathrm{d}x ,
```

where ``\phi`` solves the Poisson problem whose right-hand side is a *fresh deposit of the
positions in `q`* with the weights `params.w`. The state it reads is the one it is given:
depositing `model.distribution` in its place would repeat the defect the field update of T3.2
fixed. It leaves the model's `potential` holding that deposit, as the acceleration substeps do,
and each of them deposits from its own state before it reads the field, so this cannot change a
trajectory.

The field term is `dot(coefficients(potential), rhs(potential)) / 2`. The Galerkin system of the
spline solver is ``S c = \text{rhs}`` for the stiffness matrix ``S``, so ``c^\top S c = c^\top
\text{rhs}`` is ``\frac{1}{2}\int(\partial_x\phi)^2`` without holding ``S`` and without
allocating.
"""
struct VlasovHamiltonian{MT}
    model::MT
end

function (h::VlasovHamiltonian)(t, q, params)
    model = h.model
    update_potential!(model, q, params.w)

    e = zero(eltype(q))
    w = params.w
    for i in axes(q, 2)
        e += w[1, i] * q[2, i]^2
    end

    return e / 2 +
           dot(PoissonSolvers.coefficients(model.potential),
        PoissonSolvers.rhs(model.potential)) / 2
end

# Constructor for a splitting method from GeometricIntegrators
# The problem is setup such that one solution step pushes all particles.
# While this allows for a simple implementation, it is not well-suited
# for parallelisation.
function SplittingMethod(
        model::VlasovPoisson{1, 1, <:ParticleDistribution}, tspan::Tuple, tstep::Real)
    # The model's own particles carry the weights of this path, so the parameters are the live
    # view of them and the deposit reads what it always read.
    params = (w = model.distribution.particles.w,)

    # create geometric problem
    equ = GeometricEquations.SODEProblem(
        (v_advection!, VAcceleration(model)),
        (s_advection!, SAcceleration(model)),
        tspan, tstep, copy(model.distribution.particles.z);
        parameters = params)

    # create integrator
    int = GeometricIntegrators.GeometricIntegrator(equ, GeometricIntegrators.Strang())

    # put together splitting method
    SplittingMethod(model, equ, int)
end

@doc raw"""
    sodeensemble(model, members, tspan, tstep)

An `EnsembleProblem` that generates the training data of `model`: one sample per member of
`members`, each member a `ParticleDistribution` whose particle count is the number of particles
the model's potential is deposited from.

Every member keeps the initial condition `(q = copy(m.particles.z),)` and its own weights,
`(w = copy(m.particles.w),)`. The weights are the deposit's, and an importance-sampled member
carries different ones from the next, so holding them in the shared model would be wrong for
every member but one; they are copied so that a later change to a member's distribution does not
reach the problem. The problem also carries the Hamiltonian `VlasovHamiltonian` as its invariant,
which no `SODEProblem` of this model had: `invariants(ens)` is then a `NamedTuple` rather than
`NullInvariants()`.

`GeometricIntegrators.integrate(ens, Strang())` solves every member into one `EnsembleSolution`,
whose per-member parameters are the members' weights; the HDF5 extension of `GeometricSolutions`
round-trips it with `h5save` and `h5load`.

A member whose weights and whose state name a different particle count throws an `ArgumentError`.
The model, its potential and the state are built in one element type; mixed input is not
supported.
"""
function sodeensemble(model::VlasovPoisson{1, 1, <:ParticleDistribution},
        members::AbstractVector{<:ParticleDistribution{T, 1, 1}},
        tspan::Tuple, tstep::Real) where {T}
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
