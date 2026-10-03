
struct VlasovPoisson{XD, VD, DT <: DistributionFunction{<:Any, XD, VD}, PT <: Potential,
    WT <: AbstractVector} <: VlasovModel
    distribution::DT
    potential::PT

    # the buffer of the `local_width(basis)` basis values that the deposit and the field
    # evaluation fill, so that neither allocates
    work::WT

    function VlasovPoisson(dist::DistributionFunction{<:Any, XD, VD}, potential) where {
            XD, VD}
        work = _local_buffers(basis(potential), eltype(PoissonSolvers.rhs(potential)))
        new{XD, VD, typeof(dist), typeof(potential), typeof(work)}(dist, potential, work)
    end
end

function update_potential!(model::VlasovPoisson)
    projection!(model.potential, model.distribution)
    PoissonSolvers.update!(model.potential)
end

# The integrator steps its own copy of the particle state, so the field of a step is deposited
# from the positions in `z`, its first row, and not from `model.distribution`, which the
# integrator never writes.
function update_potential!(model::VlasovPoisson, z::AbstractMatrix)
    _deposit!(model.potential, model.work, z, model.distribution.particles.w)
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

# vector field
function lorentz_force!(ż, t, z, params)
    update_potential!(params.model)
    for i in axes(ż, 2)
        ż[1, i] = z[2, i]
        ż[2, i] = - params.ϕ(z[1, i], 1)
    end
end

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

# Vector field for acceleration
function v_acceleration!(ż, t, z, params)
    update_potential!(params.model, z)
    work = params.model.work
    for i in axes(ż, 2)
        ż[1, i] = 0
        ż[2, i] = _electric_field(params.ϕ, work, z[1, i])
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
function s_acceleration!(z, t, z̄, t̄, params)
    update_potential!(params.model, z̄)
    work = params.model.work
    for i in axes(z, 2)
        z[1, i] = z̄[1, i]
        z[2, i] = z̄[2, i] + (t - t̄) * _electric_field(params.ϕ, work, z̄[1, i])
    end
end

# Constructor for a splitting method from GeometricIntegrators
# The problem is setup such that one solution step pushes all particles.
# While this allows for a simple implementation, it is not well-suited
# for parallelisation.
function SplittingMethod(
        model::VlasovPoisson{1, 1, <:ParticleDistribution}, tspan::Tuple, tstep::Real)
    # collect parameters
    params = (ϕ = model.potential, model = model)

    # create geometric problem
    equ = GeometricEquations.SODEProblem(
        (v_advection!, v_acceleration!),
        (s_advection!, s_acceleration!),
        tspan, tstep, copy(model.distribution.particles.z);
        parameters = params)

    # create integrator
    int = GeometricIntegrators.GeometricIntegrator(equ, GeometricIntegrators.Strang())

    # put together splitting method
    SplittingMethod(model, equ, int)
end
