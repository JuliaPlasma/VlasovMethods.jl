# The `@allocated` barriers the electric-field, reduced-march, model and moment test files share.
# It is included by `test/electric_field.jl`, `test/particles/electric_field.jl`,
# `test/particles/time_marching.jl`, `test/particles/reduced_equations.jl`,
# `test/models/vlasov_poisson.jl`, `test/models/vlasov_with_external_fields.jl` and
# `test/gridbased/moments.jl`, and is listed in no `runtests.jl`.

# `@allocated` through a fixed-arity barrier whose arguments have concrete types, so that a Julia
# 1.11 closure or splat boxing is not counted. Each call runs once before it is measured.
function _alloc_update(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    (f(a, b, c, d); @allocated f(a, b, c, d))
end
function _alloc_efield(f::F, a::A, b::B, c::C) where {F, A, B, C}
    (f(a, b, c); @allocated f(a, b, c))
end
function _alloc_energy(f::F, a::A) where {F, A}
    (f(a); @allocated f(a))
end

# A call of four arguments: the leapfrog step `leapfrog_step!(IC, efield, Δt, t)` and the gradient
# `∇H!(g, t, x, params)`.
function _alloc4(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    (f(a, b, c, d); @allocated f(a, b, c, d))
end

# A call of five arguments: the advection solution `s_advection!(z, t, z̄, t̄, params)`.
function _alloc5(f::F, a::A, b::B, c::C, d::D, e::E) where {F, A, B, C, D, E}
    (f(a, b, c, d, e); @allocated f(a, b, c, d, e))
end

# A call of six arguments: the moments `_apply_∫dv!(y, f, ci, li, h₁, h₂)` and `_apply_∫dvᵀ!`.
function _alloc6(f::F, a::A, b::B, c::C, d::D, e::E, g::G) where {F, A, B, C, D, E, G}
    (f(a, b, c, d, e, g); @allocated f(a, b, c, d, e, g))
end

# A call of seven arguments: the save `save_solution(SS, IC, Ψ, efield, p, t, ts)`, and the
# moments `_apply_∫vdv!(y, f, v, ci, li, h₁, h₂)`, `_apply_∫v²dv!`, `_apply_∫dv_μ!` and
# `_apply_∫dvᵀ_μ!`.
function _alloc_save(f::F, a::A, b::B, c::C, d::D, e::E, g::G,
        h::H) where {F, A, B, C, D, E, G, H}
    (f(a, b, c, d, e, g, h); @allocated f(a, b, c, d, e, g, h))
end
