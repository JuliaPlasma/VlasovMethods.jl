# The `@allocated` barriers the electric-field and reduced-march test files share. It is included by
# `test/electric_field.jl`, `test/particles/electric_field.jl`, `test/particles/time_marching.jl`
# and `test/particles/reduced_equations.jl`, and is listed in no `runtests.jl`.

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

# The leapfrog step `leapfrog_step!(IC, efield, Δt, t)` and the gradient `∇H!(g, t, x, params)`.
function _alloc_step(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    (f(a, b, c, d); @allocated f(a, b, c, d))
end
function _alloc_gradient(f::F, a::A, b::B, c::C, d::D) where {F, A, B, C, D}
    (f(a, b, c, d); @allocated f(a, b, c, d))
end

# The save `save_solution(SS, IC, Ψ, efield, w, p, t, it, nsave, save)`, whose ten arguments are all
# concrete at the call site.
function _alloc_save(f::F, a::A, b::B, c::C, d::D, e::E, g::G, h::H, i::I, j::J,
        k::K) where {F, A, B, C, D, E, G, H, I, J, K}
    (f(a, b, c, d, e, g, h, i, j, k); @allocated f(a, b, c, d, e, g, h, i, j, k))
end
