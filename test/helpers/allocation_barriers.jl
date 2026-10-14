# The `@allocated` barriers the two electric-field test files share. It is included by
# `test/electric_field.jl` and `test/particles/electric_field.jl`, and is listed in no
# `runtests.jl`.

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
