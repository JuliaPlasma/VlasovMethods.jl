@doc raw"""
    parity_split(f::AbstractMatrix, f0::AbstractVector)

The parity parts `(f⁺, f⁻)` of the `nx × nv` grid values `f` about the background `f0`, a vector
over the `nv` velocity nodes, by the reflection ``(R f)_{ij} = f_{i, n_v + 1 - j}``:

```math
f^- = \tfrac{1}{2} (f - R f) , \qquad f^+ = \tfrac{1}{2} (f + R f) - 1 \otimes f_0 .
```

On a velocity grid symmetric about 0, ``R`` maps ``v`` to ``-v``, so ``f^-`` is odd and, for an
even `f0`, ``f^+`` is even. For a floating-point `f`, both parts have its element type.
[`parity_merge`](@ref) inverts the split. An `f0` whose length is not `size(f, 2)` throws a
`DimensionMismatch`.
"""
function parity_split(f::AbstractMatrix, f0::AbstractVector)
    Base.require_one_based_indexing(f, f0)
    length(f0) == size(f, 2) || throw(DimensionMismatch(
        "parity_split takes a background of size(f, 2) = $(size(f, 2)) values; got $(length(f0))"))
    T = eltype(f)
    Rf = f[:, end:-1:1]
    # halved before the sum, so that the sum of two values up to floatmax(T) does not overflow
    fm = f ./ 2 .- Rf ./ 2
    fp = f ./ 2 .+ Rf ./ 2 .- T.(f0)'
    return fp, fm
end

"""
    parity_merge(f0::AbstractVector, fp::AbstractMatrix, fm::AbstractMatrix)

The grid values `1 ⊗ f0 + fp + fm` of the background `f0` and the parity parts `fp` and `fm` of
[`parity_split`](@ref), which this inverts. Parts of two sizes, or an `f0` whose length is not
`size(fp, 2)`, throw a `DimensionMismatch`.
"""
function parity_merge(f0::AbstractVector, fp::AbstractMatrix, fm::AbstractMatrix)
    size(fp) == size(fm) && length(f0) == size(fp, 2) || throw(DimensionMismatch(
        "parity_merge takes two parts of one size and a background of size(fp, 2) values; " *
        "got parts of size $(size(fp)) and $(size(fm)) and $(length(f0)) values"))
    return f0' .+ fp .+ fm
end

@doc raw"""
    ParitySplitting(model::VlasovPoisson{1, 1, <:GridDistribution}, tspan::Tuple, tstep::Real;
                    background = nothing, tol = 4eps(T), max_iterations = 10)

The parity Strang splitting of the grid 1D1V Vlasov–Poisson model `model` over `tspan` with the
step `tstep`, with `T` the element type of `model.distribution.values`; [`run!`](@ref) advances
it. On a velocity grid symmetric about 0, the state is ``f = 1 \otimes f_0 + f^+ + f^-``, with
``f_0`` an even background over ``v``, ``f^+`` even and ``f^-`` odd ([`parity_split`](@ref)). With
``P(f, h)`` the model's bracket term, ``h(g) = G \rho(g) \otimes 1 + 1 \otimes v^2 / 2``,
``\rho(g) = h_v \sum_j g_{\cdot j}``, ``G`` = `model.work` and ``F_0 = 1 \otimes f_0``, the
equations are

```math
\dot{f}^+ = P(f^-, h(F_0 + f^+)) , \qquad \dot{f}^- = P(F_0 + f^+, h(F_0 + f^+)) ,
```

and one step of ``\tau`` is

1. ``f^- \leftarrow f^- + \tfrac{\tau}{2} P(F_0 + f^+, h(F_0 + f^+))``, exact for the ``f^-``
   equation, whose right-hand side does not depend on ``f^-``;
2. ``f^+ \leftarrow y``, the implicit midpoint step of the ``f^+`` equation with ``f^-`` frozen,
   the solution of ``y - f^+ - \tau P(f^-, h(F_0 + (f^+ + y) / 2)) = 0``;
3. the half step of 1 again.

``P(f^-, \cdot)`` with ``f^-`` frozen is skew-gradient for the discrete energy
``H = h_x h_v \sum f v^2 / 2 + (h_x / 2) \rho^\mathsf{T} G \rho``, which depends on the even part
of ``f`` alone, so step 2 keeps ``H`` to the tolerance of its solve, and steps 1 and 3 do not
change it. The splitting is of second order and keeps the mass; it does not keep the ``L^2`` norm.

Step 2 is linear in ``y`` and is solved by GeometricSolvers' `Newton` with the matrix-free `GMRES`,
`FixedForcing(1e-3)`, and a `JVPFunction` of the derivative
``w \mapsto w - \tfrac{\tau}{2} P(f^-, (G \rho(w)) \otimes 1)``, from ``y = f^+``, to the absolute
residual test ``\|F\|_2 \le \mathrm{tol} \, \|F_0 + f^+\|_2`` in at most `max_iterations` Newton
steps; the relative residual and step tests are off. A solve that does not return `SUCCESS` makes
[`run!`](@ref) throw an `ErrorException` that names its status code, ``\|F\|`` and the iterations.

`background` is ``f_0``, a vector over the `nv` velocity nodes; `nothing` takes the x-mean of the
even part of the initial values,
``f_{0,j} = \frac{1}{n_x} \sum_i (f_{ij} + f_{i, n_v + 1 - j}) / 2``. ``f_0`` is replaced by
``(f_0 + \mathrm{reverse}(f_0)) / 2``. An x-dependent background is not stationary, so a matrix
is refused. The method stores a copy of the initial values, from which every run starts, and
`tspan`, `tstep` and `tol` converted to `T`. It runs on the periodic velocity closure of the
model's bracket term, in `Float32` and `Float64`, and an odd `nv` is accepted.

The constructor throws

  - an `ArgumentError` when the v-grid is not symmetric,
    ``|v_j + v_{n_v + 1 - j}| > 10 \, \mathrm{eps}(T) \max(|v_1|, |v_{n_v}|)`` for some ``j``;
  - an `ArgumentError` when ``f_0`` is not even,
    ``|f_{0,j} - f_{0, n_v + 1 - j}| > 10 \, \mathrm{eps}(T) \max |f_0|`` for some ``j``;
  - an `ArgumentError` when `background` is a matrix, or neither `nothing` nor a vector;
  - a `DimensionMismatch` when `background` does not hold `nv` values;
  - an `ArgumentError` when the initial values or `background` hold a `NaN` or an `Inf`;
  - an `ArgumentError` when `tstep` is not finite and positive in `T`, when `tspan` is not two
    finite times in `T` or runs backward, when `tol` is negative or `NaN`, and when
    `max_iterations` is negative or above `typemax(Int32)`.

The time grid `GeometricIntegrators.TimeSeries(t₀, t₁, τ)` throws an `InexactError` when its
number of steps does not fit an `Int`, as for `tspan = (0, 1)` and `tstep = floatmin(T)`.
"""
struct ParitySplitting{T, MT <: VlasovPoisson, AT <: Arakawa{T}} <: VlasovMethod
    model::MT
    tspan::Tuple{T, T}
    tstep::T
    background::Vector{T}
    initial::Matrix{T}
    tol::T
    max_iterations::Int
    arakawa::AT
end

function ParitySplitting(model::VlasovPoisson{1, 1, <:GridDistribution}, tspan::Tuple,
        tstep::Real; background = nothing,
        tol::Real = 4 * eps(eltype(model.distribution.values)),
        max_iterations::Integer = 10)
    dist = model.distribution
    values = dist.values
    T = eltype(values)
    nx, nv = size(values)

    τ = T(tstep)
    isfinite(τ) && τ > zero(T) || throw(ArgumentError(
        "the step tstep is not finite and positive in the element type $T of the distribution; " *
        "got tstep = $tstep"))
    length(tspan) == 2 && all(t -> isfinite(T(t)), tspan) || throw(ArgumentError(
        "the time span tspan is not two finite times in the element type $T; got $tspan"))
    t₀, t₁ = T(tspan[1]), T(tspan[2])
    # the time grid of `run!`, which refuses a span that does not point forward
    GeometricIntegrators.TimeSeries(t₀, t₁, τ)
    tol ≥ 0 || throw(ArgumentError("the tolerance tol is negative or NaN; got tol = $tol"))
    0 ≤ max_iterations ≤ typemax(Int32) || throw(ArgumentError(
        "max_iterations is not between 0 and typemax(Int32); got $max_iterations"))

    v = dist.v
    vmax = max(abs(v[1]), abs(v[nv]))
    all(j -> abs(v[j] + v[nv + 1 - j]) ≤ 10 * eps(T) * vmax, 1:nv) || throw(ArgumentError(
        "the parity splitting needs a velocity grid symmetric about 0; got v from $(v[1]) to " *
        "$(v[nv])"))

    f0 = _parity_background(values, background)
    all(isfinite, values) && all(isfinite, f0) || throw(ArgumentError(
        "the initial values and the background of the parity splitting must be finite"))
    all(j -> abs(f0[j] - f0[nv + 1 - j]) ≤ 10 * eps(T) * maximum(abs, f0), 1:nv) ||
        throw(ArgumentError("the background f0 of the parity splitting is not even in v"))
    f0 = f0 ./ 2 .+ reverse(f0) ./ 2

    arakawa = Arakawa(nx, nv, step(dist.x), step(v))
    return ParitySplitting(model, (t₀, t₁), τ, f0, Matrix{T}(values), T(tol),
        Int(max_iterations), arakawa)
end

# The background `f0` of `ParitySplitting` in the element type of `values`: the x-mean of the
# even part of `values` for `nothing`, the given vector, converted, otherwise.
function _parity_background(values::AbstractMatrix{T}, ::Nothing) where {T}
    nx = size(values, 1)
    return vec(sum(values ./ 2 .+ values[:, end:-1:1] ./ 2; dims = 1)) ./ nx
end

function _parity_background(values::AbstractMatrix{T}, f0::AbstractVector) where {T}
    length(f0) == size(values, 2) || throw(DimensionMismatch(
        "the background of the parity splitting takes nv = $(size(values, 2)) values; got " *
        "$(length(f0))"))
    return collect(T, f0)
end

function _parity_background(values, f0::AbstractMatrix)
    throw(ArgumentError(
        "the background f0 of the parity splitting is a function of v alone, a vector of nv " *
        "values; got a $(size(f0, 1)) × $(size(f0, 2)) matrix, and an x-dependent background " *
        "is not stationary"))
end

function _parity_background(values, f0)
    throw(ArgumentError(
        "the background of the parity splitting is nothing or a vector of nv values; got a " *
        "$(typeof(f0))"))
end

# The bracket term `P(g, h) = poisson_apply(Arakawa(nx, nv, hx, hv), g, h) / (hx hv)` of the grid
# model, on states stored as vectors of `nx * nv` values, x fastest.
_parity_bracket(a::Arakawa, g, h) = poisson_apply(a, g, h) ./ (a.hx * a.hv)

# The field part `G ρ(g) ⊗ 1` of the Hamiltonian, with `ρ(g) = hv Σⱼ g`, as a vector.
function _parity_field(m::ParitySplitting, g)
    a = m.arakawa
    ρ = a.hv .* vec(sum(reshape(g, a.nx, a.nv); dims = 2))
    return repeat(m.model.work * ρ, a.nv)
end

# the half step of f⁻: f⁻ + τ/2 P(F0 + f⁺, h(F0 + f⁺)), with the kinetic part `kinetic` of h
function _parity_fminus_step(m::ParitySplitting, kinetic, F0, fp, fm, τ)
    g = F0 .+ fp
    return fm .+ (τ / 2) .* _parity_bracket(m.arakawa, g, _parity_field(m, g) .+ kinetic)
end

# The residual F(y) = y - f⁺ - τ P(f⁻, h(F0 + (f⁺ + y)/2)) of the implicit midpoint step of f⁺
# with f⁻ frozen, written into `r`.
function _parity_fplus_residual!(r, y, m::ParitySplitting, kinetic, F0, fp, fm, τ)
    g = F0 .+ (fp .+ y) ./ 2
    r .= y .- fp .- τ .* _parity_bracket(m.arakawa, fm, _parity_field(m, g) .+ kinetic)
    return r
end

# The derivative w ↦ w - τ/2 P(f⁻, G ρ(w) ⊗ 1) of the residual, which is affine in y, written
# into `Jv`.
function _parity_fplus_jvp!(Jv, w, m::ParitySplitting, fm, τ)
    Jv .= w .- (τ / 2) .* _parity_bracket(m.arakawa, fm, _parity_field(m, w))
    return Jv
end

# The implicit midpoint step of f⁺ with f⁻ frozen: the solution y of
# `_parity_fplus_residual!` = 0.
function _parity_fplus_step(m::ParitySplitting{T}, kinetic, F0, fp, fm, τ) where {T}
    residual!(r, y, p) = _parity_fplus_residual!(r, y, m, kinetic, F0, fp, fm, τ)
    jvp!(Jv, y, w, p) = _parity_fplus_jvp!(Jv, w, m, fm, τ)
    method = GeometricSolvers.Newton(;
        linsolve = GeometricSolvers.GMRES(; forcing = GeometricSolvers.FixedForcing(1e-3)),
        jacobian = GeometricSolvers.JVPFunction(jvp!))
    options = GeometricSolvers.Options(T; f_abstol = m.tol * norm(F0 .+ fp), f_reltol = 0,
        x_reltol = 0, max_iterations = m.max_iterations)
    y, status = GeometricSolvers.solve(
        GeometricSolvers.NonlinearProblem(residual!), copy(fp),
        method; options)
    status.code === GeometricSolvers.SUCCESS || error(
        "the f⁺ solve of the parity splitting returned $(status.code), with ‖F‖ = " *
        "$(status.fnorm) after $(status.iterations) iterations")
    return y
end

# one step of τ: the half step of f⁻, the implicit midpoint step of f⁺, the half step of f⁻
function _parity_step(m::ParitySplitting, kinetic, F0, fp, fm, τ)
    fm = _parity_fminus_step(m, kinetic, F0, fp, fm, τ)
    fp = _parity_fplus_step(m, kinetic, F0, fp, fm, τ)
    fm = _parity_fminus_step(m, kinetic, F0, fp, fm, τ)
    return fp, fm
end

"""
    run!(method::ParitySplitting, h5file)

Advance the parity splitting `method` over its time span from the copy of the initial values
that it stores, write every state into the HDF5 file `h5file`, copy the final state into the
model's distribution, and return that distribution.

The time grid is `GeometricIntegrators.TimeSeries(t₀, t₁, τ)`, the grid of
`run!(::GeometricIntegrator, h5file)`: one time per state, the initial one included, ending at the
last whole step at or below `t₁`. The file holds the dataset `z` of the size `(nx * nv, ntime)`,
one column `vec(parity_merge(f0, f⁺, f⁻))` per state, and the dataset `t` of the `ntime` times,
both in the element type of the distribution. The file is created, or overwritten, and closed on
return, also when a step throws; the distribution is then left as it was. A progress bar is
shown for a run of more than five seconds.
"""
function run!(method::ParitySplitting{T}, h5file) where {T}
    (; model, tspan, tstep, background) = method
    nx, nv = size(method.initial)
    times = GeometricIntegrators.TimeSeries(tspan[1], tspan[2], tstep)
    kinetic = repeat(collect(model.distribution.v) .^ 2 ./ 2; inner = nx)
    F0 = repeat(background; inner = nx)
    fp, fm = parity_split(method.initial, background)
    fp, fm = vec(fp), vec(fm)
    # the state of the parts `fp` and `fm`, as a vector of `nx * nv` values
    state(fp, fm) = vec(parity_merge(background, reshape(fp, nx, nv), reshape(fm, nx, nv)))

    h5 = h5open(h5file, "w")
    try
        h5z = create_dataset(
            h5, "z", T, ((nx * nv, length(times)), (nx * nv, -1)), chunk = (nx * nv, 1))
        h5t = create_dataset(h5, "t", T, ((length(times),), (-1,)), chunk = (1,))
        h5z[:, 1] = state(fp, fm)
        h5t[1] = tspan[1]

        @showprogress 5 for n in 1:(length(times) - 1)
            fp, fm = _parity_step(method, kinetic, F0, fp, fm, tstep)
            h5z[:, n + 1] = state(fp, fm)
            h5t[n + 1] = times[n]
        end
    finally
        close(h5)
    end

    model.distribution.values .= reshape(state(fp, fm), nx, nv)
    return model.distribution
end
