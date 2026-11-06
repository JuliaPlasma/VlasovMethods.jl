```@meta
CurrentModule = VlasovMethods
```

# The Landau collision operator

The package implements the two particle methods of the Landau manuscript for two velocity
dimensions.

  - `Landau(dist, ent; γ = -3, ν = 1.0)` is the appendix method. Its vector field is the
    least-squares solution of a system in the matrices ``\mathbb{K}`` of [`compute_K!`](@ref),
    with the right-hand side ``\mathbb{L}`` of [`compute_L!`](@ref) times ``\mathbb{L}_k`` of
    [`compute_J!`](@ref). [`Picard_iterate_Landau_nls!`](@ref) integrates it in time with the
    implicit midpoint rule.
  - [`PairwiseLandau`](@ref) is the particle-pair method of the main text. Its vector field is
    the product of the matrix ``\mathbb{G}`` of [`apply_G!`](@ref) with the gradient of the
    entropy of [`entropy_gradient!`](@ref). [`discrete_gradient_step!`](@ref) integrates it in
    time with a discrete gradient, by default the [`GonzalezGradient`](@ref).

Both methods take the kernel ``U`` of [`kernel`](@ref), with the exponent `γ`.

## The appendix method

```@docs
kernel
compute_K!
compute_L!
compute_J!(::Any, ::SplineDistribution{T, XD, VD}, ::Landau) where {T, XD, VD}
Picard_iterate_Landau_nls!
```

## The particle-pair method

```@docs
PairwiseLandau
entropy_gradient!
apply_G!
collisional_vectorfield!(::AbstractArray{ST}, ::AbstractArray{ST}, ::Any,
    ::PairwiseLandau) where {ST}
DiscreteGradient
GonzalezGradient
discrete_gradient!
discrete_gradient_step!
```
