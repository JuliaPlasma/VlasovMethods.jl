```@meta
CurrentModule = VlasovMethods
```

# VlasovMethods

Documentation for [VlasovMethods](https://github.com/JuliaPlasma/VlasovMethods.jl).

```@index
```

```@autodocs
Modules = [VlasovMethods]
Filter = t -> !(t in (VlasovMethods.build_spline_parity_factors,
    VlasovMethods.build_grid_parity_factors, VlasovMethods.tie_continuum_frontends!,
    VlasovMethods.dual_parity_modes, VlasovMethods.klimontovich_deposition,
    VlasovMethods.klimontovich_bridge_loss, VlasovMethods.kruskal_oberman_free_energy,
    VlasovMethods.Landau, VlasovMethods.kernel, VlasovMethods.compute_L!,
    VlasovMethods.compute_J!, VlasovMethods.compute_K!, VlasovMethods.PairwiseLandau,
    VlasovMethods.entropy_gradient!, VlasovMethods.apply_G!, VlasovMethods.DiscreteGradient,
    VlasovMethods.GonzalezGradient, VlasovMethods.AVFGradient,
    VlasovMethods.discrete_gradient!,
    VlasovMethods.discrete_gradient_step!, VlasovMethods.Picard_iterate_Landau_nls!,
    VlasovMethods.collisional_vectorfield!, VlasovMethods.plot_distribution,
    VlasovMethods.plot_particles, VlasovMethods.plot_series, VlasovMethods.plot_groups,
    VlasovMethods.plot_damping, VlasovMethods.plot_phase_space,
    VlasovMethods.plot_velocity_projections, VlasovMethods.plot_collision_distribution, VlasovMethods.plot_drift,
    VlasovMethods.plot_entropy, VlasovMethods.plot_error, VlasovMethods.plot_convergence,
    VlasovMethods.plot_cost))
```

The `compute_J!` method of the metriplectic Lenard–Bernstein model; the Landau method is on the
page of the Landau collision operator.

```@docs
compute_J!(::Any, ::Any, ::Any, ::SplineDistribution{T, XD, 1},
    ::MetriplecticLenardBernstein) where {T, XD}
```
