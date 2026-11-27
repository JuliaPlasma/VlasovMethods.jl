```@meta
CurrentModule = VlasovMethods
```

# Distributions and projections

## Distribution functions

```@autodocs
Modules = [VlasovMethods]
Pages = ["distributions/particle_distribution.jl", "distributions/spline_distribution.jl",
    "distributions/grid_distribution.jl", "distributions/tensor_spline_distribution.jl",
    "distributions/maxwellian.jl",
    "entropies/collision_entropy.jl"]
```

## Initial conditions

```@autodocs
Modules = [VlasovMethods]
Pages = ["examples/example.jl", "examples/normal.jl", "examples/shiftednormalv.jl",
    "examples/uniform.jl", "examples/shifteduniform.jl", "examples/doublemaxwellian.jl",
    "examples/sum_maxwellian.jl", "examples/bump.jl", "examples/bumpontail.jl",
    "examples/anisotropic_bimaxwellian.jl", "examples/double_maxwellian_2d.jl",
    "examples/bkw_distribution.jl"]
```

## Projections

```@autodocs
Modules = [VlasovMethods]
Pages = ["projections/projection.jl", "projections/potential.jl",
    "projections/distribution.jl", "projections/density.jl"]
```
