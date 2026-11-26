```@meta
CurrentModule = VlasovMethods
```

# Models and methods

## Methods

```@autodocs
Modules = [VlasovMethods]
Pages = ["methods/splitting.jl", "methods/geometric_integrator.jl"]
```

## Vlasov–Poisson

```@autodocs
Modules = [VlasovMethods]
Pages = ["models/vlasov_poisson.jl", "gridbased/vlasov_poisson_1d1v.jl", "gridbased/parity.jl",
    "splinebased/vlasov_poisson_1d1v.jl"]
```

## External electromagnetic fields

```@autodocs
Modules = [VlasovMethods]
Pages = ["models/vlasov_with_external_fields.jl"]
```

## Lenard–Bernstein collisions

```@autodocs
Modules = [VlasovMethods]
Pages = ["models/lenard_bernstein.jl", "models/lenard_bernstein_conservative.jl",
    "models/rescaled_lenard_bernstein_conservative.jl",
    "models/lenard_bernstein_metriplectic.jl"]
```
