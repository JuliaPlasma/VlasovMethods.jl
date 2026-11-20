# Reproducing the papers

Every numerical experiment of the papers has a driver script in `scripts/`. This part of the
manual says how to run the scripts, what each one writes, and what a run with the papers'
parameters gives.

| paper | page | scripts |
|:--|:--|:--|
| Symplectic model reduction methods for the Vlasov equation | [Symplectic model reduction](reduced_basis.md) | `bump_on_tail_*.jl`, `analytic_quadratic.jl` |
| PD-LB, particle discretisation of the Lenard–Bernstein operator | [The collision manuscripts](collisions.md) | `pd_lb_*.jl` |
| MP-LB, metriplectic particle methods for the Lenard–Bernstein operator | [The collision manuscripts](collisions.md) | `lenard_bernstein_metriplectic.jl` |
| MP-L, metriplectic particle methods for the Landau operator | [The collision manuscripts](collisions.md) | `landau_cases.jl`, `landau_compare.jl` |

The other scripts are examples of their own: `vlasov_poisson.jl` and `lenard_bernstein.jl` run the
models of [Getting started](../manual/getting_started.md) from the command line,
`landau_damping_1d1v.jl` fits the damping rate of the grid model, and `verify_conservation.jl`
measures the conservation properties of the conservative Lenard–Bernstein operator on the three
boundary conditions of the spline basis.

## The scripts environment

The scripts have their own environment, `scripts/Project.toml`, which takes VlasovMethods from the
checkout it is in. Instantiate it once, from the root of the repository:

```
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
```

## Options and output

Every script takes `--runs-dir` and `--results-dir`. A script writes its data, the HDF5 files, into
the first, and its figures and tables into the second. The defaults are `runs/` and `results/` at
the root of the repository, which git ignores. A directory that does not exist is created. A flag
that a script does not know stops it with an `ArgumentError`, which names the flags it takes.

```
julia --project=scripts scripts/bump_on_tail_training.jl --runs-dir /tmp/runs --results-dir /tmp/results
```

Every HDF5 file carries the root attributes `commit`, the commit of the checkout that wrote it, or
`"unknown"` outside a checkout, and `command`, the arguments of the script. The drivers of the
collision manuscripts also write the attribute `script` and one attribute per option.

## The gate

The flag `--quick` of the collision drivers replaces some of the manuscripts' parameters with
smaller ones. `scripts/gate.jl` runs every driver with it, each in a process of its own and into a
temporary directory, and checks that each exits 0 within 900 s, writes its figures and its HDF5
files, and records its provenance:

```
julia --startup-file=no --project=scripts scripts/gate.jl
```

## Run times

A run time on these pages is one run of the script in a new Julia process, after
`Pkg.instantiate()` and `Pkg.precompile()` of the scripts environment, on an Apple silicon Mac with
other jobs beside it. It includes the load of the packages and the compilation of the script's own
calls. Each page names the Julia version of its measurements.
