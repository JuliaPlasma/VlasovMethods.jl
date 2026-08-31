# VlasovMethods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPlasma.github.io/VlasovMethods.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaPlasma.github.io/VlasovMethods.jl/latest/)
[![Build Status](https://github.com/JuliaPlasma/VlasovMethods.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaPlasma/VlasovMethods.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaPlasma/VlasovMethods.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaPlasma/VlasovMethods.jl)

This package implements various numerical methods for the Vlasov-Equation, including particle methods, grid-based methods, and reduced basis methods.


## Development

### Git hooks

Two hooks live in `.githooks`. They are **not active in a fresh clone** — `core.hooksPath` is local
configuration and does not travel with a push — so enable them once per clone:

```sh
git config core.hooksPath .githooks
```

**`pre-commit`** acts on **staged `.jl` files only**, and exits immediately when a commit stages
none, so a documentation- or workflow-only commit is not slowed down by it:

- **JuliaFormatter `--check`**, honouring this repository's own `.JuliaFormatter.toml` — **blocks**
  the commit. Formatting is mechanical and always fixable.
- **`fatou lint`**, when `fatou` is installed — **advisory only**, and deliberately so: its
  `unused-import` rule does not follow `include`, so it flags the load-bearing imports of every
  module file.
- **`using <Package>`**, which catches a syntax error or a broken `include` — **blocks**. Note that
  this package does not currently load under Julia 1.13: `RungeKutta` 0.5 pulls in
  `GenericLinearAlgebra`, whose `eigencopy_oftype` method for `UpperHessenberg` overwrites the one
  Julia 1.13 added. Until the dependency is bumped, this check will block a commit that stages a
  `.jl` file.

**`pre-push`** runs the full test suite with `--check-bounds=auto`, but **only when pushing to
`main` or `master`**; a topic branch is left to CI. It prints nothing for **10–30 minutes**, which
looks exactly like a network hang and is not one. If you do interrupt it, check for an orphaned
Julia process that the killed hook left behind.

Either hook can be bypassed for a single command with `--no-verify`, for a change you know it does
not apply to:

```sh
git commit --no-verify
git push --no-verify
```

The hooks are generated from one shared copy and are byte-identical across the related
repositories, so edit them there rather than here — a local edit is silently undone by the next
install.
