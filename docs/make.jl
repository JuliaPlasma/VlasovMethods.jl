using VlasovMethods
using Documenter

DocMeta.setdocmeta!(VlasovMethods, :DocTestSetup, :(using VlasovMethods); recursive = true)

makedocs(;
    modules = [VlasovMethods],
    authors = "Michael Kraus",
    repo = Remotes.GitHub("JuliaPlasma", "VlasovMethods.jl"),
    sitename = "VlasovMethods.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuliaPlasma.github.io/VlasovMethods.jl",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting started" => "manual/getting_started.md",
            "Distribution functions" => "manual/distributions.md",
            "Initial conditions" => "manual/initial_conditions.md",
            "Projections" => "manual/projections.md",
            "Running a model" => "manual/methods.md",
            "Vlasov–Poisson" => "manual/vlasov_poisson.md",
            "External electromagnetic fields" => "manual/external_fields.md",
            "Lenard–Bernstein collisions" => "manual/lenard_bernstein.md",
            "Landau collisions" => "manual/landau.md",
            "Grid-based reduced tensors" => "manual/grid_tensors.md",
            "Reduced-order models" => "manual/reduced_order.md",
            "Multi-modal autoencoder frontends" => "manual/autoencoder.md",
            "Plotting" => "manual/plotting.md"
        ],
        "Reproducing the papers" => [
            "Overview" => "papers/index.md",
            "Symplectic model reduction" => "papers/reduced_basis.md",
            "The collision manuscripts" => "papers/collisions.md"
        ],
        "Library" => [
            "Overview" => "library/index.md",
            "Distributions and projections" => "library/distributions.md",
            "Models and methods" => "library/models.md",
            "Landau collisions" => "library/landau.md",
            "Reduction" => "library/reduction.md",
            "Plotting" => "library/plotting.md"
        ]
    ]
)

deploydocs(;
    repo = "github.com/JuliaPlasma/VlasovMethods.jl",
    devurl = "latest",
    devbranch = "main"
)
