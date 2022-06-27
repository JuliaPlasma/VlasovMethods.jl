using VlasovParticleMethods
using Documenter

DocMeta.setdocmeta!(VlasovParticleMethods, :DocTestSetup, :(using VlasovParticleMethods); recursive=true)

makedocs(;
    modules=[VlasovParticleMethods],
    authors="Michael Kraus",
    repo="https://github.com/JuliaPlasma/VlasovParticleMethods.jl/blob/{commit}{path}#{line}",
    sitename="VlasovParticleMethods.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaPlasma.github.io/VlasovParticleMethods.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPlasma/VlasovParticleMethods.jl",
    devbranch="main",
)
