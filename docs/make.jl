using VlasovMethods
using Documenter

DocMeta.setdocmeta!(VlasovMethods, :DocTestSetup, :(using VlasovMethods); recursive=true)

makedocs(;
    modules=[VlasovMethods],
    authors="Michael Kraus",
    repo="https://github.com/JuliaPlasma/VlasovMethods.jl/blob/{commit}{path}#{line}",
    sitename="VlasovMethods.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaPlasma.github.io/VlasovMethods.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPlasma/VlasovMethods.jl",
    devbranch="main",
)
