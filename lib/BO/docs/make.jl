using BO
using Documenter

DocMeta.setdocmeta!(BO, :DocTestSetup, :(using BO); recursive=true)

makedocs(;
    modules=[BO],
    authors="Beforerr <zzj956959688@gmail.com> and contributors",
    sitename="BO.jl",
    format=Documenter.HTML(;
        canonical="https://Beforerr.github.io/BO.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
    ],
)

deploydocs(;
    repo="github.com/Beforerr/BO.jl",
    devbranch="main",
)
