using Documenter
using PlasmaWaves

DocMeta.setdocmeta!(PlasmaWaves, :DocTestSetup, :(using PlasmaWaves); recursive = true)

makedocs(
    sitename = "PlasmaWaves.jl",
    format = Documenter.HTML(),
    modules = [PlasmaWaves],
    pages = [
        "Home" => "index.md",
    ],
    checkdocs = :exports,
    doctest = true
)

deploydocs(
    repo = "github.com/JuliaSpacePhysics/PlasmaWaves.jl",
)
