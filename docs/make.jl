using Documenter
using DocumenterCitations
using PlasmaWaves

DocMeta.setdocmeta!(PlasmaWaves, :DocTestSetup, :(using PlasmaWaves); recursive = true)

const bib = CitationBibliography(joinpath(@__DIR__, "PlasmaWaves.jl.bib"))

makedocs(
    sitename = "PlasmaWaves.jl",
    format = Documenter.HTML(),
    modules = [PlasmaWaves],
    pages = [
        "Home" => "index.md",
    ],
    checkdocs = :exports,
    plugins = [bib],
    doctest = true
)

deploydocs(
    repo = "github.com/JuliaSpacePhysics/PlasmaWaves.jl",
)
