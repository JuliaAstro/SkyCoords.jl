using Documenter
using SkyCoords

DocMeta.setdocmeta!(SkyCoords, :DocTestSetup, :(using SkyCoords); recursive = true)
include("pages.jl")
makedocs(
    modules = [SkyCoords],
    sitename = "SkyCoords.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://juliaastro.org/SkyCoords/stable/",
    ),
    pages = pages,
    )

deploydocs(
    repo = "github.com/JuliaAstro/SkyCoords.jl",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
