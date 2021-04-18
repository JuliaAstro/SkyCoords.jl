using Documenter
using SkyCoords

DocMeta.setdocmeta!(SkyCoords, :DocTestSetup, :(using SkyCoords); recursive = true)

makedocs(
    modules = [SkyCoords],
    sitename = "SkyCoords.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://juliaastro.github.io/SkyCoords.jl",
    ),
    pages = ["Home" => "index.md", "API/Reference" => "api.md"],
    strict = true,
)

deploydocs(repo = "github.com/JuliaAstro/SkyCoords.jl", push_preview=true)
