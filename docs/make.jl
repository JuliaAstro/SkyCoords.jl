using Documenter
using SkyCoords

DocMeta.setdocmeta!(SkyCoords, :DocTestSetup, :(using SkyCoords); recursive = true)

makedocs(
    sitename = "SkyCoords",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["Home" => "index.md", "API/Reference" => "api.md"],
    modules = [SkyCoords],
    strict = true,
)

deploydocs(repo = "github.com/JuliaAstro/SkyCoords.jl", push_preview=true)
