using Documenter
using DocumenterInterLinks
using NearestNeighbors
using SkyCoords

DocMeta.setdocmeta!(SkyCoords, :DocTestSetup, :(using SkyCoords); recursive = true)
include("pages.jl")

# TODO: replace with upstream when merged
links = InterLinks(
    "SOFA" => "https://juliaastro.org/SOFA.jl/previews/PR38/",
)

makedocs(;
    modules = [Base.get_extension(SkyCoords, :NearestNeighborsExt), SkyCoords],
    plugins = [links],
    sitename = "SkyCoords.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://juliaastro.org/SkyCoords/stable/",
    ),
    pages,
    warnonly = [:linkcheck],
)

deploydocs(
    repo = "github.com/JuliaAstro/SkyCoords.jl",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
