using Documenter
using SkyCoords

makedocs(
    sitename = "SkyCoords",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API/Reference" => "api.md"
    ],
    modules = [SkyCoords],
    # strict=true
)

deploydocs(
    repo = "github.com/juliaastro/skycoords.jl.git"
)
