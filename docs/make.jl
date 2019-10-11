using Documenter
using SkyCoords

makedocs(
    sitename = "SkyCoords",
    format = Documenter.HTML(),
    modules = [SkyCoords],
    # strict=true
)

deploydocs(
    repo = "github.com/juliaastro/skycoords.jl.git"
)
