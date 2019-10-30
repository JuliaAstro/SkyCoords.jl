SkyCoords.jl
============

Basic astronomical coordinate systems in Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/SkyCoords.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliaastro.github.io/SkyCoords.jl/latest)
[![Build Status](https://img.shields.io/travis/JuliaAstro/SkyCoords.jl.svg)](https://travis-ci.org/JuliaAstro/SkyCoords.jl)
[![Build status](https://img.shields.io/appveyor/ci/kbarbary/skycoords-jl.svg?label=windows)](https://ci.appveyor.com/project/kbarbary/skycoords-jl/branch/master)
[![codecov](https://codecov.io/gh/JuliaAstro/SkyCoords.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/SkyCoords.jl)

## Installation

```julia
julia> Pkg.add("SkyCoords")
```

## Usage

There are currently three supported coordinate systems. The following
immutable types are used to represent coordinates in each system:

- `ICRSCoords`: ICRS coordinates system
- `GalCoords`: Galactic coordinates system
- `FK5Coords`: FK5 coordinates system (with arbitrary equninox)

```julia
julia> c1 = ICRSCoords(0., 0.)  # inputs are ra, dec in radians
ICRSCoords{Float64}(0.0, 0.0)

julia> c2 = convert(GalCoords, c1) # convert to a different system
GalCoords{Float64}(1.6814027872278692, -1.0504884034813007)

julia> convert(FK5Coords{2000}, c1)
FK5Coords{2000,Float64}(1.1102233723050067e-7, 4.411803426976326e-8)

julia> separation(c1, ICRSCoords(1., 0.)) # radians
1.0

```

For more information, visit the [documentation](https://juliaastro.github.io/skycoords.jl/stable)

## License and Credits

License is MIT. This package profits from the hard work that went into
astropy.coordinates, *especially* in terms of testing and coordinate system
definitions.
