# SkyCoords.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/JuliaAstro/SkyCoords.jl)
[![CI](https://github.com/JuliaAstro/SkyCoords.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaAstro/SkyCoords.jl/actions/workflows/ci.yml)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SkyCoords.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/JuliaAstro/SkyCoords.jl/graph/badge.svg?token=0WIe7bWYFj)](https://codecov.io/gh/JuliaAstro/SkyCoords.jl)

SkyCoords.jl provides a type system for astronomical coordinate systems with appropriate conversions between them.

## Installation
From the Julia REPL

```julia-repl
pkg> add SkyCoords

julia> using SkyCoords
```

```@meta
DocTestSetup = :(using SkyCoords)
```

## Usage


There are currently five supported coordinate systems. The following
immutable types are used to represent coordinates in each system:

- [`ICRSCoords`](@ref): ICRS coordinates system
- [`GalCoords`](@ref): Galactic coordinates system
- [`SuperGalCoords`](@ref): Supergalactic coordinates system
- [`FK5Coords`](@ref): FK5 coordinates system (with arbitrary equinox)
- [`EclipticCoords`](@ref): Ecliptic coordinates system

Each type holds a longitude and latitude, and each is a subtype of
[`AbstractSkyCoords`](@ref).

```jldoctest
julia> c1 = ICRSCoords(0.0, 0.0)  # inputs are ra, dec in radians
ICRSCoords{Float64}(0.0, 0.0)

julia> c1.ra # access ra, dec individually
0.0

julia> c2 = convert(GalCoords, c1) # convert to a different system
GalCoords{Float64}(1.6814027872278692, -1.0504884034813007)

julia> c2.l # Note that galactic coordinate fields are l, b
1.6814027872278692

julia> c1 |> FK5Coords{2000} # Can use piping syntax for conversion
FK5Coords{2000, Float64}(1.1102233723050067e-7, 4.411803426976326e-8)
```
### Units

There is built-in support for units via [Unitful.jl](https://github.com/PainterQubits/Unitful.jl)

```jldoctest unitangles
julia> using Unitful

julia> c = ICRSCoords(0.11255u"°", 0.00091u"rad")
ICRSCoords{Float64}(0.0019643680731196178, 0.00091)

julia> c2 = FK5Coords{2000}(0.1u"rad", 0.5)
FK5Coords{2000, Float64}(0.1, 0.5)

julia> SkyCoords.lat(u"μrad", c)
910.0 μrad
```

### Parsing from strings

The [AstroAngles.jl](https://github.com/JuliaAstro/AstroAngles.jl) package provides convenient string parsing utilities

```jldoctest astroangles
julia> using AstroAngles

julia> c3 = ICRSCoords(hms"05:34:31.94", dms"+22:00:52.2")
ICRSCoords{Float64}(1.4596726677614607, 0.3842255081802917)
```

for example, to load coordinates from a target list

```julia-repl
julia> using CSV, DataFrames

julia> table = CSV.File("target_list.csv") |> DataFrame;

julia> [table.ra table.dec]
203×2 Matrix{String}:
 "00 05 01.42"  "40 03 35.82"
 "00 05 07.52"  "73 13 11.34"
 "00 36 01.40"  "-11 12 13.00"
[...]

julia> coords = @. ICRSCoords(hms2rad(table.ra), dms2rad(table.dec))
203-element Vector{ICRSCoords{Float64}}:
 ICRSCoords{Float64}(0.021919880964005448, 0.6991780256843024)
 ICRSCoords{Float64}(0.022363485482220672, 1.277926878539953)
 ICRSCoords{Float64}(0.15718144355252264, -0.19553990200190915)
[...]
```

### Angular Separation between Coordinates

The [`separation`](@ref) function allows you to compute the angular (great-circle)
distance between two coordinates, in radians, using
the [Vincenty formula](http://en.wikipedia.org/wiki/Great-circle_distance).  The
coordinates can be also given in different systems.  For example, according to
SIMBAD the `FK5Coords{2000}` coordinates
of [Mizar](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=MIZAR) are

```jldoctest sep
julia> mizar = FK5Coords{2000}(3.507787, 0.958628)
FK5Coords{2000, Float64}(3.507787, 0.958628)
```

while the `GalCoords` coordinates of Alcor are

```jldoctest sep
julia> alcor = GalCoords(1.968189, 1.072829)
GalCoords{Float64}(1.968189, 1.072829)
```

Their angular separation is given by

```jldoctest sep
julia> separation(mizar, alcor) # Radians
0.003435309169452965

julia> rad2deg(separation(mizar, alcor)) * 60 # Arcminutes
11.809723003934822
```

with an angle

```jldoctest sep
julia> position_angle(mizar, alcor) # radians
1.2446024012417884

julia> position_angle(mizar, alcor) |> rad2deg # degrees
71.31046476300233

```

## Catalog Matching

SkyCoords.jl offers coordinate catalog matching functionality through an extension that depends on [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl). This functionality requires Julia ≥ v1.9 and `NearestNeighbors.jl` to be loaded (e.g., `using NearestNeighbors.jl`).

The [`SkyCoords.match`](@ref) function can match two catalogs of coordinates with an interface similar to Astropy's `match_coordinates_sky`. This function operates on two arrays of coordinates, the first being the "reference" catalog that will be searched to find the closest coordinates to those in the second catalog. This function returns the indices into the reference catalog of the matches and the angular separation (in radians) between each coordinate and its match in the reference catalog.

```jldoctest matching
using NearestNeighbors # Required to use `match` method
using SkyCoords
# Generate random coordinates
N = 1000
lons = 2pi .* rand(N) # (0, 2π)
lats = pi .* (rand(N) .- 0.5) # (-π, π)
# The catalog to match against
refcat = ICRSCoords.(lons, lats)
# The catalog of coordinates for which you want to find neighbors in "refcat"
matchcat = refcat[[1,5,10]]

ids, sep = SkyCoords.match(refcat, matchcat)
ids == [1,5,10] # Indices for which `refcat[ids]` match to `matchcat`
# output
true
```

Note that [`SkyCoords.match`](@ref) is not exported (to avoid clashing with `Base.match`) and should be used via the qualified signature `SkyCoords.match` (as above) or explicitly imported (e.g., `using SkyCoords: match`).

This extension additionally supports construction of [`NearestNeighbors.KDTree`](@ref KDTree)s from `AbstractArray{<:AbstractSkyCoords}` and extends methods for general nearest neighbors queries ([`nn`](@ref), [`knn`](@ref)) and queries for all neighbors within a given separation ([`inrange`](@ref), similar to Astropy's `search_around_sky`).

More complicated catalog joins are supported by the [FlexiJoins.jl](https://github.com/JuliaAPlavin/FlexiJoins.jl) package. For example, if `L` and `R` are two catalogs with coordinate keys `:coordsL` and `:coordsR` respectively, the two catalogs can be joined based on angular separation with `FlexiJoins.innerjoin((L, R), FlexiJoins.by_distance(:coordsL, :coordsR, SkyCoords.separation, <=(0.1)))` where the final condition indicates you only want to keep matches that have separations less than or equal to 0.1 rad. See their documentation on astronomy-specific applications [here](https://aplavin.github.io/FlexiJoins.jl/notebooks/skycoords.html).

## Accuracy

All the supported conversions have been compared to the results of
[astropy.coordinates](https://docs.astropy.org/en/stable/coordinates/) (to better than 0.0001 arcsec agreement for `Float64`).
In turn, [astropy.coordinates](https://docs.astropy.org/en/stable/coordinates/) has been tested against many other tools.

## Performance

For small and moderate numbers of coordinates, conversions are much
faster than [astropy.coordinates](https://docs.astropy.org/en/stable/coordinates/) in Python. The following plot shows the
performance for converting ICRS coordinates to various other systems
(Galactic, FK5J2000 and FK5J1975), using [astropy.coordinates](https://docs.astropy.org/en/stable/coordinates/) (`py_*`
labels) and SkyCoords.jl (`jl_*` labels). The x axis denotes the
number of coordinates being simultaneously converted, with 1
corresponding to scalar coordinates.

![times](assets/bench.png)

| Specs           |                                        |
|:----------------|:---------------------------------------|
| CPU             | Intel core i5-8259U @ 2.3GHz (4 cores) |
| RAM             | 16GB                                   |
| Julia Version   | 1.2                                    |
| Python Version  | 3.7                                    |
| Astropy Version | 3.1.2                                  |

For scalar coordinates, SkyCoords.jl is up to 100,000 times
faster. For very large vectors of one million coordinates or more,
SkyCoords.jl is 2-4 times faster.  The source code for these
benchmarks can be found in `bench/`.

## Contributing

If you would like to contribute to SkyCoords please head over to the [GitHub page](https://github.com/juliaastro/skycoords.jl) and file an issue or open a pull request!
