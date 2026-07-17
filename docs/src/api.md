# API/Reference

```@meta
DocTestSetup = :(using SkyCoords)
```

## Index

```@index
Pages = ["api.md"]
```

## Types

```@docs
AbstractSkyCoords
ICRSCoords
GalCoords
SuperGalCoords
FK4Coords
FK4NoETerms
FK5Coords
EclipticCoords
AltAzCoords
CartesianCoords
ProjectedCoords
Observer
```

## Conversion

To convert between types, there are three (equivalent) methods of doing so.

```jldoctest convsetup
julia> c1 = ICRSCoords(0., 0.)
ICRSCoords{Float64}(0.0, 0.0)
```

- using `convert`
  ```jldoctest convsetup
  julia> convert(GalCoords, c1)
  GalCoords{Float64}(1.6814027872278692, -1.0504884034813007)
  ```
- using constructors
  ```jldoctest convsetup
  julia> GalCoords(c1)
  GalCoords{Float64}(1.6814027872278692, -1.0504884034813007)
  ```
- using `|>`
  ```jldoctest convsetup
  julia> c1 |> GalCoords
  GalCoords{Float64}(1.6814027872278692, -1.0504884034813007)
  ```

### Horizontal (alt/az) coordinates

Conversions to and from [`AltAzCoords`](@ref) depend on where and when the observation takes place, so they are not available through `convert`. Instead, load [Astrometry.jl](https://github.com/JuliaAstro/Astrometry.jl) (a weak dependency providing the IAU SOFA algorithms) and pass an [`Observer`](@ref) along with a UTC Julian date:

```julia
julia> using Astrometry

julia> mt_wilson = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742);

julia> jd = 2459526.5 + 4 / 24;  # 2021-11-08T04:00 UTC

julia> m13 = ICRSCoords(deg2rad(250.423475), deg2rad(36.4613194));

julia> altaz = AltAzCoords(m13, mt_wilson, jd)
AltAzCoords{Float64}(0.23313661745182013, 5.3268602406305074)

julia> ICRSCoords(altaz, mt_wilson, jd) ≈ m13
true
```

See the [`AltAzCoords`](@ref) documentation for the supported keyword arguments (Earth orientation parameters and atmospheric refraction).

## Catalog Matching

```@docs
SkyCoords.match
KDTree
nn
knn
inrange
```

## Functions

```@docs
separation
position_angle
offset
cartesian
spherical
```
