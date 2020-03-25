# API/Reference

```@meta
DocTestSetup = :(using SkyCoords)
```

## Index

```@index
pages = ["api.md"]
```

## Types

```@docs
AbstractSkyCoords
ICRSCoords
GalCoords
FK5Coords
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

## Functions

```@docs
SkyCoords.str2rad
separation
position_angle
```
