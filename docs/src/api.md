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
FK5Coords
EclipticCoords
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

## Catalog Matching

```@docs
SkyCoords.match(::KDTree, ::AbstractArray{<:AbstractSkyCoords})
KDTree(::AbstractArray{<:SkyCoords.AbstractSkyCoords})
nn(::KDTree, ::SkyCoords.AbstractSkyCoords)
knn(::KDTree, coord::SkyCoords.AbstractSkyCoords, ::Int, ::Bool)
inrange(::KDTree, ::SkyCoords.AbstractSkyCoords, ::Number)
```

## Functions

```@docs
separation
position_angle
offset
```
