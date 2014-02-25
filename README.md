# SkyCoords

Astronomical coordinate systems in Julia

## Install

```julia
julia> Pkg.add("git://github.com/kbarbary/SkyCoords.jl.git")
```

## Usage

```julia
julia> using SkyCoords

julia> c = ICRSCoords(0., 0.)
ICRSCoords(0.0,0.0)

julia> c.ra
0.0

julia> to_fk5j2000(c)
FK5J2000Coords(1.1102233723050067e-7,4.411803426976326e-8)
```
