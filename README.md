SkyCoords.jl
============

Astronomical coordinate systems in Julia

## Install

```julia
julia> Pkg.add("git://github.com/kbarbary/SkyCoords.jl.git")
```

## Usage

There are three supported coordinate systems:

- `ICRSCoords`
- `GalacticCoords`
- `FK5Coords`

```julia
julia> using SkyCoords

# create a coordinates object
julia> c1 = ICRSCoords(0., 0.)  # inputs are ra, dec in radians
ICRSCoords(0.0,0.0)

# access ra, dec individually
julia> c1.ra
0.0

# convert to a different system
julia> c2 = convert(GalacticCoords, c1)
GalacticCoords(1.681404315278054,-1.0504869904089078)

# Note that galactic coordinate fields are l, b
julia> c2.l
1.681404315278054

# Note that FK5Coords is parameterized on equinox
julia> convert(FK5Coords{2000}, c1)
FK5Coords{2000}(1.1102233710147402e-7,4.411803426976326e-8)

# Arrays of coordinates

# create an array of coordinates 
julia> c1 = [ICRSCoords(0., 0.) for i=1:3]
3-element Array{ICRSCoords,1}:
 ICRSCoords(0.0,0.0)
 ICRSCoords(0.0,0.0)
 ICRSCoords(0.0,0.0)

# convert entire array to a different system
julia> convert(Vector{GalacticCoords}, c1)
3-element Array{GalacticCoords,1}:
 GalacticCoords(1.681404315278054,-1.0504869904089078)
 GalacticCoords(1.681404315278054,-1.0504869904089078)
 GalacticCoords(1.681404315278054,-1.0504869904089078)

# There's no performance gain from using this "vectorized" convert,
# except conversions to/from FK5Coords, where the equinox precession
# can be done just once for the entire vector, leading to a modest ~2x
# speed up.
julia> convert(Vector{FK5Coords{1975}}, c1)
3-element Array{FK5Coords{1975},1}:
 FK5Coords{1975}(6.277595732508468,-0.0024292220493946897)
 FK5Coords{1975}(6.277595732508468,-0.0024292220493946897)
 FK5Coords{1975}(6.277595732508468,-0.0024292220493946897)
```


## Performance

For small numbers of coordinates, conversions are *much* faster than
astropy.coordinates in Python. The follow plot shows the performance
for converting ICRS coordinates to various other systems (Galactic,
FK5J2000 and FK5J1975), using astropy.coordinates (`py_*`) and
SkyCoords.jl (`jl_*`). The x axis denotes the number of coordinates
being simultaneously converted, with 1 cooresponding to scalar
coordinates.

![times](bench/bench.png)

For scalar coordinates, SkyCoords.jl is up to *100,000 times*
faster. Even for a vector of one million coordinates, SkyCoords.jl is
still 2-4 times faster.  The source code for these benchmarks can be
found in `bench/`.
