# SkyCoords

Astronomical coordinate systems in Julia

## Install

```julia
julia> Pkg.add("git://github.com/kbarbary/SkyCoords.jl.git")
```

## Usage

#### Static systems: `ICRS`, `Galactic`, `FK5J2000`

```julia
julia> using SkyCoords

julia> c1 = ICRS(0., 0.)  # ra, dec in radians
ICRS(0.0,0.0)

julia> c1.ra
0.0

julia> c2 = to_galactic(c1)
Galactic(1.681404315278054,-1.0504869904089078)

julia> c2.l  # galactic coordinate fields are l, b rather than ra, dec
1.681404315278054

# FK5J2000 is a version of FK5 fixed at epoch J2000
julia> to_fk5j2000(c1)
FK5J2000(1.1102233710147402e-7,4.411803426976326e-8)
```


#### Dynamic systems: `FK5`

These types include an extra field, `equinox`. The equinox must be specified
when creating the object or converting from another coordinate system:

```julia
julia> c2 = to_fk5(c1, 1970.)
FK5(6.276477890771361,-0.0029150998550493794,1970.0)

julia> to_icrs(c2)
ICRS(0.0,-8.673617379884037e-19)  # pretty close to round-tripping
```


## Speed

For small numbers of points, this package can be much faster than
`astropy.coordinates` in Python. The following is an example of
transforming points from ICRS to Galactic. The speed difference is
factor of approximately 10,000 for 10 coordinates, 80 for 1000
coordinates and 4 for 100,000 coordinates.

**astropy.coordinates**

```python
In [1]: from numpy import pi

In [2]: from numpy.random import rand

In [3]: from astropy.coordinates import SkyCoord

In [4]: for n in [10, 1000, 100000]:
   ...:     c = SkyCoord(2.*pi*rand(n), pi*rand(n)-pi/2, unit=('rad', 'rad'))
   ...:     %timeit c.galactic
   ...: 
100 loops, best of 3: 12.5 ms per loop
100 loops, best of 3: 13 ms per loop
10 loops, best of 3: 66.7 ms per loop
```

**SkyCoords.jl**

```jlcon
julia> using SkyCoords

julia> using TimeIt

julia> for n in [10, 1000, 100000]
           c = [ICRS(2pi*rand(), pi*(rand() - 0.5)) for i=1:n]
           @timeit to_galactic(c)
       end
100000 loops, best of 3: 1.33 µs per loop
1000 loops, best of 3: 163.94 µs per loop
10 loops, best of 3: 16.23 ms per loop
```
