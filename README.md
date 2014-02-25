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

For small numbers of points, this can be much faster than
`astropy.coordinates`. Here is an example transforming 100 points
from ICRS to Galactic:

**astropy.coordinates**

```python
In [1]: from astropy.coordinates import ICRS

In [2]: import numpy as np

In [3]: ra = 2.* np.pi * np.random.rand(100)

In [4]: dec = np.pi * np.random.rand(100) - np.pi/2.

In [5]: c1 = ICRS(ra, dec, unit=('rad', 'rad'))

In [6]: %timeit c1.galactic
100 loops, best of 3: 6.68 ms per loop
```

**SkyCoords**

```julia
julia> using SkyCoords

julia> using TimeIt

julia> ra = 2pi*rand(100)

julia> dec = pi*rand(100) - pi/2.

julia> c1 = [ICRS(ra[i], dec[i]) for i=1:100]

julia> @timeit to_galactic(c1)
10000 loops, best of 3: 61.16 µs per loop
```