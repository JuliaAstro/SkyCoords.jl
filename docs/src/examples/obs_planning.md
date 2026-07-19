# Observation planning with horizontal (alt/az) coordinates

```@example obs
using SkyCoords, SOFA
```

First, we create the instance of the object in the sky we want to observe:

```@example obs
m13 = ICRSCoords(deg2rad(250.423475), deg2rad(36.4613194))
```

Then we define the observation location, for example at the top of Mount Wilson:

```@example obs
mt_wilson = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742)
```
!!! note
    [`Observer`](@ref) takes the geodetic latitude and east-positive longitude in radians, and the altitude above the WGS84 ellipsoid in meters.

We will observe at 4:00 UTC on 2021 Nov 08, which corresponds to 8:00 PM the night before local time. A time package such as [AstroTime.jl](@extref AstroTime :doc:`index`) can be useful for computing the required UTC Julian date from a calendar date, but here we write it directly:

```@example obs
jd = 2459526.5 + 4 / 24
```

Now we perform the conversion to alt/az coordinates:

```@example obs
altaz = AltAzCoords(m13, mt_wilson, jd)
```

We may want to see this result in degrees, which is easy with [`Base.Math.rad2deg`](@extref). M13 is visible at an elevation of 13.4 degrees off the horizon and at a (true north) compass heading of 305.2 degrees, which is WNW.

```@example obs
rad2deg(altaz.alt), rad2deg(altaz.az)
```

For higher precision, we can supply the Earth orientation parameters published in the IERS bulletins (UT1-UTC and polar motion), and for the apparent, refracted position we can supply the ambient weather conditions:

```@example obs
refracted = AltAzCoords(
    m13, mt_wilson, jd;
    dut1 = -0.1104, pressure = 820, temperature = 10, relative_humidity = 0.4,
)
```

```@example obs
rad2deg(altaz.alt), rad2deg(altaz.az)
```

We can also plan across the night by broadcasting over a range of times, e.g., to find when the target is above 13 degrees:

```@example obs
jds = jd .+ range(-6 / 24, 6 / 24, length = 100)

alts = [AltAzCoords(m13, mt_wilson, t).alt for t in jds]

observable = jds[alts .> deg2rad(13)]

println("M13 is above 13° between JD $(first(observable)) and JD $(last(observable))")
```

And go the other way, from a local direction back to the celestial sphere:

```@example obs
zenith = AltAzCoords(π / 2, 0)

ICRSCoords(zenith, mt_wilson, jd)
```
