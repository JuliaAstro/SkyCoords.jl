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

We will observe at 4:00 UTC on 2021 Nov 08, which corresponds to 8:00 PM the night before local time. A time package such as [AstroTime.jl](@extref AstroTime :doc:`index`) can be useful for computing the required UTC Julian date from a calendar date, but here we write it directly. Together, the location and time define the observing context, which we describe with an [`AltAzFrame`](@ref):

```@example obs
jd = 2459526.5 + 4 / 24

frame = AltAzFrame(mt_wilson, jd)
```

Now we perform the conversion to alt/az coordinates:

```@example obs
altaz = AltAzCoords(m13, frame)
```

We may want to see this result in degrees, which is easy with [`Base.Math.rad2deg`](@extref). M13 is visible at an elevation of 13.4 degrees off the horizon and at a (true north) compass heading of 305.2 degrees, which is WNW.

```@example obs
rad2deg(altaz.alt), rad2deg(altaz.az)
```

For higher precision, we can supply the Earth orientation parameters published in the IERS bulletins (UT1-UTC and polar motion), and for the apparent, refracted position we can supply the ambient weather conditions:

```@example obs
frame_refracted = AltAzFrame(
    mt_wilson, jd;
    dut1 = -0.1104, pressure = 820, temperature = 10, relative_humidity = 0.4,
)

refracted = AltAzCoords(m13, frame_refracted)
```

```@example obs
rad2deg(refracted.alt), rad2deg(refracted.az)
```

We can also plan across the night by broadcasting over a range of times, e.g., to find when the target is above 13 degrees:

```@example obs
jds = jd .+ range(-6 / 24, 6 / 24, length = 100)

frames = AltAzFrame.(mt_wilson, jds)

alts = [AltAzCoords(m13, f).alt for f in frames]

observable = jds[alts .> deg2rad(13)]

println("M13 is above 13° between JD $(first(observable)) and JD $(last(observable))")
```

!!! tip
    A frame computes and caches its astrometry context (precession-nutation, Earth rotation, refraction constants) the first time it is used, so converting many coordinates against the same frame is cheap. Broadcast a whole catalog with `AltAzCoords.(catalog, frame)`, or a grid of targets and times with `AltAzCoords.(catalog, permutedims(frames))`, and the expensive setup runs only once per frame.

And go the other way, from a local direction back to the celestial sphere:

```@example obs
zenith = AltAzCoords(π / 2, 0)

ICRSCoords(zenith, frame)
```
