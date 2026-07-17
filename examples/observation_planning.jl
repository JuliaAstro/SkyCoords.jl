# Observation planning with horizontal (alt/az) coordinates.
#
# Conversions to and from `AltAzCoords` are provided by the Astrometry.jl
# package extension, so make sure Astrometry.jl is loaded.
using SkyCoords
using Astrometry

# First we create the instance of the object in the sky we want to observe
m13 = ICRSCoords(deg2rad(250.423475), deg2rad(36.4613194))

# Then we define the observation location: the top of Mount Wilson.
# `Observer` takes the geodetic latitude and east-positive longitude in
# radians, and the altitude above the WGS84 ellipsoid in meters.
mt_wilson = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742)

# We will observe at 4AM UTC on Nov 8, 2021, which corresponds to 8PM the
# night before local time. The conversion needs the UTC Julian date; a time
# package such as AstroTime.jl can compute it from a calendar date, here we
# write it directly: JD 2459526.5 is 2021-11-08T00:00 UTC.
jd = 2459526.5 + 4 / 24

# Now we perform the conversion
altaz = AltAzCoords(m13, mt_wilson, jd)

# We may want to see this result in degrees, which is easy with rad2deg.
# M13 is visible at an elevation of 13.4 degrees off the horizon and at a
# (true north) compass heading of 305.2 degrees, which is WNW.
@show rad2deg(altaz.alt) rad2deg(altaz.az)

# For higher precision we can supply the Earth orientation parameters
# published in the IERS bulletins (UT1-UTC and polar motion), and for the
# apparent, refracted position we can supply the ambient weather conditions:
refracted = AltAzCoords(
    m13, mt_wilson, jd;
    dut1 = -0.1104, pressure = 820, temperature = 10, relative_humidity = 0.4,
)
@show rad2deg(refracted.alt)

# We can also plan across the night by broadcasting over a range of times,
# e.g. to find when the target is above 13 degrees:
jds = jd .+ range(-6 / 24, 6 / 24, length = 100)
alts = [AltAzCoords(m13, mt_wilson, t).alt for t in jds]
observable = jds[alts .> deg2rad(13)]
println("M13 is above 13° between JD $(first(observable)) and JD $(last(observable))")

# And go the other way, from a local direction back to the celestial sphere:
zenith = AltAzCoords(π / 2, 0)
@show ICRSCoords(zenith, mt_wilson, jd)
