using SkyCoords, AstroTime, AstroAngles, AstroLib

# First we create the instance of the object in the sky we want to observe
M13 = ICRSCoords(250.423475 |> deg2rad, 36.4613194 |> deg2rad)

# Then we define the observation time and location
# We will use the top of Mount Wilson at 4AM UTC at Nov 8, 2021 which corresponds
# to 8PM the night before local time
mt_wilson = Observatory("Mount Wilson",32.2264,-118.0642,1693.25,-8)
time = from_utc("2021-11-08T04:00")
# We will need the julian date for the conversions, which we can do now
jd = julian(time) |> value

# Now we perform the conversion
altaz = AltAzCoords(M13,jd,mt_wilson)

# We may want to see this result in degrees, which is easy with rad2deg
# For this result, M13 is visible at an elevation of 12 degrees off the horizon
# and at a (true north) compass heading of 305.8, which is WNW.

# FIXME plot over time
# This is a few tenths of a degree off the astropy solution - why the difference?
