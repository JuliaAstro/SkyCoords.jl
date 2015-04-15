"""Script to generate reference coordinate transformations for tests, using
astropy.coordinates."""

from numpy import pi
from numpy.random import rand
from astropy.coordinates import SkyCoord, FK5, ICRS, Galactic
import astropy.units as u

# generate 1000 random coordinates and dates, save to file
n = 1000
lon = 2. * pi * rand(n)
lat = pi * rand(n) - pi / 2

inputfile = open("input_coords.csv", "w")
inputfile.write("lon,lat\n")
for i in range(n):
    inputfile.write("{0!r},{1!r}\n".format(lon[i], lat[i]))
inputfile.close()

# run conversions
coords = {'icrs': ICRS, 'fk5j2000': FK5, 'fk5j1975': FK5, 'gal': Galactic}
equinox = {'fk5j2000': 'J2000.0', 'fk5j1975': 'J1975.0'}
for insys in coords.keys():
    if insys in equinox:
        in_coords = coords[insys](lon * u.rad, lat * u.rad,
                                  equinox=equinox[insys])
    else:
        in_coords = coords[insys](lon * u.rad, lat * u.rad)

    for outsys in coords.keys():
        if outsys == insys:
            continue
        if outsys in equinox:
            out_coords = in_coords.transform_to(
                coords[outsys](equinox=equinox[outsys]))
        else:
            out_coords = in_coords.transform_to(coords[outsys])

        outlon = out_coords.spherical.lon.rad
        outlat = out_coords.spherical.lat.rad

        # write to file
        outputfile = open("{0:s}_to_{1:s}.csv".format(insys, outsys), "w")
        outputfile.write("lon,lat\n")
        for i in range(n):
            outputfile.write("{0!r},{1!r}\n".format(outlon[i], outlat[i]))
        outputfile.close()
