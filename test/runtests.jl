include("SkyCoords.jl")
using SkyCoords
using Base.Test

rad2deg = radians2degrees
deg2rad = degrees2radians

# Angular separation between two points (angles in radians)
#
#   The angular separation is calculated using the Vincenty formula [1]_,
#   which is slighly more complex and computationally expensive than
#   some alternatives, but is stable at at all distances, including the
#   poles and antipodes.
#
#   [1] http://en.wikipedia.org/wiki/Great-circle_distance
function angsep(lon1, lat1, lon2, lat2)

    sdlon = sin(lon2 - lon1)
    cdlon = cos(lon2 - lon1)
    slat1 = sin(lat1)
    slat2 = sin(lat2)
    clat1 = cos(lat1)
    clat2 = cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denom = slat1 * slat2 + clat1 * clat2 * cdlon

    atan2(sqrt(num1*num1 + num2*num2), denom)
end

# Tolerance is 0.03 arcsec
tol = deg2rad(0.03 / 3600.)

# The testdata file was generated with the ref_icrs_fk5.py script in
# astropy.coordinates.test.accuracy. The reference values were computed
# using AST.
data, hdr = readcsv("testdata/icrs_fk5.csv", has_header=true)

for i=1:size(data, 1)
    equinox = float(data[i,1][2:end])
    ra_in = deg2rad(data[i,3])
    dec_in = deg2rad(data[i,4])

    # ICRS to FK5 (assume ra_in, dec_in are ICRS)
    c1 = ICRS(ra_in, dec_in)
    c2 = to_fk5(c1, equinox)
    @test angsep(deg2rad(data[i,5]), deg2rad(data[i,6]), c2.ra, c2.dec) < tol

    # FK5 to ICRS (assume ra_in and dec_in are FK5)
    c1 = FK5(ra_in, dec_in, equinox)
    c2 = to_icrs(c1)
    @test angsep(deg2rad(data[i,7]), deg2rad(data[i,8]), c2.ra, c2.dec) < tol

end
