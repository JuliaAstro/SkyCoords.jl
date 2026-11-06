using PythonCall

## python imports
apc = pyimport("astropy.coordinates")

## data generation
N = 1000
const lons = 2pi .* rand(rng, N) # (0, 2π)
const lats = pi .* (rand(rng, N) .- 0.5) # (-π, π)

# get the astropy frames from the julia type
astropy_conversion(::Type{<:ICRSCoords}) = apc.ICRS
astropy_conversion(::Type{<:FK4Coords{F}}) where {F} = apc.FK4(equinox = "B$F")
astropy_conversion(::Type{<:FK4NoETerms{F}}) where {F} = apc.FK4NoETerms(equinox = "B$F")
astropy_conversion(::Type{<:FK5Coords{F}}) where {F} = apc.FK5(equinox = "J$F")
astropy_conversion(::Type{<:GalCoords}) = apc.Galactic
astropy_conversion(::Type{<:SuperGalCoords}) = apc.Supergalactic
astropy_conversion(::Type{<:EclipticCoords{F}}) where {F} = apc.GeocentricMeanEcliptic(equinox = "J$F")

function test_against_astropy(intype, outtype; atol = 0)
    ## get julia values
    output_coords = map((lon, lat) -> convert(outtype, intype(lon, lat)), lons, lats)
    output_lons = map(lon, output_coords)
    output_lats = map(lat, output_coords)

    ## get astropy values
    ap_input_coord = astropy_conversion(intype)
    ap_output_coord = astropy_conversion(outtype)
    input_coord_list = apc.SkyCoord(lons, lats, unit = ("rad", "rad"), frame = ap_input_coord)
    output_coord_list = input_coord_list.transform_to(ap_output_coord)
    # have to copy data from python
    ap_output_lons = pyconvert(Vector, output_coord_list.spherical.lon.rad)
    ap_output_lats = pyconvert(Vector, output_coord_list.spherical.lat.rad)

    @test output_lons ≈ ap_output_lons  atol = atol * √N
    return @test output_lats ≈ ap_output_lats  atol = atol * √N

end

# Float32 has a large tolerance compared to Float64 and BigFloat, but here we
# are more interested in making sure that the infrastructure works for different
# floating types.
@testset "Testing $F" for (F, TOL) in (
        (Float32, 0.2),
        (Float64, 0.0001),
        (BigFloat, 0.0001),
    )

    systems = (
        ICRSCoords{F},
        FK4Coords{1950, F},
        FK4Coords{1975, F},
        FK4NoETerms{1950, F},
        FK4NoETerms{1975, F},
        FK5Coords{2000, F},
        FK5Coords{1975, F},
        GalCoords{F},
        SuperGalCoords{F},
        EclipticCoords{2000, F},
        EclipticCoords{1975, F},
    )
    @testset "$IN_SYS --> $OUT_SYS" for IN_SYS in systems, OUT_SYS in systems
        # FK4Coords/FK4NoETerms need a small amount of slack against (Super)Galactic:
        # astropy routes FK4->Galactic through a *different* (legacy, B1950-based)
        # galactic pole than the J2000-based one used everywhere else in this package
        # (see comment on FK5J2000_TO_GAL), so precessing an off-equinox FK4/FK4NoETerms
        # to Galactic picks up a small (sub-arcsecond) discrepancy between the two
        # conventions. This is independent of the E-terms, so it affects both types.
        atol = if IN_SYS <: EclipticCoords || OUT_SYS <: EclipticCoords
            1.0e-3
        elseif IN_SYS <: FK4Coords || OUT_SYS <: FK4Coords || IN_SYS <: FK4NoETerms || OUT_SYS <: FK4NoETerms
            1.0e-5
        else
            0
        end
        test_against_astropy(IN_SYS, OUT_SYS; atol)
    end
end

@testset "AltAz" begin
    apt = pyimport("astropy.time")
    iers = pyimport("astropy.utils.iers")
    iers.conf.auto_download = false

    # Mount Wilson at 2021-11-08T04:00 UTC
    observer = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742.0)
    location = apc.EarthLocation.from_geodetic(
        lon = rad2deg(observer.longitude), lat = rad2deg(observer.latitude),
        height = observer.altitude,
    )
    jd = 2459526.5 + 4 / 24
    time = apt.Time(jd, format = "jd", scale = "utc", location = location)

    # Look up the Earth orientation parameters astropy uses, so that both
    # libraries run with identical inputs.
    eop = iers.earth_orientation_table.get()
    dut1 = pyconvert(Float64, eop.ut1_utc(time).to_value("s"))
    pm = eop.pm_xy(time)
    xp = deg2rad(pyconvert(Float64, pm[0].to_value("deg")))
    yp = deg2rad(pyconvert(Float64, pm[1].to_value("deg")))

    frame = apc.AltAz(obstime = time, location = location, pressure = 0)
    n = 100
    icrs_list = apc.SkyCoord(lons[1:n], lats[1:n], unit = ("rad", "rad"), frame = apc.ICRS)
    ap_altaz = icrs_list.transform_to(frame)
    ap_alts = pyconvert(Vector, ap_altaz.alt.rad)
    ap_azs = pyconvert(Vector, ap_altaz.az.rad)

    atol = 1.0e-8 # radians, ~2 μas
    skycoords_frame = AltAzFrame(observer, jd; dut1, xp, yp)
    for (lon, lat, ap_alt, ap_az) in zip(lons[1:n], lats[1:n], ap_alts, ap_azs)
        altaz = AltAzCoords(ICRSCoords(lon, lat), skycoords_frame)
        @test separation(altaz, AltAzCoords(ap_alt, ap_az)) < atol
        icrs = ICRSCoords(AltAzCoords(ap_alt, ap_az), skycoords_frame)
        @test separation(icrs, ICRSCoords(lon, lat)) < atol
    end
end
