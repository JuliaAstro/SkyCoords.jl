using PythonCall

## python imports
apc = pyimport("astropy.coordinates")

## data generation
N = 1000
const lons = 2pi .* rand(rng, N) # (0, 2π)
const lats = pi .* (rand(rng, N) .- 0.5) # (-π, π)

# get the astropy frames from the julia type
astropy_conversion(::Type{<:ICRSCoords}) = apc.ICRS
astropy_conversion(::Type{<:FK5Coords{F}}) where {F} = apc.FK5(equinox="J$F")
astropy_conversion(::Type{<:GalCoords}) = apc.Galactic
astropy_conversion(::Type{<:EclipticCoords{F}}) where {F} = apc.GeocentricMeanEcliptic(equinox="J$F")

function test_against_astropy(intype, outtype; atol=0)
    ## get julia values
    output_coords = map((lon, lat) -> convert(outtype, intype(lon, lat)), lons, lats)
    output_lons = map(lon, output_coords)
    output_lats = map(lat, output_coords)

    ## get astropy values
    ap_input_coord = astropy_conversion(intype)
    ap_output_coord = astropy_conversion(outtype)
    input_coord_list = apc.SkyCoord(lons, lats, unit=("rad", "rad"), frame=ap_input_coord)
    output_coord_list = input_coord_list.transform_to(ap_output_coord)
    # have to copy data from python
    ap_output_lons = pyconvert(Vector, output_coord_list.spherical.lon.rad)
    ap_output_lats = pyconvert(Vector, output_coord_list.spherical.lat.rad)

    @test output_lons ≈ ap_output_lons  atol=atol*√N
    @test output_lats ≈ ap_output_lats  atol=atol*√N
    
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
        FK5Coords{2000,F},
        FK5Coords{1975,F},
        GalCoords{F},
        EclipticCoords{2000,F},
        EclipticCoords{1975,F},
    )
    @testset "$IN_SYS --> $OUT_SYS" for IN_SYS in systems, OUT_SYS in systems
        atol = IN_SYS <: EclipticCoords || OUT_SYS <: EclipticCoords ? 1e-3 : 0
        test_against_astropy(IN_SYS, OUT_SYS; atol)
    end
end
