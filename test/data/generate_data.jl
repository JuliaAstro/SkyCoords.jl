# Script to generate reference coordinate transformations for tests
# using astropy.coordinates
using DelimitedFiles
using PythonCall

## python imports
apc = pyimport("astropy.coordinates")
apu = pyimport("astropy.units")

## data generation
N = 1000
lons = 2pi .* rand(rng, N) # (0, 2π)
lats = pi .* (rand(rng, N) .- 0.5) # (-π, π)

open(coords_path, "w") do fh
    writedlm(fh, [lons lats], ',')
end
@info "Input coordinates saved to file" filename=abspath(coords_path)

coordinate_systems = Dict(
    "icrs" => apc.ICRS,
    "fk5j2000" => apc.FK5(equinox="J2000.0"),
    "fk5j1975" => apc.FK5(equinox="J1975.0"),
    "gal" => apc.Galactic,
    "eclip" => apc.GeocentricMeanEcliptic
)

## run conversions
for (input_coord_key, input_coord) in coordinate_systems
    # generate inputs
    input_coord_list = apc.SkyCoord(lons, lats, unit=("rad", "rad"), frame=input_coord)
    for (output_coord_key, output_coord) in coordinate_systems
        input_coord_key == output_coord_key && continue # skip same systems

        # convert
        output_coord_list = input_coord_list.transform_to(output_coord)

        # have to copy data from python
        output_lons = pyconvert(Vector, output_coord_list.spherical.lon.rad)
        output_lats = pyconvert(Vector, output_coord_list.spherical.lat.rad)

        # save to file
        output_path = joinpath(datapath, "$(input_coord_key)_to_$(output_coord_key).csv")
        open(output_path, "w") do fh
            writedlm(fh, [output_lons output_lats], ',')
        end
        @info "Generated coords from $input_coord_key to $output_coord_key" filename=abspath(output_path)
    end
end
