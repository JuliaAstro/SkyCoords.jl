include("SkyCoords.jl")
using SkyCoords
using Base.Test

data, hdr = readcsv("../testdata/icrs_fk5.csv", has_header=true)

#for i=1:3
#    c = ICRSCoords(data[i,1], data[i,2])
#    g = to_gal(c)
#    @printf "astropy: %f,%f  julia: %f,%f\n" data[i,3] data[i,4] g.l g.b
#end
