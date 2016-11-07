#!/usr/bin/env julia

# Tests against astropy.

using SkyCoords
using Base.Test

datapath = joinpath(dirname(@__FILE__), "data")

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

lon(c::GalCoords) = c.l
lat(c::GalCoords) = c.b
lon(c::AbstractSkyCoords) = c.ra
lat(c::AbstractSkyCoords) = c.dec

function angsep{T<:AbstractSkyCoords}(c1::T, c2::T)
    angsep(lon(c1), lat(c1), lon(c2), lat(c2))
end

function angsep{T<:AbstractSkyCoords}(c1::Array{T}, c2::Array{T})
    size(c1) == size(c2) || error("size mismatch")
    result = similar(c1, SkyCoords._eltype(first(c1)))
    for i=1:length(c1)
        result[i] = angsep(c1[i], c2[i])
    end
    result
end

rad2arcsec(r) = 3600.*rad2deg(r)

# input coordinates
fname = joinpath(datapath, "input_coords.csv")
indata, inhdr = readcsv(fname; header=true)

# Float32 has a large tolerance compared to Float64 and BigFloat, but here we
# are more interested in making sure that the infrastructure works for different
# floating types.
for (F, TOL) in ((Float32, 0.2), (Float64, 0.0001), (BigFloat, 0.0001))
    println("Testing type ", F)

    for (insys, T) in (("icrs", ICRSCoords{F}), ("fk5j2000", FK5Coords{2000,F}),
                       ("fk5j1975", FK5Coords{1975,F}), ("gal", GalCoords{F}))

        c_in = T[T(indata[i, 1], indata[i, 2]) for i=1:size(indata,1)]

        for (outsys, S) in (("icrs", ICRSCoords{F}), ("fk5j2000", FK5Coords{2000,F}),
                            ("fk5j1975", FK5Coords{1975,F}), ("gal", GalCoords{F}))
            (outsys == insys) && continue
            c_out = S[convert(S, c) for c in c_in]

            # Read in reference answers.
            fname = joinpath(datapath, "$(insys)_to_$(outsys).csv")
            refdata, hdr = readcsv(fname; header=true)
            c_ref = S[S(refdata[i, 1], refdata[i, 2]) for i=1:size(refdata,1)]

            # compare
            sep = angsep(c_out, c_ref)
            maxsep = rad2arcsec(maximum(sep))
            meansep = rad2arcsec(mean(sep))
            minsep = rad2arcsec(minimum(sep))
            @printf "%8s --> %8s : max=%6.4f\"  mean=%6.4f\"   min=%6.4f\"\n" insys outsys maxsep meansep minsep
            @test maxsep < TOL
        end
    end
end

# Test conversion with mixed floating types.
c1 = ICRSCoords(e, pi/2)
for T in (GalCoords, FK5Coords{2000})
    c2 = convert(T{Float32}, c1)
    c3 = convert(T{Float64}, c1)
    c4 = convert(T{BigFloat}, c1)
    @test typeof(c2) === T{Float32}
    @test typeof(c3) === T{Float64}
    @test typeof(c4) === T{BigFloat}
    @test isapprox(lat(c2), lat(c3), rtol=sqrt(eps(Float32)))
    @test isapprox(lat(c3), lat(c4), rtol=sqrt(eps(Float64)))
    @test isapprox(lon(c2), lon(c3), rtol=sqrt(eps(Float32)))
    @test isapprox(lon(c3), lon(c4), rtol=sqrt(eps(Float64)))
end

println()
println("All tests passed.")
