#!/usr/bin/env julia

# Tests against astropy.

using SkyCoords

using Test, DelimitedFiles, Printf, Statistics, Unitful

import SkyCoords: lat, lon

datapath = joinpath(dirname(@__FILE__), "data")

rad2arcsec(r) = 3600 * rad2deg(r)

# input coordinates
fname = joinpath(datapath, "input_coords.csv")
indata, inhdr = readdlm(fname, ','; header = true)

# Float32 has a large tolerance compared to Float64 and BigFloat, but here we
# are more interested in making sure that the infrastructure works for different
# floating types.
@testset "Testing $F" for (F, TOL) in ((Float32, 0.2), (Float64, 0.0001), (BigFloat, 0.0001))

    for (insys, T) in (("icrs", ICRSCoords{F}), ("fk5j2000", FK5Coords{2000,F}),
                       ("fk5j1975", FK5Coords{1975,F}), ("gal", GalCoords{F}))

        c_in = T[T(indata[i, 1], indata[i, 2]) for i = 1:size(indata, 1)]
        for c in c_in
            @test convert(T, c) == c
        end

        for (outsys, S) in (("icrs", ICRSCoords{F}), ("fk5j2000", FK5Coords{2000,F}),
                            ("fk5j1975", FK5Coords{1975,F}), ("gal", GalCoords{F}))
            (outsys == insys) && continue
            c_out = S[convert(S, c) for c in c_in]

            # Read in reference answers.
            ref_fname = joinpath(datapath, "$(insys)_to_$(outsys).csv")
            refdata, hdr = readdlm(ref_fname, ','; header = true)
            c_ref = S[S(refdata[i, 1], refdata[i, 2]) for i = 1:size(refdata, 1)]

            # compare
            sep = separation.(c_out, c_ref)
            maxsep = rad2arcsec(maximum(sep))
            meansep = rad2arcsec(mean(sep))
            minsep = rad2arcsec(minimum(sep))
            @printf "%8s --> %8s : max=%6.4f\"  mean=%6.4f\"   min=%6.4f\"\n" insys outsys maxsep meansep minsep
            @test maxsep < TOL
        end
    end
end

# Test separation between coordinates and conversion with mixed floating types.
@testset "Separation" begin
    c1 = ICRSCoords(ℯ, pi / 2)
    c5 = ICRSCoords(ℯ, 1 + pi / 2)
    @test separation(c1, c5) ≈ separation(c5, c1) ≈
        separation(c1, convert(GalCoords, c5)) ≈
        separation(convert(FK5Coords{1980}, c5), c1) ≈ 1
    for T in (GalCoords, FK5Coords{2000})
        c2 = convert(T{Float32}, c1)
        c3 = convert(T{Float64}, c1)
        c4 = convert(T{BigFloat}, c1)
        @test typeof(c2) === T{Float32}
        @test typeof(c3) === T{Float64}
        @test typeof(c4) === T{BigFloat}
        @test isapprox(lat(c2), lat(c3), rtol = sqrt(eps(Float32)))
        @test isapprox(lat(c3), lat(c4), rtol = sqrt(eps(Float64)))
        @test isapprox(lon(c2), lon(c3), rtol = sqrt(eps(Float32)))
        @test isapprox(lon(c3), lon(c4), rtol = sqrt(eps(Float64)))
        c6 = convert(T, c5)
        @test separation(c3, c6) ≈ separation(c6, c3) ≈ 1
    end
end

@testset "Pre-conversion of inputs" begin
    # Constructor of FK5Coords{1950} with non-float arguments
    @test typeof(FK5Coords{1950}(1, big(2))) == FK5Coords{1950,BigFloat}

    # Unitful
    @test ICRSCoords(20u"°", 0u"rad") == ICRSCoords(0.3490658503988659, 0.0)
    @test GalCoords(20u"°", 0u"rad") == GalCoords(0.3490658503988659, 0.0)
    @test FK5Coords{2000}(20u"°", 0u"rad") == FK5Coords{2000}(0.3490658503988659, 0.0)
end

@testset "Coord String" begin
    @test coord"0h0m0s" == coord"0°0'0\"" == 0.0
    @test coord"12h0m0s" ≈ coord"180°0'0\"" ≈ π
    @test coord"-12h0m0s" ≈ coord"-180°0'0\"" ≈ -π
    @test coord"+12h0m0s" ≈ coord"+180°0'0\"" ≈ π
    @test coord"18h36m56.33635s" ≈ 4.873565508047745
    @test coord"38°47'01.2802\"" ≈ 0.6769030681498798
    @test coord"-  10° 02'  10.885 \"" ≈ -0.17389837681291273
    @test coord"2 h 12 m 33.33 s" ≈ 0.5783824791568358
    @test_throws LoadError eval(:(coord"02 50 1242.444"))
end
