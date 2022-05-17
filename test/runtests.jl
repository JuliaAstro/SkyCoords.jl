#!/usr/bin/env julia

# Tests against astropy.

using AstroAngles
using DelimitedFiles
using Printf
using SkyCoords
using Statistics
using Test
using ConstructionBase: setproperties

import SkyCoords: lat, lon

datapath = joinpath(dirname(@__FILE__), "data")

rad2arcsec(r) = 3600 * rad2deg(r)

# input coordinates
fname = joinpath(datapath, "input_coords.csv")
indata, inhdr = readdlm(fname, ','; header = true)

# Float32 has a large tolerance compared to Float64 and BigFloat, but here we
# are more interested in making sure that the infrastructure works for different
# floating types.
@testset "Testing $F" for (F, TOL) in (
    (Float32, 0.2),
    (Float64, 0.0001),
    (BigFloat, 0.0001),
)
    for (insys, T) in (
        ("icrs", ICRSCoords{F}),
        ("fk5j2000", FK5Coords{2000,F}),
        ("fk5j1975", FK5Coords{1975,F}),
        ("gal", GalCoords{F}),
    )
        c_in = T[T(indata[i, 1], indata[i, 2]) for i = 1:size(indata, 1)]
        for c in c_in
            @test convert(T, c) == c
        end

        for (outsys, S) in (
            ("icrs", ICRSCoords{F}),
            ("fk5j2000", FK5Coords{2000,F}),
            ("fk5j1975", FK5Coords{1975,F}),
            ("gal", GalCoords{F}),
        )
            (outsys == insys) && continue

            c_out = S[convert(S, c) for c in c_in]

            # Test pipe and constructor conversion
            @test c_out == S[S(c) for c in c_in]
            @test c_out == S[c |> S for c in c_in]

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
    @test separation(c1, c5) ≈ separation(c5, c1) ≈ separation(c1, convert(GalCoords, c5)) ≈
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

@testset "string construction" for C in [
    ICRSCoords,
    GalCoords,
    FK5Coords{2000},
    FK5Coords{1970},
]
    @test C(hms"0h0m0", dms"0d0m0") == C(0.0, 0.0)
    @test C(hms"12h0.0m0.0s", dms"90:0:0") == C(π, π / 2)
    @test C(hms"18h0:0", dms"90:0:0") == C(3π / 2, π / 2)
    @test C(hms"12:0:0", dms"90:0:0") == C(π, π / 2)
end

# Test separation between coordinates and conversion with mixed floating types.
@testset "Position Angles" begin
    c1 = ICRSCoords(0, 0)
    c2 = ICRSCoords(deg2rad(1), 0)

    # interface
    @test @inferred position_angle(c1, c2) ≈ @inferred position_angle(c1, c2 |> GalCoords)
    @test position_angle(c1, c2) ≈ position_angle(c1, c2 |> GalCoords)
    
    # accuracy
    @test position_angle(c1, c2) ≈ π / 2

    c3 = ICRSCoords(deg2rad(1), deg2rad(0.1))
    @test position_angle(c1, c3) < π / 2

    c4 = ICRSCoords(0, deg2rad(1))
    @test position_angle(c1, c4) ≈ 0

    # types
    for T in [ICRSCoords, GalCoords, FK5Coords{2000}]
        c1 = T(0, 0)
        c2 = T(deg2rad(1), 0)
        @test position_angle(c1, c2) ≈ π / 2
    end
end



@testset "Offset ($T1, $T2)" for T1 in [ICRSCoords, GalCoords, FK5Coords{2000}], T2 in [ICRSCoords, GalCoords, FK5Coords{2000}]
    # simple integration tests, depend that separation and position_angle are accurate
    c1s = [
        T1(0, -π/2), # south pole
        T1(0, π/2), # north pole
        T1(deg2rad(1), deg2rad(2))
    ]
    c2 = T2(deg2rad(5), deg2rad(10))

    for c1 in c1s
        sep, pa = @inferred offset(c1, c2)
        test_c2 = @inferred offset(c1, sep, pa)
        @test test_c2 isa T1
        test_c2 = T2(test_c2) 
        @test lon(test_c2) ≈ lon(c2)
        @test lat(test_c2) ≈ lat(c2)
    end

    # specific cases to cover special cases.
    c1 = T1(0, deg2rad(89))
    for (pa, sep) in [(0, 2), (180, 358)]
        sep = deg2rad(sep)
        pa = deg2rad(pa)
        c2 = offset(c1, sep, pa)
        @test lon(c2) |> rad2deg ≈ 180
        @test lat(c2) |> rad2deg ≈ 89

        c2 = offset(c1, 2sep, pa)
        @test lon(c2) |> rad2deg ≈ 180
        @test lat(c2) |> rad2deg ≈ 87
    end

    # verify antipode
    c1 = T1(deg2rad(10), deg2rad(47))
    for pa in range(0, stop=377, length=10)
        c2 = offset(c1, deg2rad(180), deg2rad(pa))
        @test lon(c2) |> rad2deg ≈ 190
        @test lat(c2) |> rad2deg ≈ -47

        c2 = offset(c1, deg2rad(360), deg2rad(pa))
        @test lon(c2) |> rad2deg ≈ 10
        @test lat(c2) |> rad2deg ≈ 47
    end

    c1 = T1(deg2rad(10), deg2rad(60))
    c2 = offset(c1, deg2rad(1), deg2rad(90))
    @test 11.9 < lon(c2) |> rad2deg < 12.0
    @test 59.9 < lat(c2) |> rad2deg < 60.0
end

@testset "constructionbase" begin
    @test setproperties(ICRSCoords(1, 2), ra=3) == ICRSCoords(3, 2)
    @test setproperties(GalCoords(1, 2), l=3) == GalCoords(3, 2)
    @test setproperties(FK5Coords{2000}(1, 2), ra=3) == FK5Coords{2000}(3, 2)
end
