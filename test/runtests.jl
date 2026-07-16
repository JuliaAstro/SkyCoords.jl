using AstroAngles
using Astrometry
using Accessors
using Unitful
using DynamicQuantities: @us_str, uconvert
using ConstructionBase: setproperties
using DelimitedFiles
using LinearAlgebra: normalize
using NearestNeighbors
using Random: randperm
using SkyCoords
using StableRNGs
using Statistics
using Test
import Makie

import SkyCoords: lat, lon

const rng = StableRNG(2000)
rad2arcsec(r) = 3600 * rad2deg(r)

# tests against astropy.coordinates
include("astropy.jl")

# Test separation between coordinates and conversion with mixed floating types.
@testset "Separation" begin
    c1 = ICRSCoords(ℯ, pi / 2)
    c5 = ICRSCoords(ℯ, 1 + pi / 2)
    @test separation(c1, c5) ≈ separation(c5, c1) ≈ separation(c1, convert(GalCoords, c5)) ≈
        separation(convert(FK5Coords{1980}, c5), c1) ≈ 1
    for T in (GalCoords, FK5Coords{2000}, EclipticCoords{2000})
        c2 = convert(T{Float32}, c1)
        c3 = convert(T{Float64}, c1)
        c4 = convert(T{BigFloat}, c1)
        @test typeof(c2) === T{Float32}
        @test typeof(c3) === T{Float64}
        @test typeof(c4) === T{BigFloat}
        @test isapprox(c2, c3, rtol = sqrt(eps(Float32)))
        @test isapprox(c3, c4, rtol = sqrt(eps(Float64)))
        c6 = convert(T, c5)
        @test separation(c3, c6) ≈ separation(c6, c3) ≈ 1
    end
end

@testset "string construction" for C in [
        ICRSCoords,
        GalCoords,
        FK5Coords{2000},
        FK5Coords{1970},
        EclipticCoords{2000},
        EclipticCoords{1970},
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
    for T in [ICRSCoords, GalCoords, FK5Coords{2000}, EclipticCoords{2000}]
        c1 = T(0, 0)
        c2 = T(deg2rad(1), 0)
        @test position_angle(c1, c2) ≈ π / 2
    end
end


@testset "Offset ($T1, $T2)" for T1 in [ICRSCoords, GalCoords, FK5Coords{2000}, EclipticCoords{2000}], T2 in [ICRSCoords, GalCoords, FK5Coords{2000}, EclipticCoords{2000}]
    # simple integration tests, depend that separation and position_angle are accurate
    c1s = [
        T1(0, -π / 2), # south pole
        T1(0, π / 2), # north pole
        T1(deg2rad(1), deg2rad(2)),
    ]
    c2 = T2(deg2rad(5), deg2rad(10))

    for c1 in c1s
        sep, pa = @inferred offset(c1, c2)
        test_c2 = @inferred offset(c1, sep, pa)
        @test test_c2 isa T1
        test_c2 = T2(test_c2)
        @test test_c2 ≈ c2
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
    for pa in range(0, stop = 377, length = 10)
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

@testset "cartesian" begin
    for CT in [ICRSCoords, FK5Coords{2000}, FK5Coords{1975}, EclipticCoords{2000}, EclipticCoords{1975}, GalCoords]
        @test cartesian(CT(0, 0)) |> vec ≈ [1, 0, 0]
        @test cartesian(CT(0, π / 2)) |> vec ≈ [0, 0, 1]
        @test cartesian(CT(π / 2, 0)) |> vec ≈ [0, 1, 0]
        @test spherical(CartesianCoords{CT}(1, 0, 0)) ≈ CT(0, 0)
        @test spherical(CartesianCoords{CT}(normalize([1, 2, 3]))) ≈ CT(atan(2, 1), atan(3, sqrt(5)))

        c = CT(2, 1)
        c3 = cartesian(c)
        @test c === spherical(c)
        @test c3 === cartesian(c3)
        c_conv = convert(ICRSCoords, c)
        c3_conv = convert(CartesianCoords{ICRSCoords}, c3)
        @test c3_conv == CartesianCoords{ICRSCoords}(c3)
        @test cartesian(c_conv) ≈ c3_conv
        @test c_conv ≈ spherical(c3_conv)

        c_conv3 = convert(CartesianCoords{GalCoords}, c_conv)
        c3_conv = convert(CartesianCoords{GalCoords}, c3_conv)
        @test c_conv3 == CartesianCoords{GalCoords}(c_conv)
        @test c3_conv == CartesianCoords{GalCoords}(c3_conv)
        @test c_conv3 ≈ c3_conv
    end

    a = ICRSCoords(1, 2)
    b = GalCoords(1, 2)
    a3 = cartesian(a)
    b3 = cartesian(b)
    @test separation(a, b) ≈ separation(a3, b3) ≈ separation(a, b3) ≈ separation(a3, b)
end

@testset "constructionbase" begin
    @test setproperties(ICRSCoords(1, 2), ra = 3) == ICRSCoords(3, 2)
    @test setproperties(GalCoords(1, 2), l = 3) == GalCoords(3, 2)
    @test setproperties(FK5Coords{2000}(1, 2), ra = 3) == FK5Coords{2000}(3, 2)
    @test setproperties(EclipticCoords{2000}(1, 2), lon = 3) == EclipticCoords{2000}(3, 2)
    @test setproperties(cartesian(ICRSCoords(1, 2)), vec = [1.0, 0, 0]) == cartesian(ICRSCoords(0, 0))
end

VERSION > v"1.9-DEV" && @testset "Accessors" begin
    @testset for c in (ICRSCoords(1, 0.5), FK5Coords{2000}(1, 0.5), GalCoords(1, 0.5), EclipticCoords{2000}(1, 0.5))
        Accessors.test_getset_laws(lon, c, 1.5, 0.123)
        Accessors.test_getset_laws(lat, c, 1.5, 0.123)

        cart = cartesian(c)
        cart1 = @set lat(spherical(cart)) = 0.123
        @test typeof(cart1) == typeof(cart)
        @test lat(spherical(cart1)) ≈ 0.123

        c1 = @set vec(cartesian(c)) = [1.0, 0, 0]
        @test typeof(c1) == typeof(c)
        @test lat(c1) == 0
        @test lon(c1) == 0

        Accessors.test_getset_laws(spherical, c, c1, c; cmp = ≈)
        Accessors.test_getset_laws(cartesian, c, cart1, cart; cmp = ≈)
        Accessors.test_getset_laws(spherical, cart, c1, c; cmp = ≈)
        Accessors.test_getset_laws(cartesian, cart, cart1, cart; cmp = ≈)
    end
end

VERSION > v"1.9-DEV" && @testset "Unitful" begin
    @test ICRSCoords(1u"rad", 0.5) === ICRSCoords(1, 0.5)
    @test GalCoords(1u"rad", 0.5u"rad") === GalCoords(1, 0.5)
    @test FK5Coords{2000}(1u"°", 0.5) === FK5Coords{2000}(deg2rad(1), 0.5)
    @test EclipticCoords{2000}(1u"°", 0.5u"°") === EclipticCoords{2000}(deg2rad(1), deg2rad(0.5))

    @test SkyCoords.lat(u"rad", ICRSCoords(1, 0.5)) === 0.5u"rad"
    @test SkyCoords.lon(u"°", ICRSCoords(1, 0.5)) === rad2deg(1)u"°"

    @test separation(u"rad", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) === 0.7u"rad"
    @test separation(u"°", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) === rad2deg(0.7)u"°"

    @test position_angle(u"rad", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) === Float64(π) * u"rad"
    @test position_angle(u"°", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) === 180.0u"°"

    # offset() works without special Unitful support
    @test offset(ICRSCoords(1, 0.5), 0.1u"rad", 2) === offset(ICRSCoords(1, 0.5), 0.1, 2)
    @test offset(ICRSCoords(1, 0.5), 0.1u"rad", 2u"rad") === offset(ICRSCoords(1, 0.5), 0.1, 2)
    @test offset(ICRSCoords(1, 0.5), 0.1, 100u"°") === offset(ICRSCoords(1, 0.5), 0.1, deg2rad(100))
    @test offset(ICRSCoords(1, 0.5), 0.1u"°", 100u"°") === offset(ICRSCoords(1, 0.5), deg2rad(0.1), deg2rad(100))
end

VERSION > v"1.9-DEV" && @testset "DynamicQuantities" begin
    # Construction from quantities strips to plain radians (=== holds for the underlying Float64 fields).
    @test ICRSCoords(1us"rad", 0.5) === ICRSCoords(1, 0.5)
    @test GalCoords(1us"rad", 0.5us"rad") === GalCoords(1, 0.5)
    @test FK5Coords{2000}(1us"deg", 0.5) === FK5Coords{2000}(deg2rad(1), 0.5)
    @test EclipticCoords{2000}(1us"deg", 0.5us"deg") === EclipticCoords{2000}(deg2rad(1), deg2rad(0.5))

    # Accessors return symbolic quantities in the requested units.
    # Symbolic `Quantity`s compare equal with `==` (not `===`, which is struct identity).
    @test SkyCoords.lat(us"rad", ICRSCoords(1, 0.5)) == 0.5us"rad"
    @test SkyCoords.lon(us"deg", ICRSCoords(1, 0.5)) == rad2deg(1)us"deg"
    @test SkyCoords.lonlat(us"rad", ICRSCoords(1, 0.5)) == (1us"rad", 0.5us"rad")

    @test SkyCoords.separation(us"rad", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) == 0.7us"rad"
    @test SkyCoords.separation(us"deg", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) == rad2deg(0.7)us"deg"

    @test SkyCoords.position_angle(us"rad", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) == Float64(π) * us"rad"
    @test SkyCoords.position_angle(us"deg", ICRSCoords(1, 0.5), ICRSCoords(1, -0.2)) == 180.0us"deg"

    # offset() accepts angular quantities for the separation and position angle
    @test offset(ICRSCoords(1, 0.5), 0.1us"rad", 2) === offset(ICRSCoords(1, 0.5), 0.1, 2)
    @test offset(ICRSCoords(1, 0.5), 0.1us"rad", 2us"rad") === offset(ICRSCoords(1, 0.5), 0.1, 2)
    @test offset(ICRSCoords(1, 0.5), 0.1, 100us"deg") === offset(ICRSCoords(1, 0.5), 0.1, deg2rad(100))
    @test offset(ICRSCoords(1, 0.5), 0.1us"deg", 100us"deg") === offset(ICRSCoords(1, 0.5), deg2rad(0.1), deg2rad(100))
end

@testset "equality" begin
    @testset for T in [ICRSCoords, GalCoords, FK5Coords{2000}, EclipticCoords{2000}, AltAzCoords]
        c1 = T(1.0, 2.0)
        c2 = T(1.0, 2.001)
        c3 = T{Float32}(1.0, 2.0)
        c4 = T{Float32}(1.0, 2.001)
        @test c1 == c1
        @test_broken c1 == c3
        @test c1 ≈ c1
        @test c1 ≈ c3
        @test !(c1 ≈ c2)
        @test !(c1 ≈ c4)
        @test c1 ≈ c2  rtol = 1.0e-3
        @test c1 ≈ c4  rtol = 1.0e-3
    end

    @test_broken (!(ICRSCoords(1, 2) ≈ FK5Coords{2000}(1, 2)); true)
    @test_broken (!(FK5Coords{2000}(1, 2) ≈ FK5Coords{1950}(1, 2)); true)
end

@testset "conversion" begin
    systems = (ICRSCoords, FK5Coords{2000}, FK5Coords{1975}, EclipticCoords{2000}, EclipticCoords{1975}, GalCoords)
    for IN_SYS in systems, OUT_SYS in systems
        coord_in = IN_SYS(rand(rng), rand(rng))
        coord_out = convert(OUT_SYS, coord_in)
        # Test pipe and constructor conversion
        @test coord_out == OUT_SYS(coord_in)
        @test coord_out == coord_in |> OUT_SYS
    end
end

@testset "AltAz" begin
    # construction, promotion, and azimuth normalization
    @test AltAzCoords(0.5, 1) === AltAzCoords{Float64}(0.5, 1.0)
    @test AltAzCoords(0.5f0, 1.0f0) === AltAzCoords{Float32}(0.5, 1.0)
    @test AltAzCoords(0.3, 2π + 0.1).az ≈ 0.1
    @test convert(AltAzCoords{Float32}, AltAzCoords(0.3, 0.4)) === AltAzCoords{Float32}(0.3, 0.4)

    @test Observer(1, 2) === Observer{Float64}(1.0, 2.0, 0.0)
    @test Observer(1.0f0, 2.0f0, 3.0f0) === Observer{Float32}(1.0, 2.0, 3.0)

    # conversions without an observer location and time are undefined
    @test_throws ArgumentError convert(AltAzCoords, ICRSCoords(1, 2))
    @test_throws ArgumentError convert(ICRSCoords, AltAzCoords(1, 2))
    @test_throws ArgumentError AltAzCoords(ICRSCoords(1, 2))
    @test_throws ArgumentError GalCoords(AltAzCoords(1, 2))
    @test_throws ArgumentError convert(EclipticCoords{2000}, AltAzCoords(1, 2))
    @test_throws ArgumentError convert(AltAzCoords, EclipticCoords{2000}(1, 2))

    # cartesian representations work within the horizontal frame, but changing
    # the frame still requires an observer location and time
    cc = convert(CartesianCoords{AltAzCoords}, AltAzCoords(0.3, 1.2))
    @test vec(cc) == vec(cartesian(AltAzCoords(0.3, 1.2)))
    @test convert(AltAzCoords{Float64}, cc) ≈ AltAzCoords(0.3, 1.2)
    @test_throws ArgumentError convert(CartesianCoords{ICRSCoords}, AltAzCoords(0.3, 1.2))
    @test_throws ArgumentError convert(AltAzCoords{Float64}, cartesian(ICRSCoords(1, 2)))

    # M13 observed from Mount Wilson at 2021-11-08T04:00 UTC.
    # Reference values are cross-checked against astropy (see astropy.jl for a
    # comparison with identical Earth orientation parameters).
    mt_wilson = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742)
    jd = 2459526.5 + 4 / 24
    m13 = ICRSCoords(deg2rad(250.423475), deg2rad(36.4613194))
    altaz = AltAzCoords(m13, mt_wilson, jd)
    @test rad2deg(altaz.alt) ≈ 13.357744229945306  rtol = 1.0e-9
    @test rad2deg(altaz.az) ≈ 305.2066098441702  rtol = 1.0e-9

    # round trips through every celestial frame
    @test ICRSCoords(altaz, mt_wilson, jd) ≈ m13  atol = 1.0e-9
    for T in (ICRSCoords, GalCoords, SuperGalCoords, FK5Coords{2000}, EclipticCoords{2000})
        back = T(altaz, mt_wilson, jd)
        @test back isa T
        @test AltAzCoords(back, mt_wilson, jd) ≈ altaz  atol = 1.0e-9
    end

    # refraction lifts the apparent altitude and leaves the azimuth unchanged
    refr = AltAzCoords(m13, mt_wilson, jd; pressure = 820, temperature = 10, relative_humidity = 0.4)
    @test refr.alt > altaz.alt
    @test refr.az ≈ altaz.az
    @test ICRSCoords(refr, mt_wilson, jd; pressure = 820, temperature = 10, relative_humidity = 0.4) ≈ m13  atol = 1.0e-7

    # Earth orientation parameters shift the result
    @test separation(AltAzCoords(m13, mt_wilson, jd; dut1 = 0.5), altaz) > 0

    # converting between two horizontal frames is ambiguous
    @test_throws ArgumentError AltAzCoords(altaz, mt_wilson, jd)
    @test_throws ArgumentError AltAzCoords{Float64}(altaz, mt_wilson, jd)

    # angles between simultaneous observations match the celestial frame
    # (up to differential aberration)
    c2 = ICRSCoords(deg2rad(250.0), deg2rad(36.0))
    altaz2 = AltAzCoords(c2, mt_wilson, jd)
    @test separation(altaz, altaz2) ≈ separation(m13, c2)  atol = 1.0e-6

    # offset and cartesian round trips preserve the (alt, az) field order
    sep, pa = offset(altaz, altaz2)
    @test offset(altaz, sep, pa) ≈ altaz2
    @test spherical(cartesian(altaz)) ≈ altaz
end

VERSION >= v"1.9" && @testset "plotting with Makie" begin
    coo = ICRSCoords(1, 2)

    @test Makie.convert_arguments(Makie.Scatter, coo) == ([Makie.Point(1, 2)],)
    @test Makie.convert_arguments(Makie.Scatter, [coo]) == ([Makie.Point(1, 2)],)
    @test Makie.convert_arguments(Makie.Lines, [coo, coo]) == ([Makie.Point(1, 2), Makie.Point(1, 2)],)
    @test Makie.convert_arguments(Makie.Lines, [coo][1:0]) == ([],)
end

VERSION >= v"1.9" && @testset "Matching ($T1, $T2)" for T1 in [ICRSCoords, GalCoords, FK5Coords{2000}, EclipticCoords{2000}], T2 in [ICRSCoords, GalCoords, FK5Coords{2000}, EclipticCoords{2000}]
    ## data generation
    N = 1000
    lons = 2pi .* rand(rng, N) # (0, 2π)
    lats = pi .* (rand(rng, N) .- 0.5) # (-π, π)
    refcat = T1.(lons, lats)
    tree = KDTree(refcat)
    # Test mixed coordinate input to KDTree
    @test KDTree([refcat[1], refcat[2]]).data ≈ KDTree([refcat[1], convert(T2, refcat[2])]).data
    @test_throws ArgumentError KDTree(T1[]) # empty data throws
    @test_throws ArgumentError nn(tree, T1[]) # empty coords throws
    @test_throws ArgumentError knn(tree, T1[], 2) # empty coords throws
    for n in (1, 10, N)
        # Test single coord
        @test nn(tree, convert(T2, refcat[n]))[1] == n
        # Test multiple coords
        @test nn(tree, convert.(Ref(T2), [refcat[n], refcat[2]]))[1] == [n, 2]
        # Test knn, single coord
        id, sep = knn(tree, convert(T2, refcat[n]), 2)
        @test length(id) == length(sep) == 2
        @test n in id
        # Test knn with multiple coords which returns a vector of vectors
        # First dimension is number of points=3, second is arg k=2
        id, sep = knn(tree, convert.(Ref(T2), [refcat[n], refcat[2], refcat[3]]), 2)
        @test length(id) == length(sep) == 3
        @test all(length(id[i]) .== length(sep[i]) .== 2 for i in eachindex(id, sep))
        @test all((n in id[1], 2 in id[2], 3 in id[3]))
        id, sep = knn(tree, convert.(Ref(T2), [refcat[n], refcat[2], refcat[3]]), 2, true)
        @test (id[1][1] == n) & (id[2][1] == 2) & (id[3][1] == 3) # Order guaranteed by sortres = true
        # Test match
        rr = randperm(rng, n)
        matchcat = convert.(Ref(T2), refcat)[rr]
        id, sep = SkyCoords.match(refcat, matchcat)
        @test id == rr
        @test all(isapprox.(sep, 0; atol = 1.0e-12))
        id2, sep2 = SkyCoords.match(tree, matchcat)
        @test id2 == id
        @test sep == sep2
        # Test with nthneighbor != 1
        id3, sep3 = SkyCoords.match(refcat, matchcat; nthneighbor = 2)
        @test all(id3 .!= id2)
        @test all(sep3 .> sep2)
        kid, ksep = knn(tree, matchcat[n], 2, false)
        # sortres = false does not guarantee order;
        # the second neighbor is whichever one has greater separation
        a = argmax(ksep)
        @test id3[n] == kid[a]
        @test sep3[n] == ksep[a]
        kid, ksep = knn(tree, matchcat[n], 2, true)
        @test id3[n] == kid[2]
        @test sep3[n] == ksep[2]
        # Test with CartesianCoords
        @test SkyCoords.match(refcat, convert.(Ref(CartesianCoords{T2}), refcat)[rr])[1] == rr
        @test SkyCoords.match(convert.(Ref(CartesianCoords{T1}), refcat), convert.(Ref(T2), refcat)[rr])[1] == rr
    end
    ir = inrange(tree, convert(T2, refcat[1]), 0.1)
    @test sort(ir) == sort(knn(tree, convert(T2, refcat[1]), length(ir))[1])
    @test reduce(==, inrange(tree, convert.(Ref(T2), [refcat[1], refcat[1]]), 0.1))
    @test_nowarn KDTree(reshape(refcat, (100, 10)))
    # Test non-vector input, nn
    id1, sep1 = nn(tree, reshape(refcat[1:4], (2, 2)))
    @test size(id1) == size(sep1) == (2, 2)
    id2, sep2 = nn(tree, refcat[1:4])
    @test (vec(id1) == id2) & (vec(sep1) == sep2)
    # Test non-vector input, knn
    id1, sep1 = knn(tree, reshape(refcat[1:4], (2, 2)), 3)
    @test size(id1) == size(sep1) == (2, 2)
    id2, sep2 = knn(tree, refcat[1:4], 3)
    @test (vec(id1) == id2) & (vec(sep1) == sep2)
    # Test non-vector input, match
    rr = randperm(rng, 1000)
    matchcat = refcat[rr]
    @test_throws ArgumentError SkyCoords.match(refcat, ICRSCoords[])
    @test_throws ArgumentError SkyCoords.match(ICRSCoords[], matchcat)
    id1, sep1 = SkyCoords.match(reshape(refcat, (10, 100)), reshape(matchcat, (100, 10)))
    id2, sep2 = SkyCoords.match(refcat, matchcat)
    @test size(id1) == size(sep1) == (100, 10)
    @test size(id2) == size(sep2) == (1000,)
    @test (vec(id1) == id2 == rr) & (vec(sep1) == sep2)
    id1, sep1 = SkyCoords.match(reshape(refcat, (10, 100)), reshape(matchcat, (100, 10)); nthneighbor = 2)
    id2, sep2 = SkyCoords.match(refcat, matchcat)
end
