# BenchmarkTools suite consumed by AirspeedVelocity.jl (benchpkg / CI action).
#
# When included (by AirspeedVelocity or by astropy/compare.jl) this file only
# defines `SUITE`. Run it directly to also execute the suite and print a
# summary table:
#
#     julia --project=benchmark benchmark/benchmarks.jl
#
# Two top-level groups:
#
#   "astropy"    mirrors benchmarks/coordinates.py from astropy-benchmarks
#                (https://github.com/astropy/astropy-benchmarks). Names are
#                upstream's without the `time_` prefix, and data sizes and
#                values follow upstream's `setup` methods. See README.md for
#                the full mapping and for the upstream benchmarks that have no
#                counterpart here.
#   "skycoords"  operations upstream does not cover.
#
# AirspeedVelocity does not tune the suite, so untuned benchmarks run with
# evals = 1, where the timer granularity of a few tens of ns would swamp a
# scalar kernel. Benchmarks faster than a few microseconds therefore set an
# explicit `evals` and read their arguments through Refs, which keeps the
# compiler from constant-folding the call away.

using BenchmarkTools
using DynamicQuantities: DynamicQuantities
using NearestNeighbors: NearestNeighbors # loads the extension behind SkyCoords.match
using SOFA: SOFA # loads the extension behind the AltAz transforms
using SkyCoords
using StableRNGs: StableRNG
using Unitful: Unitful

# Keep the full-suite runtime PR-friendly (about a minute of measurement per revision).
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0

const SUITE = BenchmarkGroup()

# Both unit packages export `@u_str`, so the degree unit of each is taken
# qualified. astropy's Quantity-based benchmarks are mirrored once per package,
# since each goes through its own SkyCoords extension.
const UNIT_BACKENDS = (
    "Unitful" => Unitful.u"°",
    "DynamicQuantities" => DynamicQuantities.u"°",
)

# -----------------------------------------------------------------------------
# astropy/frame: coordinates.FrameBenchmarks

const FRAME = SUITE["astropy"]["frame"] = BenchmarkGroup()

# time_init_scalar, time_init_array: FK5(ra, dec) from Quantity input
const FRAME_ARRAY_RA = collect(range(0.0, 360.0; length = 1000))
const FRAME_ARRAY_DEC = collect(range(-90.0, 90.0; length = 1000))
for (backend, deg) in UNIT_BACKENDS
    ra, dec = Ref(3.2deg), Ref(2.2deg)
    FRAME["init_scalar"][backend] = @benchmarkable FK5Coords{2000}($ra[], $dec[]) evals = 1000
    ras, decs = FRAME_ARRAY_RA .* deg, FRAME_ARRAY_DEC .* deg
    FRAME["init_array"][backend] = @benchmarkable FK5Coords{2000}.($ras, $decs)
end

# time_concatenate_scalar, time_concatenate_array. Both are bound by the
# allocator, the scalar one entirely, so it is the noisiest entry of the suite.
const ICRS_SCALAR = Ref(ICRSCoords(deg2rad(1.0), deg2rad(2.0)))
const ICRS_1E4 = let rng = StableRNG(12345)
    ICRSCoords.(deg2rad.(rand(rng, 10^4)), deg2rad.(rand(rng, 10^4)))
end
FRAME["concatenate_scalar"] = @benchmarkable [$ICRS_SCALAR[], $ICRS_SCALAR[]] evals = 1000
FRAME["concatenate_array"] = @benchmarkable vcat($ICRS_1E4, $ICRS_1E4)

# time_coord_match_{uniform,sphere,clusters}: nearest-neighbor matching of one
# catalog against another, on the three point distributions upstream uses
# (motivated by astropy/astropy#7324). Upstream builds its catalogs from
# Cartesian points; here they are stored the way catalogs usually are, as
# spherical ICRS coordinates in the direction of each point.
direction(x, y, z) = ICRSCoords(atan(y, x), atan(z, hypot(x, y)))
direction(p) = direction(p[1], p[2], p[3])

# Points uniform in the unit cube, i.e., directions filling one octant
uniform_cube(rng, n) = [direction(rand(rng, 3)) for _ in 1:n]

# Points uniform on the sphere
function uniform_sphere(rng, n)
    ϕ = 2π .* rand(rng, n)
    θ = acos.(2 .* rand(rng, n) .- 1)
    return direction.(cos.(ϕ) .* sin.(θ), sin.(ϕ) .* sin.(θ), cos.(θ))
end

# `n` points scattered with unit variance around each of the `centers`
clustered(rng, centers, n) =
    [direction(c .+ randn(rng, 3)) for _ in 1:n for c in eachcol(centers)]

const MATCH_CATALOGS = let rng = StableRNG(42)
    centers = 200 .* rand(rng, 3, 8) .- 100
    (
        "coord_match_uniform" => (uniform_cube(rng, 10^4), uniform_cube(rng, 10^4)),
        "coord_match_sphere" => (uniform_sphere(rng, 10^4), uniform_sphere(rng, 10^4)),
        "coord_match_clusters" => (clustered(rng, centers, 10^4), clustered(rng, centers, 10^4)),
    )
end
for (name, (coords, catalog)) in MATCH_CATALOGS
    # Upstream's match_coordinates_sky(coords, catalog) finds, for each element
    # of `coords`, its nearest neighbor in `catalog`
    FRAME[name] = @benchmarkable SkyCoords.match($catalog, $coords)
end

# -----------------------------------------------------------------------------
# astropy/skycoord: coordinates.SkyCoordBenchmarks

const SKYCOORD = SUITE["astropy"]["skycoord"] = BenchmarkGroup()

# The benchmarks over 10^6 coordinates write into these preallocated arrays.
# Whether a fresh 16 MB result gets page-faulted in depends on the state of the
# allocator, which made the cheaper of them bimodal (3 ms or 11 ms for the same
# code), so the allocation is kept out of what is timed.
const OUT_ICRS_1E6 = fill(ICRSCoords(0.0, 0.0), 10^6)
const OUT_GAL_1E6 = fill(GalCoords(0.0, 0.0), 10^6)

# time_init_scalar, time_init_array: plain numbers, which SkyCoords takes in radians
const SCALAR_RA = Ref(deg2rad(1.0))
const SCALAR_DEC = Ref(deg2rad(2.0))
const LON_1E6 = fill(deg2rad(1.0), 10^6)
const LAT_1E6 = fill(deg2rad(1.0), 10^6)
SKYCOORD["init_scalar"] = @benchmarkable ICRSCoords($SCALAR_RA[], $SCALAR_DEC[]) evals = 1000
SKYCOORD["init_array"] = @benchmarkable $OUT_ICRS_1E6 .= ICRSCoords.($LON_1E6, $LAT_1E6)

# time_init_quantity_{scalar,array}_positional, and the FK5 J2010 variants.
# Upstream's `frame="fk5", equinox=J2010` and `frame=FK5(equinox=J2010)` are the
# same call here, since the equinox is a type parameter.
const ARRAY_RA_DEG, ARRAY_DEC_DEG = let rng = StableRNG(12345)
    360 .* rand(rng, 10^6), 180 .* rand(rng, 10^6) .- 90
end
for (backend, deg) in UNIT_BACKENDS
    ra, dec = Ref(1.0deg), Ref(2.0deg)
    SKYCOORD["init_quantity_scalar_positional"][backend] =
        @benchmarkable ICRSCoords($ra[], $dec[]) evals = 1000
    SKYCOORD["init_quantity_scalar_positional_fk5"][backend] =
        @benchmarkable FK5Coords{2010}($ra[], $dec[]) evals = 1000
    ras, decs = ARRAY_RA_DEG .* deg, ARRAY_DEC_DEG .* deg
    SKYCOORD["init_quantity_array_positional"][backend] =
        @benchmarkable $OUT_ICRS_1E6 .= ICRSCoords.($ras, $decs)
end

# time_init_repr_{scalar,array}_noframe: a coordinate built from a
# representation object. The representation SkyCoords has is Cartesian.
const CART_SCALAR = Ref(cartesian(ICRS_SCALAR[]))
const CART_1E6 = cartesian.(ICRSCoords.(deg2rad.(ARRAY_RA_DEG), deg2rad.(ARRAY_DEC_DEG)))
SKYCOORD["init_repr_scalar_noframe"] = @benchmarkable ICRSCoords($CART_SCALAR[]) evals = 1000
SKYCOORD["init_repr_array_noframe"] = @benchmarkable $OUT_ICRS_1E6 .= ICRSCoords.($CART_1E6)

# time_repr_scalar, time_repr_array
const ICRS_1E3 = ICRSCoords.(LON_1E6[1:1000], LAT_1E6[1:1000])
SKYCOORD["repr_scalar"] = @benchmarkable repr($ICRS_SCALAR[]) evals = 10
SKYCOORD["repr_array"] = @benchmarkable repr($ICRS_1E3)

# time_icrs_to_galactic_scalar, time_icrs_to_galactic_array
const ICRS_1E6 = ICRSCoords.(LON_1E6, LAT_1E6)
SKYCOORD["icrs_to_galactic_scalar"] = @benchmarkable convert(GalCoords, $ICRS_SCALAR[]) evals = 1000
SKYCOORD["icrs_to_galactic_array"] = @benchmarkable $OUT_GAL_1E6 .= convert.(GalCoords, $ICRS_1E6)

# -----------------------------------------------------------------------------
# skycoords: operations upstream does not cover

# A catalog spread over the whole sphere, and two nearby scalar coordinates
const CATALOG_1E5, CATALOG2_1E5 = let rng = StableRNG(2718)
    uniform_sphere(rng, 10^5), uniform_sphere(rng, 10^5)
end
const C1 = Ref(ICRSCoords(0.3, 0.2))
const C2 = Ref(ICRSCoords(0.31, 0.21))

# Frame transforms. ICRS to Galactic is covered by the astropy group. FK5 at
# J2000 and at another equinox, supergalactic, and ecliptic are all rotations;
# FK4 folds in the E-terms of aberration, the one transform that is not, with
# an iterative solve on the way in.
#
# The target types are written out rather than looped over and interpolated: a
# type passed into the benchmark as a value is not specialized on, which would
# add a dynamic dispatch to every call.
const TRANSFORM = SUITE["skycoords"]["transform"] = BenchmarkGroup()
TRANSFORM["scalar"]["icrs_to_fk5j2000"] = @benchmarkable convert(FK5Coords{2000}, $C1[]) evals = 1000
TRANSFORM["scalar"]["icrs_to_fk5j1975"] = @benchmarkable convert(FK5Coords{1975}, $C1[]) evals = 1000
TRANSFORM["scalar"]["icrs_to_supergalactic"] = @benchmarkable convert(SuperGalCoords, $C1[]) evals = 1000
TRANSFORM["scalar"]["icrs_to_ecliptic"] = @benchmarkable convert(EclipticCoords{2000}, $C1[]) evals = 1000
TRANSFORM["scalar"]["icrs_to_fk4b1950"] = @benchmarkable convert(FK4Coords{1950}, $C1[]) evals = 1000
const FK4_SCALAR = Ref(convert(FK4Coords{1950}, C1[]))
TRANSFORM["scalar"]["fk4b1950_to_icrs"] = @benchmarkable convert(ICRSCoords, $FK4_SCALAR[]) evals = 1000
TRANSFORM["array"]["icrs_to_fk5j1975"] = @benchmarkable convert.(FK5Coords{1975}, $CATALOG_1E5)
TRANSFORM["array"]["icrs_to_fk4b1950"] = @benchmarkable convert.(FK4Coords{1950}, $CATALOG_1E5)

# Separations, position angles, and offsets
const GEOMETRY = SUITE["skycoords"]["geometry"] = BenchmarkGroup()
const GAL2 = Ref(convert(GalCoords, C2[]))
GEOMETRY["separation"] = @benchmarkable separation($C1[], $C2[]) evals = 1000
GEOMETRY["separation_mixed_frames"] = @benchmarkable separation($C1[], $GAL2[]) evals = 1000
GEOMETRY["separation_array"] = @benchmarkable separation.($CATALOG_1E5, $CATALOG2_1E5)
GEOMETRY["position_angle"] = @benchmarkable position_angle($C1[], $C2[]) evals = 1000
const OFFSET_SEP = Ref(0.01)
const OFFSET_PA = Ref(0.5)
GEOMETRY["offset"] = @benchmarkable offset($C1[], $OFFSET_SEP[], $OFFSET_PA[]) evals = 1000

# Cartesian and projected representations
const REPRESENTATION = SUITE["skycoords"]["representation"] = BenchmarkGroup()
const PROJECTED = Ref(SkyCoords.project(C1[], C2[]))
REPRESENTATION["cartesian"] = @benchmarkable cartesian($C1[]) evals = 1000
REPRESENTATION["project"] = @benchmarkable SkyCoords.project($C1[], $C2[]) evals = 1000
REPRESENTATION["unproject"] = @benchmarkable convert(ICRSCoords, $PROJECTED[]) evals = 1000

# Horizontal coordinates (SOFA extension). The first conversion against a frame
# computes its star-independent astrometry context and caches it there, so the
# fresh-frame benchmark pays for that setup and the others, which reuse a frame
# warmed here, do not.
const ALTAZ = SUITE["skycoords"]["altaz"] = BenchmarkGroup()
const OBSERVER = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742)
const JD = 2459526.6667
const WARM_FRAME = AltAzFrame(OBSERVER, JD)
const ALTAZ_SCALAR = Ref(AltAzCoords(C1[], WARM_FRAME))
const CATALOG_1E4 = CATALOG_1E5[1:(10^4)]
ALTAZ["icrs_to_altaz_fresh_frame"] = @benchmarkable AltAzCoords($C1[], AltAzFrame($OBSERVER, $JD))
ALTAZ["icrs_to_altaz"] = @benchmarkable AltAzCoords($C1[], $WARM_FRAME) evals = 10
ALTAZ["altaz_to_icrs"] = @benchmarkable ICRSCoords($ALTAZ_SCALAR[], $WARM_FRAME) evals = 10
ALTAZ["icrs_to_altaz_array"] = @benchmarkable AltAzCoords.($CATALOG_1E4, $WARM_FRAME)

# -----------------------------------------------------------------------------
# Run directly: execute the suite and show a summary table

if abspath(PROGRAM_FILE) == @__FILE__
    using PrettyTables: pretty_table

    results = run(SUITE; verbose = true)
    rows = sort(BenchmarkTools.leaves(results); by = first)
    data = permutedims(
        hcat(
            [
                [
                    join(path, "/"),
                    BenchmarkTools.prettytime(median(t).time),
                    BenchmarkTools.prettymemory(median(t).memory),
                    median(t).allocs,
                ] for (path, t) in rows
            ]...
        )
    )
    pretty_table(
        data;
        column_labels = ["Benchmark", "Median Time", "Memory", "Allocs"],
        alignment = [:l, :r, :r, :r],
        fit_table_in_display_vertically = false,
        fit_table_in_display_horizontally = false,
    )
end
