# Cross-language benchmark: SkyCoords.jl vs astropy.coordinates, on the
# benchmarks of astropy-benchmarks' coordinates.py that have a counterpart here.
#
#     julia -t 1 --project=benchmark/astropy -e 'using Pkg; Pkg.instantiate()'
#     julia -t 1 --project=benchmark/astropy benchmark/astropy/compare.jl
#
# Writes comparison.csv and provenance.toml next to this script (gitignored).
#
# The Julia side is the "astropy" group of ../benchmarks.jl itself, so the
# numbers are those of the suite CI tracks. The astropy side (astropy_cases.py)
# times the corresponding upstream statements on the same inputs.

# Single-threaded numerics on both sides. NumPy reads these at import time, so
# they must be set before PythonCall initializes Python.
ENV["OPENBLAS_NUM_THREADS"] = "1"
ENV["OMP_NUM_THREADS"] = "1"
ENV["MKL_NUM_THREADS"] = "1"
ENV["NUMEXPR_NUM_THREADS"] = "1"
# Keep ~/.local site-packages out of sys.path: a user-site numpy or astropy
# would silently shadow the pinned CondaPkg environment.
ENV["PYTHONNOUSERSITE"] = "1"

using BenchmarkTools
using Dates: now
using Printf: @printf
using PythonCall
using Statistics: median
using TOML: TOML

# Defines SUITE and the input data
include(joinpath(@__DIR__, "..", "benchmarks.jl"))

const np = pyimport("numpy")
pyimport("sys").path.insert(0, @__DIR__)
const astropy_cases = pyimport("astropy_cases")

# -----------------------------------------------------------------------------
# Provenance

const versions = pyconvert(Dict{String, String}, astropy_cases.versions())
const provenance = Dict(
    "date" => string(now()),
    "julia" => string(VERSION),
    "julia_threads" => Threads.nthreads(),
    "cpu" => Sys.cpu_info()[1].model,
    "SkyCoords" => string(pkgversion(SkyCoords)),
    "Unitful" => string(pkgversion(Unitful)),
    "DynamicQuantities" => string(pkgversion(DynamicQuantities)),
    "NearestNeighbors" => string(pkgversion(NearestNeighbors)),
    versions...,
)
# provenance.toml is written at the very end, after the CSV, so a partial run
# can never pair fresh metadata with stale numbers from an earlier run.

# -----------------------------------------------------------------------------
# Inputs, handed to Python as NumPy-owned copies. Angles go over in degrees and
# catalogs as unit vectors, the forms upstream builds its objects from.

pyarray(x) = np.array(x; dtype = np.float64)
lon_deg(cs) = pyarray(rad2deg.(SkyCoords.lon.(cs)))
lat_deg(cs) = pyarray(rad2deg.(SkyCoords.lat.(cs)))
unit_vectors(cs) = pyarray(stack(c -> vec(cartesian(c)), cs))

const data = pydict(
    Dict(
        "frame_array_ra_deg" => pyarray(FRAME_ARRAY_RA),
        "frame_array_dec_deg" => pyarray(FRAME_ARRAY_DEC),
        "icrs_1e4_ra_deg" => lon_deg(ICRS_1E4),
        "icrs_1e4_dec_deg" => lat_deg(ICRS_1E4),
        "lon_1e6_deg" => pyarray(rad2deg.(LON_1E6)),
        "lat_1e6_deg" => pyarray(rad2deg.(LAT_1E6)),
        "array_ra_deg" => pyarray(ARRAY_RA_DEG),
        "array_dec_deg" => pyarray(ARRAY_DEC_DEG),
        (
            name => pytuple((unit_vectors(coords), unit_vectors(catalog)))
                for (name, (coords, catalog)) in MATCH_CATALOGS
        )...,
    )
)

# -----------------------------------------------------------------------------
# Cross-validation: both sides must agree before timing means anything

# ICRS --> Galactic, to the 0.0001 arcsec the documentation claims
let coords = CATALOG_1E4
    l, b = astropy_cases.galactic(lon_deg(coords), lat_deg(coords))
    theirs = GalCoords.(pyconvert(Vector{Float64}, l), pyconvert(Vector{Float64}, b))
    worst = maximum(separation.(convert.(GalCoords, coords), theirs))
    @assert worst < deg2rad(1.0e-4 / 3600) "ICRS --> Galactic differs from astropy by $(rad2deg(worst) * 3600) arcsec"
    @info "Cross-validation passed: ICRS --> Galactic" max_difference_arcsec = rad2deg(worst) * 3600
end

# Catalog matching: the same neighbors at the same separations
for (name, (coords, catalog)) in MATCH_CATALOGS
    id, sep = SkyCoords.match(catalog, coords)
    idx, their_sep = astropy_cases.match(data[name]...)
    @assert id == pyconvert(Vector{Int}, idx) .+ 1 "$name: matched indices differ from astropy"
    worst = maximum(abs, sep .- pyconvert(Vector{Float64}, their_sep))
    @assert worst < 1.0e-10 "$name: separations differ from astropy by $worst rad"
    @info "Cross-validation passed: $name" max_difference_rad = worst
end

# -----------------------------------------------------------------------------
# Timing

# Mirrored benchmarks whose two sides do different work, so that a ratio would
# say nothing: upstream wraps the angles of a spherical representation, where
# SkyCoords converts from its Cartesian one.
const NOT_COMPARED = ("init_repr_scalar_noframe", "init_repr_array_noframe")

@info "Benchmarking SkyCoords.jl"
const julia_results = let suite = deepcopy(SUITE["astropy"])
    foreach(name -> delete!(suite["skycoord"], name), NOT_COMPARED)
    run(suite; verbose = true)
end

@info "Benchmarking astropy"
const astropy_ns = pyconvert(Dict{String, Float64}, astropy_cases.run(data))

# Number of coordinates each benchmark works on, where it is not a single one
const SIZES = Dict(
    "frame/init_array" => 10^3,
    "frame/concatenate_array" => 10^4,
    "frame/coord_match_uniform" => 10^4,
    "frame/coord_match_sphere" => 10^4,
    "frame/coord_match_clusters" => 8 * 10^4,
    "skycoord/init_array" => 10^6,
    "skycoord/init_quantity_array_positional" => 10^6,
    "skycoord/repr_array" => 10^3,
    "skycoord/icrs_to_galactic_array" => 10^6,
)

# A benchmark that goes through a units package has that package's name as a
# third path component. Each such variant is compared with the one astropy case.
const rows = map(sort(BenchmarkTools.leaves(julia_results); by = first)) do (path, trial)
    case = join(path[1:2], "/")
    julia_ns = median(trial).time
    (;
        benchmark = join(path, "/"), n = get(SIZES, case, 1),
        julia_ns, astropy_ns = astropy_ns[case], ratio = astropy_ns[case] / julia_ns,
    )
end

# The two files define the same set of cases, so neither can drift unnoticed
let ours = Set(join(path[1:2], "/") for (path, _) in BenchmarkTools.leaves(julia_results))
    @assert ours == keys(astropy_ns) "benchmarks.jl and astropy_cases.py disagree on: $(symdiff(ours, keys(astropy_ns)))"
end

open(joinpath(@__DIR__, "comparison.csv"), "w") do io
    println(io, "benchmark,n,julia_ns,astropy_ns,ratio")
    for r in rows
        println(io, join((r.benchmark, r.n, r.julia_ns, r.astropy_ns, r.ratio), ','))
    end
end
open(joinpath(@__DIR__, "provenance.toml"), "w") do io
    TOML.print(io, provenance; sorted = true)
end

# -----------------------------------------------------------------------------
# Console summary

# Three significant digits, without a trailing ".0" on the large ratios
function prettyratio(x)
    r = round(x; sigdigits = 3)
    return string(r >= 100 ? round(Int, r) : r, "×")
end

println("\nMedian times. A ratio above 1 means SkyCoords.jl is faster.\n")
width = maximum(r -> length(r.benchmark), rows)
@printf "%-*s  %9s  %12s  %12s  %8s\n" width "benchmark" "n" "SkyCoords.jl" "astropy" "ratio"
for r in rows
    @printf(
        "%-*s  %9d  %12s  %12s  %8s\n", width, r.benchmark, r.n,
        BenchmarkTools.prettytime(r.julia_ns), BenchmarkTools.prettytime(r.astropy_ns),
        prettyratio(r.ratio),
    )
end
println("\nastropy $(versions["astropy"]), numpy $(versions["numpy"]), Python $(versions["python"]), Julia $VERSION")
