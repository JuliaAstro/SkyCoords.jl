# SkyCoords.jl benchmarks

Two independent layers:

1. Regression tracking (`benchmarks.jl`): a [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl) suite run by [AirspeedVelocity.jl](https://github.com/MilesCranmer/AirspeedVelocity.jl) on every pull request.
2. Cross-language comparison (`astropy/`): the same benchmarks timed against [astropy.coordinates](https://docs.astropy.org/en/stable/coordinates/), run by hand.

The `bench/` folder at the repository root is unrelated to both. It holds the older standalone comparison behind the plot in the documentation.

## The suite

`benchmarks.jl` defines `SUITE` with two top-level groups.

### `astropy`

Mirrors [`benchmarks/coordinates.py`](https://github.com/astropy/astropy-benchmarks/blob/main/benchmarks/coordinates.py) from astropy-benchmarks, the source of the coordinates section of the [astropy benchmark dashboard](https://www.astropy.org/astropy-benchmarks/). Names are upstream's without the `time_` prefix, so `astropy/frame/init_scalar` corresponds to `coordinates.FrameBenchmarks.time_init_scalar`. Input sizes and values follow upstream's `setup` methods.

| Upstream benchmark | Here | n |
|:--|:--|--:|
| `FrameBenchmarks.time_init_scalar` | `FK5Coords{2000}(ra, dec)` from quantities | 1 |
| `FrameBenchmarks.time_init_array` | `FK5Coords{2000}.(ras, decs)` from quantities | 10³ |
| `FrameBenchmarks.time_concatenate_scalar` | `[c, c]` | 1 |
| `FrameBenchmarks.time_concatenate_array` | `vcat(cs, cs)` | 10⁴ |
| `FrameBenchmarks.time_coord_match_uniform` | `SkyCoords.match(catalog, coords)`, directions uniform in the unit cube | 10⁴ |
| `FrameBenchmarks.time_coord_match_sphere` | same, uniform on the sphere | 10⁴ |
| `FrameBenchmarks.time_coord_match_clusters` | same, eight Gaussian clusters | 8×10⁴ |
| `SkyCoordBenchmarks.time_init_scalar` | `ICRSCoords(ra, dec)` from plain radians | 1 |
| `SkyCoordBenchmarks.time_init_array` | `ICRSCoords.(lon, lat)` from plain radians | 10⁶ |
| `SkyCoordBenchmarks.time_init_quantity_scalar_positional` | `ICRSCoords(ra, dec)` from quantities | 1 |
| `SkyCoordBenchmarks.time_init_quantity_array_positional` | `ICRSCoords.(ras, decs)` from quantities | 10⁶ |
| `SkyCoordBenchmarks.time_init_quantity_scalar_positional_fk5_kwarg` and `..._fk5_frame` | `FK5Coords{2010}(ra, dec)` from quantities | 1 |
| `SkyCoordBenchmarks.time_init_repr_scalar_noframe` | `ICRSCoords(c)` with `c` a `CartesianCoords` | 1 |
| `SkyCoordBenchmarks.time_init_repr_array_noframe` | `ICRSCoords.(cs)` with `cs` a vector of `CartesianCoords` | 10⁶ |
| `SkyCoordBenchmarks.time_repr_scalar` | `repr(c)` | 1 |
| `SkyCoordBenchmarks.time_repr_array` | `repr(cs)` | 10³ |
| `SkyCoordBenchmarks.time_icrs_to_galactic_scalar` | `convert(GalCoords, c)` | 1 |
| `SkyCoordBenchmarks.time_icrs_to_galactic_array` | `convert.(GalCoords, cs)` | 10⁶ |

Some notes on the mapping:

- Benchmarks that take quantities appear twice, as `.../Unitful` and `.../DynamicQuantities`, since each package goes through its own extension.
- Upstream's two FK5 J2010 benchmarks differ in how the frame is passed (`frame="fk5", equinox=J2010` versus an `FK5` instance). Here the equinox is a type parameter, so both are the same call.
- The `init_repr_*_noframe` pair is the loosest match. Upstream wraps the angles of a spherical representation, while the only representation object here is Cartesian, so these convert from Cartesian to spherical.
- The match catalogs are stored as vectors of `ICRSCoords`, the usual form of a catalog here, rather than in upstream's Cartesian representation.
- The benchmarks over 10⁶ coordinates write their result into a preallocated array (`out .= ICRSCoords.(lon, lat)`). Whether a fresh 16 MB result gets page-faulted in depends on the state of the allocator, which made the cheaper of these benchmarks bimodal (3 ms or 11 ms for the same code).
- `concatenate_scalar` builds a two-element `Vector`, so it times the allocator more than anything in this package. Expect it to be the noisiest row.
- AirspeedVelocity adds a `time_to_load` measurement on its own, the counterpart of upstream's `imports.timeraw_import_astropy_coordinates`.

Twelve of the 31 upstream benchmarks have no counterpart:

| Upstream benchmark | Why not |
|:--|:--|
| `time_latitude`, `time_angle_array_repr`, `time_angle_array_str`, `time_angle_array_repr_latex` | No angle types here. Angles are plain numbers, and the latitude check runs inside every constructor timed above. |
| `RepresentationBenchmarks.time_with_differentials_scalar`, `..._array`, `FrameBenchmarks.time_init_scalar_diff` | No velocities or proper motions. |
| `FrameBenchmarks.time_init_nodata` | Frames are types, so there is no empty frame to construct. |
| `SkyCoordBenchmarks.time_init_quantity_scalar_positional_fk5_frame_extra_kwargs` | No observation time attribute. |
| `SkyCoordBenchmarks.time_init_quantity_scalar_keyword`, `..._array_keyword` | No keyword constructors. |
| `SkyCoordBenchmarks.time_iter_array` | Iterating a `Vector` of coordinates compiles down to nothing. |

### `skycoords`

Operations upstream does not cover:

- `transform`: ICRS to FK5 (J2000 and J1975), supergalactic, ecliptic, and FK4, plus FK4 back to ICRS, as scalars and for two of them over 10⁵ coordinates.
- `geometry`: `separation` (same frame, mixed frames, and over 10⁵ pairs), `position_angle`, and `offset`.
- `representation`: `cartesian`, `SkyCoords.project`, and the conversion back out of a projected coordinate.
- `altaz`: ICRS to horizontal coordinates and back through the SOFA.jl extension, against a fresh frame, against a frame whose astrometry context is already cached, and over 10⁴ coordinates.

### Adding benchmarks

- AirspeedVelocity does not tune the suite, so anything faster than a few microseconds needs an explicit `evals` and should read its arguments through `Ref`s (`$(Ref(c))[]`). Otherwise the timer granularity swamps the measurement, or the compiler folds the call away.
- Write target types out in the benchmarked expression (`convert(GalCoords, $c[])`) instead of interpolating them. An interpolated type costs a dynamic dispatch per call.
- CI runs the suite of the pull request against the default branch as well. A benchmark of a function that does not exist there yet needs a guard such as `isdefined(SkyCoords, :newfunction)`, or the base run fails.
- A package the suite needs goes into `benchmark/Project.toml`.

## CI

`.github/workflows/Benchmark.yml` runs the suite on the pull request head and on `main`, and writes the comparison tables (time and memory, with ratios) to the summary page of the workflow run. A regression shows up as a ratio above 1.

The suite and its `Project.toml` are taken from the pull request (`bench-on`), so a change to the benchmarks is exercised by the pull request that makes it. Pull requests from forks are benchmarked as well.

## Running the suite locally

Quick run with a summary table. This directory's environment is a workspace member, so the parent SkyCoords checkout is used automatically:

```sh
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/benchmarks.jl
```

To compare the working tree, uncommitted changes included, against `main`, use the `benchpkg` command that AirspeedVelocity installs:

```sh
julia -e 'using Pkg; Pkg.add("AirspeedVelocity"); Pkg.build("AirspeedVelocity")'   # once
mkdir -p benchmark/results
benchpkg SkyCoords --rev=main,dirty --path=. --script=benchmark/benchmarks.jl \
    --add=NearestNeighbors,Unitful,DynamicQuantities,SOFA,StableRNGs \
    --output-dir=benchmark/results/
```

`--add` is needed here because an explicit `--script` makes AirspeedVelocity skip `benchmark/Project.toml`. Results land in `benchmark/results/`, which is gitignored.

## Comparison with astropy

```sh
julia -t 1 --project=benchmark/astropy -e 'using Pkg; Pkg.instantiate()'
julia -t 1 --project=benchmark/astropy benchmark/astropy/compare.jl
```

The first run resolves a pinned Conda environment (see `astropy/CondaPkg.toml`). The script then:

1. Includes `benchmarks.jl` and runs its `astropy` group, so the Julia numbers are those of the suite CI tracks.
2. Hands the same inputs to Python and checks that both sides agree before timing anything: ICRS to Galactic to within 0.0001 arcsec, and identical matches for the three catalog pairs.
3. Times the upstream statement of each benchmark in Python (`astropy/astropy_cases.py`).
4. Writes `astropy/comparison.csv` and `astropy/provenance.toml` (versions, CPU, date) and prints the table. Both files are gitignored. A ratio above 1 means SkyCoords.jl is faster.

Three things to keep in mind when reading the ratios:

- Catalog matching is timed with `storekdtree=False`. Upstream's statement stores the KD-tree on the catalog after the first call, so its steady state times the query alone, whereas `SkyCoords.match(catalog, coords)` builds the tree on every call. Turning the cache off makes both sides do the same work.
- The `init_repr_*_noframe` pair is left out, since the two sides do different work there (see the notes on the mapping above).
- For the benchmarks over 10⁶ coordinates, the Julia side writes into a preallocated array (also explained above) while astropy allocates its result. Those rows favor SkyCoords.jl by the cost of that allocation, which measured anywhere from nothing to about 20 ms here.

Methodology: single-threaded on both sides (`julia -t 1`, and `OPENBLAS`/`OMP`/`MKL`/`NUMEXPR_NUM_THREADS=1`). Julia is timed with BenchmarkTools, taking the median. Python is timed with `time.perf_counter` around an inner loop sized to at least 1 ms per trial, with 2 warmup trials discarded and the median of the next 7 taken.

This environment needs Julia 1.11 or later for the `[sources]` entry that points it at the parent checkout.
