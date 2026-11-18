"""astropy side of compare.jl.

Each case is the statement timed by the benchmark of the same name in
benchmarks/coordinates.py of astropy-benchmarks
(https://github.com/astropy/astropy-benchmarks), limited to the ones that have
a SkyCoords.jl counterpart in ../benchmarks.jl. The inputs are handed over from
Julia, so both sides work on the same numbers.
"""

import sys
import time

import astropy
import numpy as np
import scipy
from astropy import units as u
from astropy.coordinates import (
    FK5,
    ICRS,
    CartesianRepresentation,
    SkyCoord,
    concatenate,
    match_coordinates_sky,
)
from astropy.time import Time

WARMUP, TRIALS = 2, 7


def bench(call):
    """Median time of `call` in ns, over trials of at least about 1 ms each."""
    call()  # first-call caches
    t0 = time.perf_counter()
    call()
    estimate = time.perf_counter() - t0
    k = max(1, int(1e-3 / max(estimate, 1e-9)))
    times = []
    for i in range(WARMUP + TRIALS):
        t0 = time.perf_counter()
        for _ in range(k):
            call()
        dt = (time.perf_counter() - t0) / k
        if i >= WARMUP:
            times.append(dt)
    return 1e9 * float(np.median(times))


def cartesian_icrs(xyz):
    """ICRS frame in a Cartesian representation, as upstream builds its catalogs."""
    return ICRS(xyz * u.kpc, representation_type=CartesianRepresentation)


def match(coords_xyz, catalog_xyz):
    """Index into the catalog, and separation in radians, of each coordinate's match.

    Upstream calls match_coordinates_sky(coords, catalog), which stores the
    KD-tree on the catalog after the first call, so its steady state times the
    query alone. Turning the cache off makes astropy build the tree on every
    call, which is what SkyCoords.match(catalog, coords) does.
    """
    idx, sep, _ = match_coordinates_sky(
        cartesian_icrs(coords_xyz), cartesian_icrs(catalog_xyz), storekdtree=False
    )
    return idx, sep.rad


def galactic(ra_deg, dec_deg):
    """Galactic (l, b) in radians of the given ICRS coordinates."""
    g = SkyCoord(ra_deg, dec_deg, unit="deg", frame="icrs").transform_to("galactic")
    return g.l.rad, g.b.rad


def versions():
    return {
        "astropy": astropy.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "python": sys.version.split()[0],
    }


def make_cases(data):
    j2010 = Time("J2010")

    # FrameBenchmarks
    scalar_ra, scalar_dec = 3.2 * u.deg, 2.2 * u.deg
    array_ra = data["frame_array_ra_deg"] * u.deg
    array_dec = data["frame_array_dec_deg"] * u.deg
    icrs_scalar = ICRS(ra=1 * u.deg, dec=2 * u.deg)
    icrs_array = ICRS(
        ra=data["icrs_1e4_ra_deg"] * u.deg, dec=data["icrs_1e4_dec_deg"] * u.deg
    )
    cases = {
        "frame/init_scalar": lambda: FK5(scalar_ra, scalar_dec),
        "frame/init_array": lambda: FK5(array_ra, array_dec),
        "frame/concatenate_scalar": lambda: concatenate((icrs_scalar, icrs_scalar)),
        "frame/concatenate_array": lambda: concatenate((icrs_array, icrs_array)),
    }
    for name in ("coord_match_uniform", "coord_match_sphere", "coord_match_clusters"):
        coords, catalog = (cartesian_icrs(xyz) for xyz in data[name])
        # See `match` for why the KD-tree cache is turned off
        cases["frame/" + name] = lambda coords=coords, catalog=catalog: (
            match_coordinates_sky(coords, catalog, storekdtree=False)
        )

    # SkyCoordBenchmarks
    coord_scalar = SkyCoord(1, 2, unit="deg", frame="icrs")
    lon_1e6, lat_1e6 = data["lon_1e6_deg"], data["lat_1e6_deg"]
    coord_array_1e3 = SkyCoord(lon_1e6[:1000], lat_1e6[:1000], unit="deg", frame="icrs")
    coord_array_1e6 = SkyCoord(lon_1e6, lat_1e6, unit="deg", frame="icrs")
    scalar_q_ra, scalar_q_dec = 1 * u.deg, 2 * u.deg
    array_q_ra = data["array_ra_deg"] * u.deg
    array_q_dec = data["array_dec_deg"] * u.deg
    # The two init_repr_*_noframe benchmarks are left out: see NOT_COMPARED in
    # compare.jl.
    cases.update(
        {
            "skycoord/init_scalar": lambda: SkyCoord(1, 2, unit="deg", frame="icrs"),
            "skycoord/init_array": lambda: SkyCoord(
                lon_1e6, lat_1e6, unit="deg", frame="icrs"
            ),
            "skycoord/init_quantity_scalar_positional": lambda: SkyCoord(
                scalar_q_ra, scalar_q_dec
            ),
            "skycoord/init_quantity_array_positional": lambda: SkyCoord(
                array_q_ra, array_q_dec
            ),
            "skycoord/init_quantity_scalar_positional_fk5": lambda: SkyCoord(
                scalar_q_ra, scalar_q_dec, frame="fk5", equinox=j2010
            ),
            "skycoord/repr_scalar": lambda: repr(coord_scalar),
            "skycoord/repr_array": lambda: repr(coord_array_1e3),
            "skycoord/icrs_to_galactic_scalar": lambda: coord_scalar.transform_to(
                "galactic"
            ),
            "skycoord/icrs_to_galactic_array": lambda: coord_array_1e6.transform_to(
                "galactic"
            ),
        }
    )
    return cases


def run(data):
    """Median time in ns of every case, keyed by benchmark name."""
    results = {}
    for name, call in make_cases(data).items():
        results[name] = bench(call)
        print(f"  {name:45s} {results[name]:16.1f} ns", flush=True)
    return results
