module SkyCoords

import ConstructionBase: constructorof
using LinearAlgebra: I, dot, norm, normalize
using Rotations: RotX, RotXYZ, RotZYZ
using StaticArrays: SA, SVector

export AbstractSkyCoords,
    ICRSCoords,
    GalCoords,
    SuperGalCoords,
    FK4Coords,
    FK4NoETerms,
    FK5Coords,
    EclipticCoords,
    CartesianCoords,
    ProjectedCoords,
    separation,
    position_angle,
    offset,
    cartesian,
    spherical

include("types.jl")
include("cartesian.jl")
include("projected.jl")

# -----------------------------------------------------------------------------
# Helper functions: Create rotation matrix about a given axis (x, y, z)

# (lon, lat) -> [x, y, z] unit vector
function coords2cart(lon, lat)
    sinlon, coslon = sincos(lon)
    sinlat, coslat = sincos(lat)
    return SVector{3}(coslat * coslon, coslat * sinlon, sinlat)
end

# [x, y, z] unit vector -> (lon, lat)
function cart2coords(v)
    x, y, z = v[begin], v[begin + 1], v[begin + 2]
    lon = atan(y, x)
    xy_norm = hypot(x, y)
    lat = atan(z, xy_norm)
    return lon, lat
end

# -----------------------------------------------------------------------------
# Constant rotation matrices and precession matrix function

# ICRS --> FK5 at J2000 (See USNO Circular 179, section 3.5)
const ICRS_TO_FK5J2000 = let
    eta0 = deg2rad(-19.9 / 3600_000)
    xi0 = deg2rad(9.1 / 3600_000)
    da0 = deg2rad(-22.9 / 3600_000)
    RotXYZ(eta0, -xi0, -da0)
end
const FK5J2000_TO_ICRS = ICRS_TO_FK5J2000'

# FK5J2000 --> Gal
# Note that galactic pole and zeropoint of l are somewhat arbitrary
# and not officially defined (as far as I know). The values below are
# from astropy.coordinates, which includes the following comment:
# | "This gives better consistency with other codes than using the values
# |  from Reid & Brunthaler 2004 and the best self-consistency between FK5
# |  -> Galactic and FK5 -> FK4 -> Galactic. The lon0 angle was found by
# |  optimizing the self-consistency."
const FK5J2000_TO_GAL = let
    ngp_fk5j2000_ra = deg2rad(192.8594812065348)
    ngp_fk5j2000_dec = deg2rad(27.12825118085622)
    lon0_fk5j2000 = deg2rad(122.9319185680026)
    RotZYZ(lon0_fk5j2000 - π, ngp_fk5j2000_dec - π / 2, -ngp_fk5j2000_ra)
end
const GAL_TO_FK5J2000 = FK5J2000_TO_GAL'

# Gal --> ICRS: simply chain through FK5J2000
const GAL_TO_ICRS = FK5J2000_TO_ICRS * GAL_TO_FK5J2000
const ICRS_TO_GAL = GAL_TO_ICRS'


# Gal --> SuperGal
# we use the same parameters as in astropy
sgp_l = deg2rad(47.37)
sgp_b = deg2rad(6.32)
# rotation matrix compared to astropy
const GAL_TO_SUPERGAL = RotZYZ(π / 2, π / 2 - sgp_b, π - sgp_l)
const SUPERGAL_TO_GAL = GAL_TO_SUPERGAL'
# SuperGal --> ICRS: chain through GAL
const SUPERGAL_TO_ICRS = GAL_TO_ICRS * SUPERGAL_TO_GAL
const ICRS_TO_SUPERGAL = SUPERGAL_TO_ICRS'
# SuperGal --> FK5J2000: chain through GAL
const SUPERGAL_TO_FK5J2000 = GAL_TO_FK5J2000 * SUPERGAL_TO_GAL
const FK5J2000_TO_SUPERGAL = SUPERGAL_TO_FK5J2000'
# FK5J2000 --> FK5{epoch}
# Computes the precession matrix from J2000 to the given Julian equinox.
# Expression from from Capitaine et al. 2003 as expressed in the USNO
# Circular 179.  This should match the IAU 2006 standard from SOFA.
const pzeta = SA[2.650545, 2306.083227, 0.2988499, 0.01801828, -0.000005971, -0.0000003173]
const pz = SA[-2.650545, 2306.077181, 1.0927348, 0.01826837, -0.000028596, -0.0000002904]
const ptheta = SA[0.0, 2004.191903, -0.4294934, -0.04182264, -0.000007089, -0.0000001274]

function precess_from_j2000(equinox)
    t = (equinox - 2000) / 100
    tn = 1.0
    zeta = z = theta = 0.0
    for i in eachindex(pzeta, pz, ptheta)
        zeta += pzeta[i] * tn
        z += pz[i] * tn
        theta += ptheta[i] * tn
        tn *= t
    end
    return RotZYZ(deg2rad(z / 3600), -deg2rad(theta / 3600), deg2rad(zeta / 3600))
end

function ecliptic_obliquity(y)
    # https://github.com/JuliaAstro/AstroBase.jl/blob/master/src/EarthAttitude/obliquity.jl
    T = (y - 2000.0) / 100
    obl = @evalpoly(T, 84381.406, -46.836769, -0.0001831, 0.0020034, -0.000000576, -0.0000000434)
    return deg2rad(obl / 3600)
end

# -----------------------------------------------------------------------------
# Type-dependent methods

lon(c::GalCoords) = c.l
lat(c::GalCoords) = c.b
lon(c::SuperGalCoords) = c.l
lat(c::SuperGalCoords) = c.b
lat(c::EclipticCoords) = c.lat
lon(c::EclipticCoords) = c.lon
lon(c::AbstractSkyCoords) = c.ra
lat(c::AbstractSkyCoords) = c.dec

lonlat(c::AbstractSkyCoords) = (lon(c), lat(c))

# Abstract away specific field names (ra, dec vs l, b)
coords2cart(c::AbstractSkyCoords) = coords2cart(lon(c), lat(c))

# Frame transformation of a Cartesian unit vector: `frame_transform(to, from, v)`.
# This is the single primitive behind every frame change in the package
# (`convert` reduces to it, for spherical and Cartesian representations alike).
# For purely rotational systems — every system except FK4Coords — it is simply
# multiplication by the `rotmat` rotation matrix below. Systems whose transform
# is not a rotation (FK4Coords with its position-dependent E-terms) override
# `frame_transform` directly instead of providing a `rotmat`.
function frame_transform(::Type{T}, ::Type{S}, v) where {T <: AbstractSkyCoords, S <: AbstractSkyCoords}
    _checkframe(T)
    return rotmat(T, S) * v
end

# A type names a coordinate frame when supplying an element type makes it
# concrete: `ICRSCoords`, `FK5Coords{2000}`, and `ICRSCoords{Float32}` all do;
# a bare `FK5Coords` or `FK4Coords` (missing its equinox) does not — it
# satisfies no `T <: FK5Coords{e}` for any single `e`, so it would silently
# fall past every `rotmat`/`frame_transform` method and fail deep inside
# `rotmat` with a confusing `MethodError`. Check before the lookup instead.
# Both checks fold away at compile time on the specialized (valid) paths.
_isframe(::Type{TC}) where {TC <: AbstractSkyCoords} =
    isconcretetype(TC) || (TC isa UnionAll && isconcretetype(TC{Float64}))
_checkframe(::Type{TC}) where {TC <: AbstractSkyCoords} =
    _isframe(TC) || throw(ArgumentError(
        "$TC does not fully specify a coordinate system (a required equinox/epoch " *
        "is missing); fully specify its type parameters, e.g. `$TC{1950}`",
    ))

# Rotation matrix between coordinate systems: `rotmat(to, from)`
# Note that all of these return SMatrix{3,3}{Float64}, regardless of
# element type of input coordinates.
rotmat(::Type{<:GalCoords}, ::Type{<:ICRSCoords}) = ICRS_TO_GAL
rotmat(::Type{<:ICRSCoords}, ::Type{<:GalCoords}) = GAL_TO_ICRS
rotmat(::Type{<:ICRSCoords}, ::Type{<:ICRSCoords}) = I
rotmat(::Type{<:GalCoords}, ::Type{<:GalCoords}) = I
rotmat(::Type{<:FK5Coords{e}}, ::Type{<:FK5Coords{e}}) where {e} = I
rotmat(::Type{<:SuperGalCoords}, ::Type{<:SuperGalCoords}) = I
rotmat(::Type{<:GalCoords}, ::Type{<:SuperGalCoords}) = SUPERGAL_TO_GAL
rotmat(::Type{<:SuperGalCoords}, ::Type{<:GalCoords}) = GAL_TO_SUPERGAL
rotmat(::Type{<:SuperGalCoords}, ::Type{<:ICRSCoords}) = ICRS_TO_SUPERGAL
rotmat(::Type{<:ICRSCoords}, ::Type{<:SuperGalCoords}) = SUPERGAL_TO_ICRS


@generated rotmat(::Type{<:EclipticCoords{e}}, ::Type{<:FK5Coords{e}}) where {e} = RotX(-ecliptic_obliquity(e))
@generated rotmat(::Type{<:FK5Coords{e}}, ::Type{<:EclipticCoords{e}}) where {e} = RotX(ecliptic_obliquity(e))
@generated rotmat(::Type{<:EclipticCoords{e}}, ::Type{T}) where {e, T <: AbstractSkyCoords} = rotmat(EclipticCoords{e}, FK5Coords{e}) * rotmat(FK5Coords{e}, T)
@generated rotmat(::Type{T}, ::Type{<:EclipticCoords{e}}) where {e, T <: AbstractSkyCoords} = rotmat(T, FK5Coords{e}) * rotmat(FK5Coords{e}, EclipticCoords{e})
# disambiguation:
@generated rotmat(::Type{<:EclipticCoords{e_to}}, ::Type{<:EclipticCoords{e_from}}) where {e_to, e_from} = rotmat(EclipticCoords{e_to}, FK5Coords{e_to}) * rotmat(FK5Coords{e_to}, EclipticCoords{e_from})

@generated rotmat(::Type{<:FK5Coords{e1}}, ::Type{<:ICRSCoords}) where {e1} =
    precess_from_j2000(e1) * ICRS_TO_FK5J2000
@generated rotmat(::Type{<:FK5Coords{e1}}, ::Type{<:GalCoords}) where {e1} =
    precess_from_j2000(e1) * GAL_TO_FK5J2000
@generated rotmat(::Type{<:FK5Coords{e1}}, ::Type{<:SuperGalCoords}) where {e1} =
    precess_from_j2000(e1) * SUPERGAL_TO_FK5J2000
@generated rotmat(::Type{<:ICRSCoords}, ::Type{<:FK5Coords{e2}}) where {e2} =
    FK5J2000_TO_ICRS * precess_from_j2000(e2)'
@generated rotmat(::Type{<:GalCoords}, ::Type{<:FK5Coords{e2}}) where {e2} =
    FK5J2000_TO_GAL * precess_from_j2000(e2)'
@generated rotmat(::Type{<:SuperGalCoords}, ::Type{<:FK5Coords{e2}}) where {e2} =
    FK5J2000_TO_SUPERGAL * precess_from_j2000(e2)'
@generated rotmat(::Type{<:FK5Coords{e1}}, ::Type{<:FK5Coords{e2}}) where {e1, e2} =
    precess_from_j2000(e1) * precess_from_j2000(e2)'

# -----------------------------------------------------------------------------
# FK4 (with E-terms of aberration)

# `FK4NoETerms{e}` (see types.jl) is a purely-rotational system, so it
# participates in the `rotmat` network like any other type. `FK4Coords` itself
# cannot, because folding the E-terms in/out is a position-dependent correction
# rather than a rotation; it overrides `frame_transform` below instead.

# Newcomb precession, valid within the FK4/Besselian system only
# (FK5 uses the IAU 1976/2000 Julian precession implemented by `precess_from_j2000`).
# Explanatory Supplement to the Astronomical Almanac (Seidelmann, 1992), as
# implemented in astropy's `FK4NoETerms._precession_matrix`.
function newcomb_precession(byear1, byear2)
    t1 = (byear1 - 1850.0) / 1000.0
    dt = (byear2 - 1850.0) / 1000.0 - t1
    dt_over_3600 = dt / 3600

    zeta1 = (0.060 * t1 + 139.720) * t1 + 23035.545
    zeta2 = -0.27 * t1 + 30.240
    zeta = ((17.995 * dt + zeta2) * dt + zeta1) * dt_over_3600

    z2 = 109.480 + 0.39 * t1
    z = ((18.325 * dt + z2) * dt + zeta1) * dt_over_3600

    theta1 = (-0.37 * t1 - 85.29) * t1 + 20051.12
    theta2 = -0.37 * t1 - 42.65
    theta = ((-41.8 * dt + theta2) * dt + theta1) * dt_over_3600

    return RotZYZ(deg2rad(z), -deg2rad(theta), deg2rad(zeta))
end

# Besselian year -> Julian date (365.2421988-day Besselian year)
besselian_to_jd(byear) = 2433282.4235 + (byear - 1950) * 365.2421988

# Mean obliquity of the ecliptic, IAU 1980 (erfa's `obl80`)
function mean_obliquity80(jd)
    t = (jd - 2451545.0) / 36525.0
    obl = @evalpoly(t, 84381.448, -46.8150, -0.00059, 0.001813)
    return deg2rad(obl / 3600)
end

# E-terms of elliptic aberration for a given FK4 equinox (Besselian year).
# A small (dimensionless, ~1e-4) vector correction in Cartesian space arising
# from folding the aberration due to the eccentricity of Earth's orbit
# directly into the FK4 catalog positions. See astropy's `fk4_e_terms`.
function fk4_eterms(byear)
    k = deg2rad(0.0056932)
    jd = besselian_to_jd(byear)
    t = (jd - besselian_to_jd(1950)) / 36525.0
    ek = k * ((-0.000000126 * t - 0.00004193) * t + 0.01673011)
    g = deg2rad((((0.012 * t + 1.65) * t + 6190.67) * t + 1015489.951) / 3600)
    mekcosg = -ek * cos(g)
    o = mean_obliquity80(jd)
    return SA[ek * sin(g), mekcosg * cos(o), mekcosg * sin(o)]
end

# FK4 (with E-terms) -> FK4NoETerms: closed-form removal.
remove_eterms(v, byear) = (e = fk4_eterms(byear); normalize(v - e + dot(e, v) * v))

# FK4NoETerms -> FK4 (with E-terms): fixed-point iteration (matches astropy).
function add_eterms(v, byear)
    e = fk4_eterms(byear)
    rhs = e + v
    w = v
    for _ in 1:10
        w = rhs / (1 + dot(e, w))
    end
    return normalize(w)
end

# B1950 -> J2000 frame rotation, Murray 1989 A&A 218,325 eqn 28.
const B1950_TO_J2000 = SA[
    0.9999256794956877 -0.0111814832204662 -0.0048590038153592
    0.0111814832391717 0.9999374848933135 -0.0000271625947142
    0.0048590037723143 -0.0000271702937440 0.9999881946023742
]

# Correction accounting for FK4 being a slowly-rotating system, per Julian
# century of obstime from B1950 (Murray 1989 eqn 29). We take obstime to equal
# the FK4 equinox, matching astropy's default when no separate obstime is given.
const FK4_CORR = SA[
    -0.0026455262 -1.1539918689 2.1111346190
    1.1540628161 -0.0129042997 0.0236021478
    -2.1112979048 -0.0056024448 0.0102587734
] * 1.0e-6

fk4_B_matrix(byear) = B1950_TO_J2000 + FK4_CORR * ((byear - 1950.0) / 100.0)

rotmat(::Type{<:FK4NoETerms{e}}, ::Type{<:FK4NoETerms{e}}) where {e} = I

@generated rotmat(::Type{<:FK5Coords{2000}}, ::Type{<:FK4NoETerms{e}}) where {e} =
    fk4_B_matrix(e) * newcomb_precession(e, 1950)
@generated rotmat(::Type{<:FK4NoETerms{e}}, ::Type{<:FK5Coords{2000}}) where {e} =
    rotmat(FK5Coords{2000}, FK4NoETerms{e})'

@generated rotmat(::Type{T}, ::Type{<:FK4NoETerms{e}}) where {e, T <: AbstractSkyCoords} =
    rotmat(T, FK5Coords{2000}) * rotmat(FK5Coords{2000}, FK4NoETerms{e})
@generated rotmat(::Type{<:FK4NoETerms{e}}, ::Type{T}) where {e, T <: AbstractSkyCoords} =
    rotmat(FK4NoETerms{e}, FK5Coords{2000}) * rotmat(FK5Coords{2000}, T)
# Disambiguation, and physically correct: precess directly within the
# Besselian/FK4 system rather than round-tripping through FK5
# (the B-matrix correction above is a one-time system-level offset,
# not something to apply twice).
@generated rotmat(::Type{<:FK4NoETerms{e_to}}, ::Type{<:FK4NoETerms{e_from}}) where {e_to, e_from} =
    newcomb_precession(e_from, e_to)

# EclipticCoords and FK4NoETerms are the only two types with a generic
# "fall back to any T <: AbstractSkyCoords" `rotmat` method (both routing
# through FK5Coords). Without an explicit rule for this specific pair, their
# two fallbacks are ambiguous with each other.
@generated rotmat(::Type{<:EclipticCoords{e1}}, ::Type{<:FK4NoETerms{e2}}) where {e1, e2} =
    rotmat(EclipticCoords{e1}, FK5Coords{e1}) * rotmat(FK5Coords{e1}, FK4NoETerms{e2})
@generated rotmat(::Type{<:FK4NoETerms{e1}}, ::Type{<:EclipticCoords{e2}}) where {e1, e2} =
    rotmat(FK4NoETerms{e1}, FK5Coords{2000}) * rotmat(FK5Coords{2000}, EclipticCoords{e2})

# FK4Coords is FK4NoETerms with the E-terms folded in, so its frame transform
# removes/adds the E-terms around a transform through FK4NoETerms. Because
# these methods define the `frame_transform` primitive itself, every `convert`
# — spherical or Cartesian, in either direction — works generically.
frame_transform(::Type{TO}, ::Type{<:FK4Coords{e}}, v) where {TO <: AbstractSkyCoords, e} =
    frame_transform(TO, FK4NoETerms{e}, remove_eterms(v, e))
frame_transform(::Type{<:FK4Coords{e}}, ::Type{FROM}, v) where {e, FROM <: AbstractSkyCoords} =
    add_eterms(frame_transform(FK4NoETerms{e}, FROM, v), e)
# disambiguation for the FK4 -> FK4 pair; same-equinox is the identity
frame_transform(::Type{<:FK4Coords{e_to}}, ::Type{<:FK4Coords{e_from}}, v) where {e_to, e_from} =
    add_eterms(rotmat(FK4NoETerms{e_to}, FK4NoETerms{e_from}) * remove_eterms(v, e_from), e_to)
frame_transform(::Type{<:FK4Coords{e}}, ::Type{<:FK4Coords{e}}, v) where {e} = v

# ------------------------------------------------------------------------------
# Distance between coordinates

function _separation(λ_1, ϕ_1, λ_2, ϕ_2)
    Δλ = λ_2 - λ_1
    sin_Δλ, cos_Δλ = sincos(Δλ)
    sin_ϕ1, cos_ϕ1 = sincos(ϕ_1)
    sin_ϕ2, cos_ϕ2 = sincos(ϕ_2)
    return atan(
        hypot(cos_ϕ2 * sin_Δλ, cos_ϕ1 * sin_ϕ2 - sin_ϕ1 * cos_ϕ2 * cos_Δλ),
        sin_ϕ1 * sin_ϕ2 + cos_ϕ1 * cos_ϕ2 * cos_Δλ,
    )
end

"""
    separation(c1::AbstractSkyCoords, c2::AbstractSkyCoords) -> distance

Return angular separation between two sky coordinates, in radians.

The angular separation is calculated using the [Vincenty formula](http://en.wikipedia.org/wiki/Great-circle_distance),
which is slightly more complex and computationally expensive than some alternatives,
but is stable at at all distances, including the poles and antipodes.
"""
separation(c1::T, c2::T) where {T <: AbstractSkyCoords} =
    _separation(lon(c1), lat(c1), lon(c2), lat(c2))

separation(c1::CartesianCoords{T}, c2::CartesianCoords{T}) where {T <: AbstractSkyCoords} =
    2 * asin(norm(vec(c1) - vec(c2)) / 2)

separation(c1::T1, c2::T2) where {T1 <: AbstractSkyCoords, T2 <: AbstractSkyCoords} =
    separation(c1, convert(T1, c2))


"""
    position_angle(c1::AbstractSkyCoords, c2::AbstractSkyCoords) -> angle

Return position angle between two sky coordinates, in positive radians.

### Examples
```jldoctest
julia> c1 = ICRSCoords(0, 0); c2 = ICRSCoords(deg2rad(1), 0);

julia> position_angle(c1, c2) |> rad2deg
90.0
```
"""
position_angle(c1::T, c2::T) where {T <: AbstractSkyCoords} = _position_angle(lon(c1), lat(c1), lon(c2), lat(c2))
position_angle(c1::T1, c2::T2) where {T1 <: AbstractSkyCoords, T2 <: AbstractSkyCoords} = position_angle(c1, convert(T1, c2))


function _position_angle(λ1, ϕ1, λ2, ϕ2)
    sin_Δλ, cos_Δλ = sincos(λ2 - λ1)
    sin_ϕ1, cos_ϕ1 = sincos(ϕ1)
    sin_ϕ2, cos_ϕ2 = sincos(ϕ2)

    x = sin_ϕ2 * cos_ϕ1 - cos_ϕ2 * sin_ϕ1 * cos_Δλ
    y = sin_Δλ * cos_ϕ2

    return mod2pi(atan(y, x))
end

"""
    offset(::AbstractSkyCoords, separation, pa) -> coordinate

Offset a coordinate by a given angular separation, `separation`, in radians and position angle, `pa`, in radians.

Uses the sine and cosine rules in spherical coordinates with corrections for the antipodes.
Returns a sky coordinate of the same type as input.

### Examples
```jldoctest
julia> c1 = ICRSCoords(0, 0);

julia> c2 = offset(c1, deg2rad(1), deg2rad(90))
ICRSCoords{Float64}(0.017453292519943295, 1.0686516840418957e-18)

julia> offset(c1, c2) .|> rad2deg
(1.0, 90.0)
```

### See Also
[`separation`](@ref), [`position_angle`](@ref)
"""
offset(c::T, sep, pa) where {T <: AbstractSkyCoords} = T(_offset(lon(c), lat(c), sep, pa)...)

"""
    offset(::AbstractSkyCoords, AbstractSkyCoords) -> angle, angle

Return the separation and position angle in radians between two sky coordinates.

### Examples
```jldoctest
julia> c1 = ICRSCoords(0, 0); c2 = ICRSCoords(deg2rad(1), 0);

julia> offset(c1, c2) .|> rad2deg
(1.0, 90.0)
```

### See Also
[`separation`](@ref), [`position_angle`](@ref)
"""
offset(c1::AbstractSkyCoords, c2::AbstractSkyCoords) = separation(c1, c2), position_angle(c1, c2)

#= use the cosine rule in spherical geometry with three points, the north pole, the starting point,
and the final point.
angles: (change in lon), (position angle), (-1/position angle)
sides: (separation), (final co-latitude), (starting co-latitude)
=#
function _offset(λ, ϕ, separation, pa)
    sin_a, cos_a = sincos(separation)
    cos_c, sin_c = sincos(ϕ)
    sin_B, cos_B = sincos(pa)

    # solving cosine rule
    cos_b = cos_c * cos_a + sin_c * sin_a * cos_B

    # solving sine rule
    xsin_A = sin_a * sin_B * sin_c
    xcos_A = cos_a - cos_b * cos_c

    # correction for antipodes, otherwise atan2
    ang = sin_c < 1.0e-12 ? π / 2 + cos_c * (π / 2 - pa) : atan(xsin_A, xcos_A)

    return mod2pi(λ + ang), asin(cos_b)
end

# Stub to extend in NearestNeighbors extension
function match end

end # module
