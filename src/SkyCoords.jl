module SkyCoords
using Compat

export AbstractSkyCoords,
       ICRSCoords,
       GalCoords,
       FK5Coords,
       separation

import Base: convert, *, transpose

# -----------------------------------------------------------------------------
# Types

abstract AbstractSkyCoords

immutable ICRSCoords{T<:AbstractFloat} <: AbstractSkyCoords
    ra::T
    dec::T
    ICRSCoords(ra, dec) = new(mod(ra, 2pi), dec)
end
ICRSCoords{T<:AbstractFloat}(ra::T, dec::T) = ICRSCoords{T}(ra, dec)
ICRSCoords(ra::Real, dec::Real) = ICRSCoords(promote(float(ra), float(dec))...)

immutable GalCoords{T<:AbstractFloat} <: AbstractSkyCoords
    l::T
    b::T
    GalCoords(l, b) = new(mod(l, 2pi), b)
end
GalCoords{T<:AbstractFloat}(l::T, b::T) = GalCoords{T}(l, b)
GalCoords(l::Real, b::Real) = GalCoords(promote(float(l), float(b))...)

# FK5 is parameterized by equinox (e)
immutable FK5Coords{e, T<:AbstractFloat} <: AbstractSkyCoords
    ra::T
    dec::T
    FK5Coords(ra, dec) = new(mod(ra, 2pi), dec)
end
@compat (::Type{FK5Coords{e}}){e,T}(ra::T, dec::T) = FK5Coords{e, T}(ra, dec)

# We'd like to define this promotion constructor, but in Julia 0.5,
# the typing algorithm can't figure out that the previous method is
# more specific, so this promotion constructor calls itself, resulting in
# stack overflow.
#(::Type{FK5Coords{e}}){e}(ra::Real, dec::Real) =
#    FK5Coords{e}(promote(float(ra), float(dec))...)

# -----------------------------------------------------------------------------
# Helper functions: Immutable array operations

# We use immutable array operations to avoid allocating memory for
# small arrays. Eventually, base Julia should support immutable arrays, at
# which point we should switch to using that functionality in base.

immutable Matrix33{T<:AbstractFloat}
    a11::T
    a12::T
    a13::T
    a21::T
    a22::T
    a23::T
    a31::T
    a32::T
    a33::T
end

@compat (::Type{Matrix33{T}}){T}(m::Matrix33{T}) = m
@compat (::Type{Matrix33{T}}){T}(m::Matrix33) =
    Matrix33(T(m.a11), T(m.a12), T(m.a13),
             T(m.a21), T(m.a22), T(m.a23),
             T(m.a31), T(m.a32), T(m.a33))

immutable Vector3{T<:AbstractFloat}
    x::T
    y::T
    z::T
end

function *(x::Matrix33, y::Matrix33)
    Matrix33(x.a11 * y.a11 + x.a12 * y.a21 + x.a13 * y.a31,
             x.a11 * y.a12 + x.a12 * y.a22 + x.a13 * y.a32,
             x.a11 * y.a13 + x.a12 * y.a23 + x.a13 * y.a33,
             x.a21 * y.a11 + x.a22 * y.a21 + x.a23 * y.a31,
             x.a21 * y.a12 + x.a22 * y.a22 + x.a23 * y.a32,
             x.a21 * y.a13 + x.a22 * y.a23 + x.a23 * y.a33,
             x.a31 * y.a11 + x.a32 * y.a21 + x.a33 * y.a31,
             x.a31 * y.a12 + x.a32 * y.a22 + x.a33 * y.a32,
             x.a31 * y.a13 + x.a32 * y.a23 + x.a33 * y.a33)
end

transpose(m::Matrix33) = Matrix33(m.a11, m.a21, m.a31,
                                  m.a12, m.a22, m.a32,
                                  m.a13, m.a23, m.a33)

function *(m::Matrix33, v::Vector3)
    Vector3(m.a11 * v.x + m.a12 * v.y + m.a13 * v.z,
            m.a21 * v.x + m.a22 * v.y + m.a23 * v.z,
            m.a31 * v.x + m.a32 * v.y + m.a33 * v.z)
end


# -----------------------------------------------------------------------------
# Helper functions: Create rotation matrix about a given axis (x, y, z)

function xrotmat(angle)
    s = sin(angle)
    c = cos(angle)
    Matrix33(1., 0., 0.,
             0.,  c,  s,
             0., -s,  c)
end

function yrotmat(angle)
    s = sin(angle)
    c = cos(angle)
    Matrix33(c,  0., -s,
             0., 1., 0.,
             s,  0.,  c)
end

function zrotmat(angle)
    s = sin(angle)
    c = cos(angle)
    Matrix33(c,   s,  0.,
             -s,  c,  0.,
             0., 0.,  1.)
end

# (lon, lat) -> [x, y, z] unit vector
function coords2cart(lon, lat)
    coslat = cos(lat)
    Vector3(coslat*cos(lon), coslat*sin(lon), sin(lat))
end

# [x, y, z] unit vector -> (lon, lat)
cart2coords(v) = atan2(v.y, v.x), atan2(v.z, sqrt(v.x*v.x + v.y*v.y))

# -----------------------------------------------------------------------------
# Constant rotation matricies and precession matrix function

# ICRS --> FK5 at J2000 (See USNO Circular 179, section 3.5)
eta0 = deg2rad(-19.9 / 3600000.)
xi0 = deg2rad(9.1 / 3600000.)
da0 = deg2rad(-22.9 / 3600000.)
const ICRS_TO_FK5J2000 = xrotmat(-eta0) * yrotmat(xi0) * zrotmat(da0)
const FK5J2000_TO_ICRS = ICRS_TO_FK5J2000'

# FK5J2000 --> Gal
# Note that galactic pole and zeropoint of l are somewhat arbitrary
# and not officially defined (as far as I know). The values below are
# from astropy.coordinates, which includes the following comment:
# | "This gives better consistency with other codes than using the values
# |  from Reid & Brunthaler 2004 and the best self-consistency between FK5
# |  -> Galactic and FK5 -> FK4 -> Galactic. The lon0 angle was found by
# |  optimizing the self-consistency."
ngp_fk5j2000_ra = deg2rad(192.8594812065348)
ngp_fk5j2000_dec = deg2rad(27.12825118085622)
lon0_fk5j2000 = deg2rad(122.9319185680026)
const FK5J2000_TO_GAL = (zrotmat(pi - lon0_fk5j2000) *
                         yrotmat(pi/2. - ngp_fk5j2000_dec) *
                         zrotmat(ngp_fk5j2000_ra))
const GAL_TO_FK5J2000 = FK5J2000_TO_GAL'

# Gal --> ICRS: simply chain through FK5J2000
const GAL_TO_ICRS = FK5J2000_TO_ICRS * GAL_TO_FK5J2000
const ICRS_TO_GAL = GAL_TO_ICRS'

# FK5J2000 --> FK5{epoch}
# Computes the precession matrix from J2000 to the given Julian equinox.
# Expression from from Capitaine et al. 2003 as expressed in the USNO
# Circular 179.  This should match the IAU 2006 standard from SOFA.
const pzeta = [2.650545, 2306.083227, 0.2988499, 0.01801828, -0.000005971,
               -0.0000003173]
const pz = [-2.650545, 2306.077181, 1.0927348, 0.01826837, -0.000028596,
            -0.0000002904]
const ptheta = [0.0, 2004.191903, -0.4294934, -0.04182264, -0.000007089,
                -0.0000001274]
function precess_from_j2000(equinox)
    t = (equinox - 2000.0) / 100.0
    tn = 1.0
    zeta = pzeta[1]
    z = pz[1]
    theta = ptheta[1]
    for i = 2:6
        tn *= t
        zeta += pzeta[i] * tn
        z += pz[i] * tn
        theta += ptheta[i] * tn
    end
    zeta = deg2rad(zeta/3600.0)
    z = deg2rad(z/3600.0)
    theta = deg2rad(theta/3600.0)
    zrotmat(-z) * yrotmat(theta) * zrotmat(-zeta)
end

# -----------------------------------------------------------------------------
# Type-dependent methods

lon(c::GalCoords) = c.l
lat(c::GalCoords) = c.b
lon(c::AbstractSkyCoords) = c.ra
lat(c::AbstractSkyCoords) = c.dec

# Abstract away specific field names (ra, dec vs l, b)
coords2cart(c::AbstractSkyCoords) = coords2cart(lon(c), lat(c))

# Rotation matrix between coordinate systems: `rotmat(to, from)`
# Note that all of these return Matrix33{Float64}, regardless of
# element type of input coordinates.
rotmat{T1<:GalCoords, T2<:ICRSCoords}(::Type{T1}, ::Type{T2}) = ICRS_TO_GAL
rotmat{T1<:ICRSCoords, T2<:GalCoords}(::Type{T1}, ::Type{T2}) = GAL_TO_ICRS

# Define both these so that `convert(FK5Coords{e}, ...)` and
# `convert(FK5Coords{e,T}, ...)` both work. Similar with other
# FK5Coords rotmat methods below.
rotmat{e1, T2<:ICRSCoords}(::Type{FK5Coords{e1}}, ::Type{T2}) =
    precess_from_j2000(e1) * ICRS_TO_FK5J2000
rotmat{e1, T1, T2<:ICRSCoords}(::Type{FK5Coords{e1,T1}}, ::Type{T2}) =
    precess_from_j2000(e1) * ICRS_TO_FK5J2000

rotmat{e1, T2<:GalCoords}(::Type{FK5Coords{e1}}, ::Type{T2}) =
    precess_from_j2000(e1) * GAL_TO_FK5J2000
rotmat{e1, T1, T2<:GalCoords}(::Type{FK5Coords{e1,T1}}, ::Type{T2}) =
    precess_from_j2000(e1) * GAL_TO_FK5J2000

rotmat{T1<:ICRSCoords, e2, T2}(::Type{T1}, ::Type{FK5Coords{e2,T2}}) =
    FK5J2000_TO_ICRS * precess_from_j2000(e2)'

rotmat{T1<:GalCoords, e2, T2}(::Type{T1}, ::Type{FK5Coords{e2,T2}}) =
    FK5J2000_TO_GAL * precess_from_j2000(e2)'

rotmat{e1, e2, T2}(::Type{FK5Coords{e1}}, ::Type{FK5Coords{e2,T2}}) =
    precess_from_j2000(e1) * precess_from_j2000(e2)'
rotmat{e1, T1, e2, T2}(::Type{FK5Coords{e1,T1}}, ::Type{FK5Coords{e2,T2}}) =
    precess_from_j2000(e1) * precess_from_j2000(e2)'

# get floating point type in coordinates
_eltype{e,T}(::Type{FK5Coords{e,T}}) = T
_eltype{T}(::Type{GalCoords{T}}) = T
_eltype{T}(::Type{ICRSCoords{T}}) = T
_eltype(c::AbstractSkyCoords) = _eltype(typeof(c))

# Scalar coordinate conversions
convert{T<:AbstractSkyCoords}(::Type{T}, c::T) = c
function convert{T<:AbstractSkyCoords, S<:AbstractSkyCoords}(::Type{T}, c::S)
    r = Matrix33{_eltype(c)}(rotmat(T, S)) * coords2cart(c)
    lon, lat = cart2coords(r)
    T(lon, lat)
end

# Vector coordinate conversions
# This is useful for FK5Coords so that we only compute the rotation
# matrix (a function of epoch) once and apply it to the entire
# vector. The compiler doesn't seem to figure out that rotmat(T, S) is
# a constant for FK5Coords{e} types.
convert{T<:AbstractSkyCoords,n}(::Type{Array{T,n}}, c::Array{T,n}) = c
function convert{T<:AbstractSkyCoords, n, S<:AbstractSkyCoords}(
    ::Type{Array{T,n}}, c::Array{S, n})
    m = Matrix33{_eltype(S)}(rotmat(T, S))
    result = similar(c, T)
    for i in 1:length(c)
        r = m * coords2cart(c[i])
        lon, lat = cart2coords(r)
        result[i] = T(lon, lat)
    end
    result
end

# ------------------------------------------------------------------------------
# Distance between coordinates

function _separation(λ_1, ϕ_1, λ_2, ϕ_2)
    Δλ = λ_2 - λ_1
    sin_Δλ = sin(Δλ)
    cos_Δλ = cos(Δλ)
    sin_ϕ1 = sin(ϕ_1)
    sin_ϕ2 = sin(ϕ_2)
    cos_ϕ1 = cos(ϕ_1)
    cos_ϕ2 = cos(ϕ_2)
    return atan2(hypot(cos_ϕ2 * sin_Δλ,
                       cos_ϕ1 * sin_ϕ2 - sin_ϕ1 * cos_ϕ2 * cos_Δλ),
                 sin_ϕ1 * sin_ϕ2 + cos_ϕ1 * cos_ϕ2 * cos_Δλ)
end

"""
    separation(c1::AbstractSkyCoords, c2::AbstractSkyCoords) -> distance

Return angular separation between two sky coordinates, in radians.

The angular separation is calculated using the Vincenty formula
(http://en.wikipedia.org/wiki/Great-circle_distance), which is slightly more
complex and computationally expensive than some alternatives, but is stable at
at all distances, including the poles and antipodes.
"""
separation{T<:AbstractSkyCoords}(c1::T, c2::T) =
    _separation(lon(c1), lat(c1), lon(c2), lat(c2))

separation{T1<:AbstractSkyCoords,T2<:AbstractSkyCoords}(c1::T1, c2::T2) =
    separation(c1, convert(T1, c2))

function separation{T<:AbstractSkyCoords}(c1::AbstractArray{T},
                                          c2::AbstractArray{T})
    @assert size(c1) == size(c2) "Size mismatch"
    result = similar(c1, _eltype(first(c1)))
    for i in eachindex(c1)
        result[i] = separation(c1[i], c2[i])
    end
    return result
end


end # module
