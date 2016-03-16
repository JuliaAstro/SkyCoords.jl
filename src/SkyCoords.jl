module SkyCoords
using Compat

export AbstractSkyCoords,
       ICRSCoords,
       GalCoords,
       FK5Coords

import Base: convert, *, transpose

# -----------------------------------------------------------------------------
# Types

abstract AbstractSkyCoords

immutable ICRSCoords <: AbstractSkyCoords
    ra::Float64
    dec::Float64
    ICRSCoords(ra::Real, dec::Real) = @compat new(mod(Float64(ra), 2pi),
                                                  Float64(dec))
end

immutable GalCoords <: AbstractSkyCoords
    l::Float64
    b::Float64
    GalCoords(l::Real, b::Real) = @compat new(mod(Float64(l), 2pi),
                                              Float64(b))
end

# FK5 is parameterized by equinox (e)
immutable FK5Coords{e} <: AbstractSkyCoords
    ra::Float64
    dec::Float64
    FK5Coords(ra::Real, dec::Real) = @compat new(mod(Float64(ra), 2pi),
                                                 Float64(dec))
end

# -----------------------------------------------------------------------------
# Helper functions: Immutable array operations

# We use immutable array operations to avoid allocating memory for
# small arrays. Eventually, base Julia should support immutable arrays, at
# which point we should switch to using that functionality in base.

immutable Matrix33
    a11::Float64
    a12::Float64
    a13::Float64
    a21::Float64
    a22::Float64
    a23::Float64
    a31::Float64
    a32::Float64
    a33::Float64
end

immutable Vector3
    x::Float64
    y::Float64
    z::Float64
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

# Abstract away specific field names (ra, dec vs l, b)
coords2cart(c::ICRSCoords) = coords2cart(c.ra, c.dec)
coords2cart(c::GalCoords) = coords2cart(c.l, c.b)
coords2cart{e}(c::FK5Coords{e}) = coords2cart(c.ra, c.dec)

# Rotation matrix between coordinate systems: `rotmat(to, from)`
rotmat(::Type{GalCoords}, ::Type{ICRSCoords}) = ICRS_TO_GAL
rotmat(::Type{ICRSCoords}, ::Type{GalCoords}) = GAL_TO_ICRS
rotmat{e}(::Type{FK5Coords{e}}, ::Type{ICRSCoords}) =
    precess_from_j2000(e) * ICRS_TO_FK5J2000
rotmat{e}(::Type{FK5Coords{e}}, ::Type{GalCoords}) =
   precess_from_j2000(e) * GAL_TO_FK5J2000
rotmat{e}(::Type{ICRSCoords}, ::Type{FK5Coords{e}}) = 
    FK5J2000_TO_ICRS * precess_from_j2000(e)'
rotmat{e}(::Type{GalCoords}, ::Type{FK5Coords{e}}) =
    FK5J2000_TO_GAL * precess_from_j2000(e)'
rotmat{e1, e2}(::Type{FK5Coords{e1}}, ::Type{FK5Coords{e2}}) =
    precess_from_j2000(e1) * precess_from_j2000(e2)'

# Scalar coordinate conversions
convert{T<:AbstractSkyCoords}(::Type{T}, c::T) = c
function convert{T<:AbstractSkyCoords, S<:AbstractSkyCoords}(::Type{T}, c::S)
    r = rotmat(T, S) * coords2cart(c)
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
    m = rotmat(T, S)
    result = similar(c, T)
    for i in 1:length(c)
        r = m * coords2cart(c[i])
        lon, lat = cart2coords(r)
        result[i] = T(lon, lat)
    end
    result
end

end # module
