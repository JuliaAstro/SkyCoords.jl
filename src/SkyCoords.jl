module SkyCoords

import ConstructionBase: constructorof
using LinearAlgebra: I, norm
using Rotations
using StaticArrays

export AbstractSkyCoords, 
       ICRSCoords,
       GalCoords,
       SuperGalCoords,
       FK5Coords,
       EclipticCoords,
       CartesianCoords,
       separation,
       position_angle,
       offset,
       cartesian,
       spherical,
       match_coords,
       CoordsKDTree

include("types.jl")
include("cartesian.jl")

# -----------------------------------------------------------------------------
# Helper functions: Create rotation matrix about a given axis (x, y, z)

# (lon, lat) -> [x, y, z] unit vector
function coords2cart(lon, lat)
    sinlon, coslon = sincos(lon)
    sinlat, coslat = sincos(lat)
    SVector{3}(coslat * coslon, coslat * sinlon, sinlat)
end

# [x, y, z] unit vector -> (lon, lat)
function cart2coords(v)
    lon = atan(v[begin + 1], v[begin])
    xy_norm = hypot(v[begin], v[begin + 1])
    lat = atan(v[begin + 2], xy_norm)
    return lon, lat
end

# -----------------------------------------------------------------------------
# Constant rotation matricies and precession matrix function

# ICRS --> FK5 at J2000 (See USNO Circular 179, section 3.5)
const ICRS_TO_FK5J2000 = let
    eta0 = deg2rad(-19.9 / 3600000)
    xi0 = deg2rad(9.1 / 3600000)
    da0 = deg2rad(-22.9 / 3600000)
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
    RotZYZ(lon0_fk5j2000 - π, ngp_fk5j2000_dec - π/2, -ngp_fk5j2000_ra)
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
const GAL_TO_SUPERGAL = RotZYZ(π/2, π/2 - sgp_b, π - sgp_l)
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
    obl = @evalpoly(T, 84381.406, -46.836769, -0.0001831, 0.00200340, -0.000000576, -0.0000000434)
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
@generated rotmat(::Type{<:EclipticCoords{e}}, ::Type{T}) where {e,T<:AbstractSkyCoords} = rotmat(EclipticCoords{e}, FK5Coords{e}) * rotmat(FK5Coords{e}, T)
@generated rotmat(::Type{T}, ::Type{<:EclipticCoords{e}}) where {e,T<:AbstractSkyCoords} = rotmat(T, FK5Coords{e}) * rotmat(FK5Coords{e}, EclipticCoords{e})
# disambiguation:
@generated rotmat(::Type{<:EclipticCoords{e_to}}, ::Type{<:EclipticCoords{e_from}}) where {e_to,e_from} = rotmat(EclipticCoords{e_to}, FK5Coords{e_to}) * rotmat(FK5Coords{e_to}, EclipticCoords{e_from})

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
@generated rotmat(::Type{<:FK5Coords{e1}}, ::Type{<:FK5Coords{e2}}) where {e1,e2} =
    precess_from_j2000(e1) * precess_from_j2000(e2)'

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

The angular separation is calculated using the [Vincenty formula](http://en.wikipedia.org/wiki/Great-circle_distance), which is slightly more
complex and computationally expensive than some alternatives, but is stable at
at all distances, including the poles and antipodes.
"""
separation(c1::T, c2::T) where {T<:AbstractSkyCoords} =
    _separation(lon(c1), lat(c1), lon(c2), lat(c2))

separation(c1::CartesianCoords{T}, c2::CartesianCoords{T}) where {T<:AbstractSkyCoords} =
    2 * asin(norm(vec(c1) - vec(c2)) / 2)

separation(c1::T1, c2::T2) where {T1<:AbstractSkyCoords,T2<:AbstractSkyCoords} =
    separation(c1, convert(T1, c2))


"""
    position_angle(c1::AbstractSkyCoords, c2::AbstractSkyCoords) -> angle

Return position angle between two sky coordinates, in positive radians.

# Examples
```jldoctest
julia> c1 = ICRSCoords(0, 0); c2 = ICRSCoords(deg2rad(1), 0);

julia> position_angle(c1, c2) |> rad2deg
90.0
```
"""
position_angle(c1::T, c2::T) where {T <: AbstractSkyCoords} = _position_angle(lon(c1), lat(c1), lon(c2), lat(c2))
position_angle(c1::T1, c2::T2) where {T1 <: AbstractSkyCoords,T2 <: AbstractSkyCoords} = position_angle(c1, convert(T1, c2))


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

Uses the sine and cosine rules in spherical coordinates with corrections for the antipodes. Returns a sky coordinate of the same type as input.

# Examples
```jldoctest
julia> c1 = ICRSCoords(0, 0);

julia> c2 = offset(c1, deg2rad(1), deg2rad(90))
ICRSCoords{Float64}(0.017453292519943295, 1.0686516840418957e-18)

julia> offset(c1, c2) .|> rad2deg
(1.0, 90.0)
```

# See Also
* [`separation`](@ref), [`position_angle`](@ref)
"""
offset(c::T, sep, pa) where T <: AbstractSkyCoords = T(_offset(lon(c), lat(c), sep, pa)...)

"""
    offset(::AbstractSkyCoords, AbstractSkyCoords) -> angle, angle

Return the separation and position angle in radians between two sky coordinates.

# Examples
```jldoctest
julia> c1 = ICRSCoords(0, 0); c2 = ICRSCoords(deg2rad(1), 0);

julia> offset(c1, c2) .|> rad2deg
(1.0, 90.0)
```

# See Also
* [`separation`](@ref), [`position_angle`](@ref)
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
    ang = sin_c < 1e-12 ? π/2 + cos_c * (π/2 - pa) : atan(xsin_A, xcos_A)

    return mod2pi(λ + ang), asin(cos_b)
end

"""
    match_coords(refcoords::AbstractArray{<:AbstractSkyCoords}, 
                 matchcoords::AbstractArray{<:AbstractSkyCoords};
                 nthneighbor::Int = 1)
*Requires Julia ≥ 1.9 and NearestNeighbors.jl to be loaded (e.g., `using NearestNeighbors`).*

Finds the nearest entries in `refcoords` to the coordinates contained in `matchcoords`. The keyword argument `nthneighbor` determines *which* nearest neighbor to search for; typically this should be `1` when matching one set of coordinates to another. Another common use case is setting `nthneighbor = 2` when matching a catalog against itself to find the nearest neighbor of each coordinate in the same catalog. 

Returns `(id, sep)`, where
 - `id` is an array containing indices of the coordinates in `refcoords` that matched with the elements of `matchcoords`, and 
 - `sep` is an array giving the angular separation between the elements of `matchcoords` and the above matches.

Note that this method creates a [`CoordsKDTree`](@ref) from `refcoords` and then calls the method below. If you plan to use the same `refcoords` to match to many different `matchcoords`, then you should directly construct `CoordsKDTree(refcoords)` and call the method below.

    match_coords(tree::CoordsKDTree, matchcoords::AbstractArray{<:AbstractSkyCoords};
                 nthneighbor::Int = 1)
As above, but uses a pre-constructed `tree::CoordsKDTree` rather than creating one from a reference catalog of coordinates.
"""
function match_coords end

"""
    CoordsKDTree(data::AbstractArray{<:AbstractSkyCoords}; kws...)
*Requires Julia ≥ 1.9 and NearestNeighbors.jl to be loaded (e.g., `using NearestNeighbors`).*

A wrapper for a NearestNeighbors.jl `KDTree` that also includes the type of coordinate that it was constructed from as the parameter `TC <: AbstractSkyCoords`. The provided `data` are used to construct the tree, and the `kws...` are passed to `NearestNeighbors.KDTree`. The only field is `tree`, which provides access to the underlying `KDTree`. An internal Cartesian coordinate representation is used, so nearest neighbor queries on the `KDTree` must first be converted into the coordinate `TC`, then to Cartesian coordinates. This conversion is applied automatically in the following methods, which extend those from `NearestNeighbors.jl`.
 - `nn(tree::CoordsKDTree, coord::AbstractSkyCoords)` queries the `tree` for the entry nearest the provided `coord`.
 - `nn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords})` queries the `tree` for the entry nearest each coordinate in `coords`.
 - `knn(tree::CoordsKDTree, coord::AbstractSkyCoords, k::Int)` queries the `tree` for the `k` coordinates nearest the provided `coord`.
 - `knn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords}, k::Int)` queries the `tree` for the `k` coordinates nearest each coordinate in `coords`.
"""
struct CoordsKDTree{TC <: AbstractSkyCoords, K}
    tree::K
end

end # module
