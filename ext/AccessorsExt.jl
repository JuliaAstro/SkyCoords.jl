module AccessorsExt
using Accessors
import Accessors: set
using SkyCoords
using SkyCoords: lat, lon

set(x::GalCoords, ::typeof(lon), v) = @set x.l = v
set(x::GalCoords, ::typeof(lat), v) = @set x.b = v
set(x::EclipticCoords, ::typeof(lon), v) = @set x.lon = v
set(x::EclipticCoords, ::typeof(lat), v) = @set x.lat = v
set(x::ICRSCoords, ::typeof(lon), v) = @set x.ra = v
set(x::ICRSCoords, ::typeof(lat), v) = @set x.dec = v
set(x::FK5Coords, ::typeof(lon), v) = @set x.ra = v
set(x::FK5Coords, ::typeof(lat), v) = @set x.dec = v

set(x::CartesianCoords, ::typeof(vec), v) = @set x.vec = v

set(x::CartesianCoords, ::typeof(cartesian), v::CartesianCoords) = cartesian(v)
set(x::CartesianCoords, ::typeof(spherical), v::AbstractSkyCoords) = cartesian(v)
set(x::AbstractSkyCoords, ::typeof(cartesian), v::CartesianCoords) = spherical(v)
set(x::AbstractSkyCoords, ::typeof(spherical), v::AbstractSkyCoords) = spherical(v)

end
