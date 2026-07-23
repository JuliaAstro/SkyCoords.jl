module DynamicQuantitiesExt

using DynamicQuantities
using SkyCoords

# The constructors below strip units positionally, so AltAzCoords fits in even
# though its two angle arguments are (alt, az) rather than (lon, lat)
_COORDTYPES_LATLON = Union{ICRSCoords, GalCoords, FK4Coords, FK4NoETerms, FK5Coords, EclipticCoords, AltAzCoords}

(::Type{T})(lon::UnionAbstractQuantity, lat) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), lat)
(::Type{T})(lon, lat::UnionAbstractQuantity) where {T <: _COORDTYPES_LATLON} = T(lon, ustrip(u"rad", lat))
(::Type{T})(lon::UnionAbstractQuantity, lat::UnionAbstractQuantity) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), ustrip(u"rad", lat))

SkyCoords.lon(u::UnionAbstractQuantity, c) = SkyCoords.lon(c) * u"rad" |> u
SkyCoords.lat(u::UnionAbstractQuantity, c) = SkyCoords.lat(c) * u"rad" |> u
SkyCoords.lonlat(u::UnionAbstractQuantity, c) = SkyCoords.lonlat(c) .* u"rad" .|> u

SkyCoords.separation(u::UnionAbstractQuantity, c1, c2) = SkyCoords.separation(c1, c2) * u"rad" |> u
SkyCoords.position_angle(u::UnionAbstractQuantity, c1, c2) = SkyCoords.position_angle(c1, c2) * u"rad" |> u

# `offset` does trigonometry on `sep`/`pa`. DynamicQuantities rejects trig on
# angular quantities (symbolic radians are not dimensionless), so strip the
# angular arguments to plain radians first.
SkyCoords.offset(c::AbstractSkyCoords, sep::UnionAbstractQuantity, pa) = SkyCoords.offset(c, ustrip(u"rad", sep), pa)
SkyCoords.offset(c::AbstractSkyCoords, sep, pa::UnionAbstractQuantity) = SkyCoords.offset(c, sep, ustrip(u"rad", pa))
SkyCoords.offset(c::AbstractSkyCoords, sep::UnionAbstractQuantity, pa::UnionAbstractQuantity) = SkyCoords.offset(c, ustrip(u"rad", sep), ustrip(u"rad", pa))

end
