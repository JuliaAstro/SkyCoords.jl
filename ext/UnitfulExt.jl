module UnitfulExt

using Unitful
using SkyCoords

_COORDTYPES_LATLON = Union{ICRSCoords, GalCoords, FK5Coords, EclipticCoords}

(::Type{T})(lon::Quantity, lat) where {T<:_COORDTYPES_LATLON} = T(ustrip(u"rad", lon), lat)
(::Type{T})(lon, lat::Quantity) where {T<:_COORDTYPES_LATLON} = T(lon, ustrip(u"rad", lat))
(::Type{T})(lon::Quantity, lat::Quantity) where {T<:_COORDTYPES_LATLON} = T(ustrip(u"rad", lon), ustrip(u"rad", lat))

SkyCoords.lon(u::Unitful.Units, c) = SkyCoords.lon(c) * u"rad" |> u
SkyCoords.lat(u::Unitful.Units, c) = SkyCoords.lat(c) * u"rad" |> u
SkyCoords.lonlat(u::Unitful.Units, c) = SkyCoords.lonlat(c) .* u"rad" .|> u

SkyCoords.separation(u::Unitful.Units, c1, c2) = SkyCoords.separation(c1, c2) * u"rad" |> u
SkyCoords.position_angle(u::Unitful.Units, c1, c2) = SkyCoords.position_angle(c1, c2) * u"rad" |> u

end
