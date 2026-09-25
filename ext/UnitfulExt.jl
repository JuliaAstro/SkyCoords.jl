module UnitfulExt

using Unitful
using SkyCoords

_COORDTYPES_LATLON = Union{ICRSCoords, GalCoords, SuperGalCoords, FK5Coords, EclipticCoords}

# Every slot is typed (`Union{Real, Quantity}` rather than `Any`) so these
# methods stay disjoint from the Real-typed coordinate constructors and from
# the DynamicQuantities extension, keeping dispatch free of ambiguities.
(::Type{T})(lon::Quantity, lat::Union{Real, Quantity}) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), lat)
(::Type{T})(lon::Union{Real, Quantity}, lat::Quantity) where {T <: _COORDTYPES_LATLON} = T(lon, ustrip(u"rad", lat))
(::Type{T})(lon::Quantity, lat::Quantity) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), ustrip(u"rad", lat))

SkyCoords.lon(u::Unitful.Units, c) = SkyCoords.lon(c) * u"rad" |> u
SkyCoords.lat(u::Unitful.Units, c) = SkyCoords.lat(c) * u"rad" |> u
SkyCoords.lonlat(u::Unitful.Units, c) = SkyCoords.lonlat(c) .* u"rad" .|> u

SkyCoords.separation(u::Unitful.Units, c1, c2) = SkyCoords.separation(c1, c2) * u"rad" |> u
SkyCoords.position_angle(u::Unitful.Units, c1, c2) = SkyCoords.position_angle(c1, c2) * u"rad" |> u

end
