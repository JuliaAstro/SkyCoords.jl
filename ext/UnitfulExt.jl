module UnitfulExt

using Unitful
using SkyCoords

# The constructors below strip units positionally, so AltAzCoords fits in even
# though its two angle arguments are (alt, az) rather than (lon, lat)
_COORDTYPES_LATLON = Union{ICRSCoords, GalCoords, SuperGalCoords, FK4Coords, FK4NoETerms, FK5Coords, EclipticCoords, AltAzCoords}

# Every slot is typed (`Union{Real, Quantity}` rather than `Any`) so these
# methods stay disjoint from the Real-typed coordinate constructors and from
# the DynamicQuantities extension, keeping dispatch free of ambiguities.
(::Type{T})(lon::Quantity, lat::Union{Real, Quantity}) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), lat)
(::Type{T})(lon::Union{Real, Quantity}, lat::Quantity) where {T <: _COORDTYPES_LATLON} = T(lon, ustrip(u"rad", lat))
(::Type{T})(lon::Quantity, lat::Quantity) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), ustrip(u"rad", lat))

# `Observer` takes its latitude/longitude as angles and its altitude as a length.
# Quantities and plain numbers can be mixed; all-plain calls dispatch to the
# base constructor, so the plain slots here just pass through unchanged.
_strip(u, x) = x isa Quantity ? ustrip(u, x) : x
SkyCoords.Observer(latitude::Union{Real, Quantity}, longitude::Union{Real, Quantity}, altitude::Union{Real, Quantity} = 0) =
    Observer(_strip(u"rad", latitude), _strip(u"rad", longitude), _strip(u"m", altitude))

SkyCoords.lon(u::Unitful.Units, c) = SkyCoords.lon(c) * u"rad" |> u
SkyCoords.lat(u::Unitful.Units, c) = SkyCoords.lat(c) * u"rad" |> u
SkyCoords.lonlat(u::Unitful.Units, c) = SkyCoords.lonlat(c) .* u"rad" .|> u

SkyCoords.separation(u::Unitful.Units, c1, c2) = SkyCoords.separation(c1, c2) * u"rad" |> u
SkyCoords.position_angle(u::Unitful.Units, c1, c2) = SkyCoords.position_angle(c1, c2) * u"rad" |> u

end
