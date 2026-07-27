module DynamicQuantitiesExt

using DynamicQuantities
using SkyCoords

# The constructors below strip units positionally, so AltAzCoords fits in even
# though its two angle arguments are (alt, az) rather than (lon, lat)
_COORDTYPES_LATLON = Union{ICRSCoords, GalCoords, SuperGalCoords, FK4Coords, FK4NoETerms, FK5Coords, EclipticCoords, AltAzCoords}

# `AbstractRealQuantity` is deliberately excluded from the quantity slots: it
# subtypes `Real`, so accepting it here would make these methods ambiguous
# against the Real-typed coordinate constructors. `RealQuantity` angles are
# not supported for construction; use `Quantity` or plain radians.
const _QUANTITY_NOT_REAL = Union{AbstractQuantity, AbstractGenericQuantity}

# Every slot is typed (rather than `Any`) so these methods stay disjoint from
# the Real-typed coordinate constructors and from the Unitful extension,
# keeping dispatch free of ambiguities.
(::Type{T})(lon::_QUANTITY_NOT_REAL, lat::Union{Real, UnionAbstractQuantity}) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), lat)
(::Type{T})(lon::Union{Real, UnionAbstractQuantity}, lat::_QUANTITY_NOT_REAL) where {T <: _COORDTYPES_LATLON} = T(lon, ustrip(u"rad", lat))
(::Type{T})(lon::_QUANTITY_NOT_REAL, lat::_QUANTITY_NOT_REAL) where {T <: _COORDTYPES_LATLON} = T(ustrip(u"rad", lon), ustrip(u"rad", lat))

# `Observer` takes its latitude/longitude as angles and its altitude as a length.
# Quantities and plain numbers can be mixed; all-plain calls dispatch to the
# base constructor, so the plain slots here just pass through unchanged.
_strip(u, x) = x isa UnionAbstractQuantity ? ustrip(u, x) : x
SkyCoords.Observer(latitude::Union{Real, UnionAbstractQuantity}, longitude::Union{Real, UnionAbstractQuantity}, altitude::Union{Real, UnionAbstractQuantity} = 0) =
    Observer(_strip(u"rad", latitude), _strip(u"rad", longitude), _strip(u"m", altitude))

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
