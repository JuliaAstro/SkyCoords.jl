abstract type AbstractProjectedCoords <: AbstractSkyCoords end

"""
    ProjectedCoords{TC <: AbstractSkyCoords, T <: Real} <: AbstractSkyCoords

Represent coordinate `c` as a flat-sky offset from `origin` using a
small-field-of-view approximation.

`offset[1]` is the longitude difference scaled by `cos(lat(origin))`, and
`offset[2]` is the latitude difference, both in radians.

A projected coordinate is a _representation_ around a point, not a coordinate frame of its own.
Its frame is the origin's. Conversions out of a `ProjectedCoords` (spherical or [`CartesianCoords`](@ref) targets alike)
work like any other coordinate. Conversions into one are not possible via `convert`, since they require an origin value.
Use [`project`](@ref) instead.
"""
struct ProjectedCoords{TC <: AbstractSkyCoords, T <: Real} <: AbstractProjectedCoords
    origin::TC
    offset::SVector{2, T}
end

origin(c::ProjectedCoords) = c.origin
lon(c::AbstractProjectedCoords) = lon(origin(c)) + c.offset[1] / cos(lat(origin(c)))
lat(c::AbstractProjectedCoords) = lat(origin(c)) + c.offset[2]

"""
    project(origin::AbstractSkyCoords, c::AbstractSkyCoords) -> ProjectedCoords

Represent `c` as a [`ProjectedCoords`](@ref) offset from `origin`.

`c` is first converted to `origin`'s frame, so the offset is expressed in that frame.
"""
function project(origin::AbstractSkyCoords, c::AbstractSkyCoords)
    cc = convert(typeof(origin), c)
    Δlon = rem2pi(lon(cc) - lon(origin), RoundNearest)
    offset = SVector(Δlon * cos(lat(origin)), lat(cc) - lat(origin))
    return ProjectedCoords(origin, offset)
end

# A projected coordinate's frame is its origin's, so the frame-change
# primitive delegates to it; every `convert` out of a ProjectedCoords then
# works generically (`coords2cart` sees the projected lon/lat above).
frame_transform(::Type{T}, ::Type{<:ProjectedCoords{TC}}, v) where {T <: AbstractSkyCoords, TC <: AbstractSkyCoords} =
    frame_transform(T, TC, v)
# Disambiguate against FK4Coords's own generic-source frame_transform method
frame_transform(::Type{<:FK4Coords{e}}, ::Type{<:ProjectedCoords{TC}}, v) where {e, TC <: AbstractSkyCoords} =
    frame_transform(FK4Coords{e}, TC, v)

# A ProjectedCoords target can never be satisfied by `convert`. The origin
# is a value, not a type parameter. Every conversion pathway checks its
# target frame through `_checkframe`, so one override turns them all into a
# clear error (instances already matching the target still convert by identity without reaching this).
_checkframe(::Type{TC}) where {TC <: AbstractProjectedCoords} = throw(ArgumentError(
    "Cannot `convert` into a projected coordinate type; construct one with " *
    "`project(origin, c)` instead.",
))

# The Cartesian representation is tagged with the origin's frame.
# `constructorof(typeof(c))` (the generic tag choice) must keep returning
# `ProjectedCoords` for `setproperties` reconstruction, but the frame this vector lives in is the origin's.
cartesian(c::AbstractProjectedCoords) = CartesianCoords{constructorof(typeof(origin(c)))}(coords2cart(c))

# Equality cannot fall back to the generic frame-tag + lon/lat comparison in
# types.jl: the bare `ProjectedCoords` tag erases the origin's frame, so raw
# lon/lat numbers from different origin frames would compare equal. Compare
# the origin (frame included) and the offset instead, with a matching `hash`.
Base.:(==)(a::AbstractProjectedCoords, b::AbstractProjectedCoords) =
    origin(a) == origin(b) && a.offset == b.offset
Base.hash(c::AbstractProjectedCoords, h::UInt) = hash(c.offset, hash(origin(c), h))
Base.isapprox(a::ProjectedCoords, b::ProjectedCoords; kwargs...) =
    isapprox(origin(a), origin(b); kwargs...) && isapprox(a.offset, b.offset; kwargs...)
