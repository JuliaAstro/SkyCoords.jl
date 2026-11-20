abstract type AbstractProjectedCoords <: AbstractSkyCoords end

"""
    ProjectedCoords{TC <: AbstractSkyCoords, T <: Number} <: AbstractSkyCoords

Represent a coordinate as an offset `(east, north)` from `origin`, in the azimuthal equidistant projection around it:
the point lies `norm(offset)` along the great circle that leaves `origin` at position angle `atan(offset[1], offset[2])`.
Up to `norm(offset) = π`, that is its angular distance from `origin`; longer offsets continue past the antipode.
The projection is applicable at any distance and at any origin.
For small offsets, `offset[1] ≈ Δlon * cos(lat)` and `offset[2] ≈ Δlat`.

Offsets are in radians, or angle quantities: Unitful ones are kept as given, DynamicQuantities ones converted to radians.

A projected coordinate is a _representation_ around a point, not a coordinate frame of its own.
Its frame is the origin's. Conversions out of a `ProjectedCoords` (spherical or [`CartesianCoords`](@ref) targets alike)
work like any other coordinate. Conversions into one are not possible via `convert`, since they require an origin value.
Use [`project`](@ref) instead.
"""
struct ProjectedCoords{TC <: AbstractSkyCoords, T <: Number} <: AbstractProjectedCoords
    origin::TC
    offset::SVector{2, T}

    # Typed like the spherical types' inner constructors so the two-argument
    # form stays disjoint from extension constructors taking other argument
    # pairs
    ProjectedCoords{TC, T}(origin::AbstractSkyCoords, offset::AbstractVector{<:Number}) where {TC <: AbstractSkyCoords, T <: Number} =
        new(origin, offset)
end
# Floated like the spherical types' angles, so `offset` results fit the same type
# A Cartesian origin is held in spherical form
function ProjectedCoords(origin::AbstractSkyCoords, offset::SVector{2, T}) where {T <: Number}
    foreach(_radians, offset) # Non-angle quantities throw at compile time
    o = spherical(origin)
    return ProjectedCoords{typeof(o), float(T)}(o, offset)
end

# An offset in radians. The Unitful extension converts angle quantities, which are otherwise kept as given.
_radians(x::Real) = x

"""
    origin(c::ProjectedCoords) -> coordinate

The coordinate that `c` is an offset from.
"""
origin(c::ProjectedCoords) = c.origin

# The point at distance ρ = |offset| along the great circle leaving the origin at position angle
# atan(east, north): cos ρ k̂ + (sin ρ / ρ) (east ê + north n̂), with k̂ the origin's unit vector and
# ê, n̂ its East and North unit vectors.
function coords2cart(c::AbstractProjectedCoords)
    east, north = map(_radians, c.offset) # Needed to keep the branch below type-stable
    # ρ is not differentiable at zero offset: use the series in ρ² there, so derivatives and
    # uncertainties stay finite (`float`: integers have no `eps`)
    ρ² = east^2 + north^2
    if ρ² < eps(float(one(ρ²)))^2
        cosρ, w = 1 - ρ² / 2, 1 - ρ² / 6
    else
        ρ = hypot(east, north)
        sinρ, cosρ = sincos(ρ)
        w = sinρ / ρ
    end
    # The origin at the offset's precision, which may exceed the origin's own
    λ₀, φ₀, _ = promote(lonlat(origin(c))..., w)
    k̂, ê, n̂ = _local_frame(λ₀, φ₀)
    return cosρ * k̂ + w * (east * ê + north * n̂)
end

# Longitude in [0, 2π), like the spherical types
function lonlat(c::AbstractProjectedCoords)
    λ, φ = cart2coords(coords2cart(c))
    return mod2pi(λ), φ
end
lon(c::AbstractProjectedCoords) = lonlat(c)[1]
lat(c::AbstractProjectedCoords) = lonlat(c)[2]

"""
    project(origin::AbstractSkyCoords, c::AbstractSkyCoords) -> ProjectedCoords

Represent `c` as a [`ProjectedCoords`](@ref) offset from `origin`: its separation from `origin`
times the unit vector `(sin, cos)` of its position angle.

`c` is first converted to `origin`'s frame, so the offset is expressed in that frame.

### Examples
```jldoctest
julia> o = ICRSCoords(0.5, 0.3); c = offset(o, 0.05, 0.6);

julia> SkyCoords.project(o, c).offset ≈ 0.05 .* [sin(0.6), cos(0.6)]
true

julia> convert(ICRSCoords, SkyCoords.project(o, c)) ≈ c
true
```
"""
project(origin::AbstractSkyCoords, c::AbstractSkyCoords) =
    ProjectedCoords(origin, _project_offset(origin, convert(typeof(origin), c)))
project(origin::CartesianCoords, c::AbstractSkyCoords) = project(spherical(origin), c)

# The offset from `origin` of `c`, which is already in `origin`'s frame: the separation θ times the
# unit vector (sin, cos) of the position angle.
# here, (s_E, s_N) = sin θ times that unit vector
function _project_offset(origin, c)
    s_E, s_N, cosθ = _vincenty_terms(promote(lonlat(origin)..., lonlat(c)...)...)
    # sin θ is not differentiable at the origin: use the series of θ / sin θ in sin²θ there, as in
    # `coords2cart`. At the antipode sin θ is tiny but nonzero (π rounds); `cosθ > 0` keeps it exact.
    sin²θ = s_E^2 + s_N^2
    f = if sin²θ < eps(one(sin²θ))^2 && cosθ > 0
        1 + sin²θ / 6
    else
        sinθ = hypot(s_E, s_N)
        atan(sinθ, cosθ) / sinθ
    end
    return SVector(f * s_E, f * s_N)
end

# A projected coordinate's frame is its origin's, so the frame-change
# primitive delegates to it; every `convert` out of a ProjectedCoords then
# works generically
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

# The type `convert` can target in place of a projected coordinate: its origin's
# frame, followed through nested projections, without an element type so offsets
# more precise than the origin are not rounded to it. Every other coordinate type is its own.
_unprojected_type(::Type{T}) where {T <: AbstractSkyCoords} = T
_unprojected_type(::Type{<:ProjectedCoords{TC}}) where {TC <: AbstractSkyCoords} = constructorof(_unprojected_type(TC))

# Express `r`, a coordinate in `c`'s underlying frame, the same way `c` is
# expressed. A projected coordinate keeps its origin, and so its type.
_reproject(::AbstractSkyCoords, r) = r
_reproject(c::T, r) where {T <: ProjectedCoords} = T(origin(c), _project_offset(origin(c), r))

# The Cartesian representation is tagged with the origin's frame.
# `constructorof(typeof(c))` (the generic tag choice) must keep returning
# `ProjectedCoords` for `setproperties` reconstruction, but the frame this vector lives in is the origin's.
cartesian(c::AbstractProjectedCoords) = CartesianCoords{_unprojected_type(typeof(c))}(coords2cart(c))

# Equality cannot fall back to the generic frame-tag + lon/lat comparison in
# types.jl: the bare `ProjectedCoords` tag erases the origin's frame, so raw
# lon/lat numbers from different origin frames would compare equal. Compare
# the origin (frame included) and the offset instead, with a matching `hash`.
Base.:(==)(a::AbstractProjectedCoords, b::AbstractProjectedCoords) =
    origin(a) == origin(b) && a.offset == b.offset
Base.hash(c::AbstractProjectedCoords, h::UInt) = hash(c.offset, hash(origin(c), h))
Base.isapprox(a::ProjectedCoords, b::ProjectedCoords; kwargs...) =
    isapprox(origin(a), origin(b); kwargs...) && isapprox(a.offset, b.offset; kwargs...)
