# -----------------------------------------------------------------------------
# Types

"""
The supertype for all sky coordinate systems.
"""
abstract type AbstractSkyCoords end

"""
    ICRSCoords(ra, dec)

[International Celestial Reference System](https://en.wikipedia.org/wiki/International_Celestial_Reference_System)

This is the current standard adopted by the International Astronomical Union notably due to its high level of accuracy compared to standard equatorial coordinate systems. What sets this apart from [`FK5Coords`](@ref) is that it is completely defined using extragalactic radio sources rather than a geocentric frame, which means the reference frame will not change due to Earth's motion.

# Coordinates
- `ra` - Right ascension in radians (0, 2π)
- `dec` - Declination in radians (-π/2, π/2)
"""
struct ICRSCoords{T<:AbstractFloat} <: AbstractSkyCoords
    ra::T
    dec::T
    ICRSCoords{T}(ra, dec) where {T<:AbstractFloat} = new(mod2pi(ra), dec)
end
ICRSCoords(ra::T, dec::T) where {T<:AbstractFloat} = ICRSCoords{T}(ra, dec)
ICRSCoords(ra::Real, dec::Real) = ICRSCoords(promote(float(ra), float(dec))...)
ICRSCoords(c::T) where {T<:AbstractSkyCoords} = convert(ICRSCoords, c)
ICRSCoords{F}(c::T) where {F,T<:AbstractSkyCoords} = convert(ICRSCoords{F}, c)

"""
    GalCoords(l, b)

[Galactic Coordinate System](https://en.wikipedia.org/wiki/Galactic_coordinate_system)

This coordinate system is defined based on the projection of the Milky Way galaxy onto our celestial sphere, with (0, 0) being approximately the center of our galaxy.

# Coordinates
- `l` - Galactic longitude in radians (-π, π)
- `b` - Galactic latitude in radians (-π/2, π/2)
"""
struct GalCoords{T<:AbstractFloat} <: AbstractSkyCoords
    l::T
    b::T
    GalCoords{T}(l, b) where {T<:AbstractFloat} = new(mod2pi(l), b)
end
GalCoords(l::T, b::T) where {T<:AbstractFloat} = GalCoords{T}(l, b)
GalCoords(l::Real, b::Real) = GalCoords(promote(float(l), float(b))...)
GalCoords(c::T) where {T<:AbstractSkyCoords} = convert(GalCoords, c)
GalCoords{F}(c::T) where {F,T<:AbstractSkyCoords} = convert(GalCoords{F}, c)

"""
    FK5Coords{equinox}(ra, dec)

[Equatorial Coordinate System](https://en.wikipedia.org/wiki/Equatorial_coordinate_system)

This coordinate system maps the celestial sphere based on a geocentric observer. Historically the oldest, this coordinate system has been shown to be inaccurate due to its definitions based on the Earth, which has long-scale precession causing the reference frame to change. Because of this, an equinox must be provided (typically 2000, commonly known as J2000) which defines the reference frame.

# Coordinates
- `ra` - Right ascension in radians (0, 2π)
- `dec` - Declination in radians (-π/2, π/2)
"""
struct FK5Coords{e,T<:AbstractFloat} <: AbstractSkyCoords
    ra::T
    dec::T
    FK5Coords{e,T}(ra, dec) where {T<:AbstractFloat,e} = new(mod2pi(ra), dec)
end
FK5Coords{e}(ra::T, dec::T) where {e,T<:AbstractFloat} = FK5Coords{e,T}(ra, dec)
FK5Coords{e}(ra::Real, dec::Real) where {e} =
    FK5Coords{e}(promote(float(ra), float(dec))...)
FK5Coords{e}(c::T) where {e,T<:AbstractSkyCoords} = convert(FK5Coords{e}, c)
FK5Coords{e,F}(c::T) where {e,F,T<:AbstractSkyCoords} = convert(FK5Coords{e,F}, c)
constructorof(::Type{<:FK5Coords{e}}) where {e} = FK5Coords{e}


# Scalar coordinate conversions
Base.convert(::Type{T}, c::T) where {T<:AbstractSkyCoords} = c

function Base.convert(::Type{T}, c::S) where {T<:AbstractSkyCoords,S<:AbstractSkyCoords}
    r = rotmat(T, S) * coords2cart(c)
    lon, lat = cart2coords(r)
    T(lon, lat)
end

Base.:(==)(a::T, b::T) where {T<:AbstractSkyCoords} = lon(a) == lon(b) && lat(a) == lat(b)
Base.isapprox(a::ICRSCoords, b::ICRSCoords; kwargs...) = isapprox(
    SVector(lon(a), lat(a)),
    SVector(lon(a) + rem2pi(lon(b) - lon(a), RoundNearest), lat(b));
    kwargs...)
Base.isapprox(a::GalCoords, b::GalCoords; kwargs...) = isapprox(
    SVector(lon(a), lat(a)),
    SVector(lon(a) + rem2pi(lon(b) - lon(a), RoundNearest), lat(b));
    kwargs...)
Base.isapprox(a::FK5Coords{e}, b::FK5Coords{e}; kwargs...) where {e} = isapprox(
    SVector(lon(a), lat(a)),
    SVector(lon(a) + rem2pi(lon(b) - lon(a), RoundNearest), lat(b));
    kwargs...)
