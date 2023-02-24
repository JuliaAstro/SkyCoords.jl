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
    SuperGalCoords(l, b)

[Supergalactic Coordinate System](https://en.wikipedia.org/wiki/Supergalactic_coordinate_system)

The supergalactic plane is part of a reference frame for the supercluster of galaxies that contains the Milky Way galaxy.
The supergalactic plane as so-far observed is more or less perpendicular to the plane of the Milky Way, the angle is 84.5 degrees. Viewed from the Earth, the plane traces a great circle across the sky through the constellations 

# Coordinates
- `l` - SuperGalCoords longitude in radians (-π, π)
- `b` - SuperGalCoords latitude in radians (-π/2, π/2)
"""
struct SuperGalCoords{T<:AbstractFloat} <: AbstractSkyCoords
    l::T
    b::T
    SuperGalCoords{T}(l, b) where {T<:AbstractFloat} = new(mod2pi(l), b)
end
SuperGalCoords(l::T, b::T) where {T<:AbstractFloat} = SuperGalCoords{T}(l, b)
SuperGalCoords(l::Real, b::Real) = SuperGalCoords(promote(float(l), float(b))...)
SuperGalCoords(c::T) where {T<:AbstractSkyCoords} = convert(SuperGalCoords, c)
SuperGalCoords{F}(c::T) where {F,T<:AbstractSkyCoords} = convert(SuperGalCoords{F}, c)


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


"""
    EclipticCoords{equinox}(lon, lat)

[Ecliptic Coordinate System](https://en.wikipedia.org/wiki/Ecliptic_coordinate_system)

This coordinate system is geocentric with the ecliptic plane as the xy-plane with x oriented according to the equinox specified by `equinox`.

# Coordinates
- `lon` - Longitude in radians (0, 2π)
- `lat` - Latitude in radians (-π/2, π/2)
"""
struct EclipticCoords{e,T<:AbstractFloat} <: AbstractSkyCoords
    lon::T
    lat::T
    EclipticCoords{e,T}(lon, lat) where {e,T<:AbstractFloat} = new(mod2pi(lon), lat)
end
EclipticCoords{e}(lon::T, lat::T) where {e,T<:AbstractFloat} = EclipticCoords{e,T}(lon, lat)
EclipticCoords{e}(lon::Real, lat::Real) where {e} = EclipticCoords(promote(float(lon), float(lat))...)
EclipticCoords{e}(c::AbstractSkyCoords) where {e} = convert(EclipticCoords{e}, c)
EclipticCoords{e,F}(c::AbstractSkyCoords) where {e,F} = convert(EclipticCoords{e,F}, c)
constructorof(::Type{<:EclipticCoords{e}}) where {e} = EclipticCoords{e}


# Scalar coordinate conversions
Base.convert(::Type{T}, c::T) where {T<:AbstractSkyCoords} = c

function Base.convert(::Type{T}, c::S) where {T<:AbstractSkyCoords,S<:AbstractSkyCoords}
    r = rotmat(T, S) * coords2cart(c)
    lon, lat = cart2coords(r)
    T(lon, lat)
end

Base.:(==)(a::T, b::T) where {T<:AbstractSkyCoords} = lon(a) == lon(b) && lat(a) == lat(b)
Base.isapprox(a::ICRSCoords, b::ICRSCoords; kwargs...) = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::GalCoords, b::GalCoords; kwargs...) = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::SuperGalCoords, b::SuperGalCoords; kwargs...) = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::FK5Coords{e}, b::FK5Coords{e}; kwargs...) where {e} = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::EclipticCoords{e}, b::EclipticCoords{e}; kwargs...) where {e} = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
