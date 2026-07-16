# -----------------------------------------------------------------------------
# Types

"""
    abstract type AbstractSkyCoords

The supertype for all sky coordinate systems.
"""
abstract type AbstractSkyCoords end

"""
    ICRSCoords <: AbstractSkyCoords
    ICRSCoords(ra, dec)

[International Celestial Reference System](https://en.wikipedia.org/wiki/International_Celestial_Reference_System)

This is the current standard adopted by the International Astronomical Union notably due to its high level of accuracy compared to standard equatorial coordinate systems. What sets this apart from [`FK5Coords`](@ref) is that it is completely defined using extragalactic radio sources rather than a geocentric frame, which means the reference frame will not change due to Earth's motion.

### Coordinates
- `ra` - Right ascension in radians (0, 2π)
- `dec` - Declination in radians (-π/2, π/2)
"""
struct ICRSCoords{T <: Real} <: AbstractSkyCoords
    ra::T
    dec::T
    ICRSCoords{T}(ra, dec) where {T <: Real} = new(mod2pi(ra), dec)
end
ICRSCoords(ra::T, dec::T) where {T <: Real} = ICRSCoords{float(T)}(ra, dec)
ICRSCoords(ra::Real, dec::Real) = ICRSCoords(promote(ra, dec)...)
ICRSCoords(c::T) where {T <: AbstractSkyCoords} = convert(ICRSCoords, c)
ICRSCoords{F}(c::T) where {F, T <: AbstractSkyCoords} = convert(ICRSCoords{F}, c)

"""
    GalCoords <: AbstractSkyCoords
    GalCoords(l, b)

[Galactic Coordinate System](https://en.wikipedia.org/wiki/Galactic_coordinate_system)

This coordinate system is defined based on the projection of the Milky Way galaxy onto our celestial sphere, with (0, 0) being approximately the center of our galaxy.

### Coordinates
- `l` - Galactic longitude in radians (-π, π)
- `b` - Galactic latitude in radians (-π/2, π/2)
"""
struct GalCoords{T <: Real} <: AbstractSkyCoords
    l::T
    b::T
    GalCoords{T}(l, b) where {T <: Real} = new(mod2pi(l), b)
end
GalCoords(l::T, b::T) where {T <: Real} = GalCoords{float(T)}(l, b)
GalCoords(l::Real, b::Real) = GalCoords(promote(l, b)...)
GalCoords(c::T) where {T <: AbstractSkyCoords} = convert(GalCoords, c)
GalCoords{F}(c::T) where {F, T <: AbstractSkyCoords} = convert(GalCoords{F}, c)

"""
    SuperGalCoords(l, b)

[Supergalactic Coordinate System](https://en.wikipedia.org/wiki/Supergalactic_coordinate_system)

The supergalactic plane is part of a reference frame for the supercluster of galaxies that contains the Milky Way galaxy.
The supergalactic plane as so-far observed is more or less perpendicular to the plane of the Milky Way, the angle is 84.5 degrees. Viewed from the Earth, the plane traces a great circle across the sky through the constellations

### Coordinates
- `l` - SuperGalCoords longitude in radians (-π, π)
- `b` - SuperGalCoords latitude in radians (-π/2, π/2)
"""
struct SuperGalCoords{T <: Real} <: AbstractSkyCoords
    l::T
    b::T
    SuperGalCoords{T}(l, b) where {T <: Real} = new(mod2pi(l), b)
end
SuperGalCoords(l::T, b::T) where {T <: Real} = SuperGalCoords{float(T)}(l, b)
SuperGalCoords(l::Real, b::Real) = SuperGalCoords(promote(l, b)...)
SuperGalCoords(c::T) where {T <: AbstractSkyCoords} = convert(SuperGalCoords, c)
SuperGalCoords{F}(c::T) where {F, T <: AbstractSkyCoords} = convert(SuperGalCoords{F}, c)


"""
    FK5Coords{e, T <: Real} <: AbstractSkyCoords
    FK5Coords{equinox}(ra, dec)

[Equatorial Coordinate System](https://en.wikipedia.org/wiki/Equatorial_coordinate_system)

This coordinate system maps the celestial sphere based on a geocentric observer.
Historically the oldest, this coordinate system has been shown to be inaccurate
due to its definitions based on the Earth, which has long-scale precession
causing the reference frame to change. Because of this, an equinox must be
provided (typically 2000, commonly known as J2000) which defines the reference
frame.

### Coordinates
- `ra` - Right ascension in radians (0, 2π)
- `dec` - Declination in radians (-π/2, π/2)
"""
struct FK5Coords{e, T <: Real} <: AbstractSkyCoords
    ra::T
    dec::T
    FK5Coords{e, T}(ra, dec) where {T <: Real, e} = new(mod2pi(ra), dec)
end
FK5Coords{e}(ra::T, dec::T) where {e, T <: Real} = FK5Coords{e, float(T)}(ra, dec)
FK5Coords{e}(ra::Real, dec::Real) where {e} = FK5Coords{e}(promote(ra, dec)...)
FK5Coords{e}(c::T) where {e, T <: AbstractSkyCoords} = convert(FK5Coords{e}, c)
FK5Coords{e, F}(c::T) where {e, F, T <: AbstractSkyCoords} = convert(FK5Coords{e, F}, c)
constructorof(::Type{<:FK5Coords{e}}) where {e} = FK5Coords{e}


"""
    EclipticCoords{equinox}(lon, lat)

[Ecliptic Coordinate System](https://en.wikipedia.org/wiki/Ecliptic_coordinate_system)

This coordinate system is geocentric with the ecliptic plane as the xy-plane with x oriented according to the equinox specified by `equinox`.

### Coordinates
- `lon` - Longitude in radians (0, 2π)
- `lat` - Latitude in radians (-π/2, π/2)
"""
struct EclipticCoords{e, T <: Real} <: AbstractSkyCoords
    lon::T
    lat::T
    EclipticCoords{e, T}(lon, lat) where {e, T <: Real} = new(mod2pi(lon), lat)
end
EclipticCoords{e}(lon::T, lat::T) where {e, T <: Real} = EclipticCoords{e, float(T)}(lon, lat)
EclipticCoords{e}(lon::Real, lat::Real) where {e} = EclipticCoords{e}(promote(lon, lat)...)
EclipticCoords{e}(c::AbstractSkyCoords) where {e} = convert(EclipticCoords{e}, c)
EclipticCoords{e, F}(c::AbstractSkyCoords) where {e, F} = convert(EclipticCoords{e, F}, c)
constructorof(::Type{<:EclipticCoords{e}}) where {e} = EclipticCoords{e}

"""
    AltAzCoords <: AbstractSkyCoords
    AltAzCoords(alt, az)

[Horizontal Coordinate System](https://en.wikipedia.org/wiki/Horizontal_coordinate_system)

This coordinate system uses the observer's local horizon as the fundamental plane to define an object in the local sky, making it useful for planning observations with azimuthal mount telescopes.

Unlike the other coordinate systems, mapping to and from horizontal coordinates requires knowing where and when the observation takes place. Conversions are therefore not available through `convert`; instead, load [Astrometry.jl](https://github.com/JuliaAstro/Astrometry.jl) and construct coordinates with an [`Observer`](@ref) and a UTC Julian date:

```julia
using SkyCoords, Astrometry

m13 = ICRSCoords(4.3713, 0.6364)
mt_wilson = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742)
jd = 2459526.6667  # 2021-11-08T04:00 UTC

altaz = AltAzCoords(m13, mt_wilson, jd)
icrs = ICRSCoords(altaz, mt_wilson, jd)
```

These conversions are performed with the IAU SOFA algorithms ([`Astrometry.SOFA.atco13`](https://juliaastro.org/Astrometry/stable/api/#Astrometry.SOFA.atco13-NTuple{18,%20AbstractFloat}) and [`Astrometry.SOFA.atoc13`](https://juliaastro.org/Astrometry/stable/api/#Astrometry.SOFA.atoc13-Tuple{Char,%20Vararg{AbstractFloat,%2014}})) and accept the following keyword arguments:

- `dut1 = 0` - UT1 - UTC in seconds (from IERS bulletins)
- `xp = 0`, `yp = 0` - Polar motion coordinates in radians (from IERS bulletins)
- `pressure = 0` - Atmospheric pressure at the observer in hPa. The default of `0` disables atmospheric refraction; set to the ambient pressure (~1000 hPa at sea level) to include it
- `temperature = 0` - Ambient temperature at the observer in °C (used for refraction)
- `relative_humidity = 0` - Relative humidity at the observer, 0-1 (used for refraction)
- `wavelength = 0.55` - Observing wavelength in μm (used for refraction)

### Coordinates
- `alt` - Altitude (elevation) angle above the observer's local horizon in radians (-π/2, π/2)
- `az` - Azimuth angle of the object around the horizon, increasing eastward from true north in radians (0, 2π)
"""
struct AltAzCoords{T <: Real} <: AbstractSkyCoords
    alt::T
    az::T
    AltAzCoords{T}(alt, az) where {T <: Real} = new(alt, mod2pi(az))
end
AltAzCoords(alt::T, az::T) where {T <: Real} = AltAzCoords{float(T)}(alt, az)
AltAzCoords(alt::Real, az::Real) = AltAzCoords(promote(alt, az)...)
AltAzCoords(c::T) where {T <: AbstractSkyCoords} = convert(AltAzCoords, c)
AltAzCoords{F}(c::T) where {F, T <: AbstractSkyCoords} = convert(AltAzCoords{F}, c)

"""
    Observer(latitude, longitude, altitude = 0)

The geodetic (WGS84) location of an observer on Earth, used for conversions to and from [`AltAzCoords`](@ref).

### Coordinates
- `latitude` - Geodetic latitude in radians (-π/2, π/2)
- `longitude` - Longitude in radians, east-positive (-π, π)
- `altitude` - Height above the WGS84 reference ellipsoid in meters
"""
struct Observer{T <: Real}
    latitude::T
    longitude::T
    altitude::T
    Observer{T}(latitude, longitude, altitude) where {T <: Real} = new(latitude, longitude, altitude)
end
Observer(latitude::T, longitude::T, altitude::T) where {T <: Real} =
    Observer{float(T)}(latitude, longitude, altitude)
Observer(latitude::Real, longitude::Real, altitude::Real = 0) =
    Observer(promote(latitude, longitude, altitude)...)

# Scalar coordinate conversions
Base.convert(::Type{T}, c::T) where {T <: AbstractSkyCoords} = c

# Changing the float type within horizontal coordinates is fine (and preserves the
# (alt, az) field order); any other conversion requires an observer location and
# time, and raises an informative error in `rotmat`.
Base.convert(::Type{T}, c::AltAzCoords) where {T <: AltAzCoords} = T(c.alt, c.az)

function Base.convert(::Type{T}, c::S) where {T <: AbstractSkyCoords, S <: AbstractSkyCoords}
    r = rotmat(T, S) * coords2cart(c)
    lon, lat = cart2coords(r)
    return T(lon, lat)
end

Base.:(==)(a::T, b::T) where {T <: AbstractSkyCoords} = lon(a) == lon(b) && lat(a) == lat(b)
Base.isapprox(a::ICRSCoords, b::ICRSCoords; kwargs...) = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::GalCoords, b::GalCoords; kwargs...) = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::SuperGalCoords, b::SuperGalCoords; kwargs...) = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::FK5Coords{e}, b::FK5Coords{e}; kwargs...) where {e} = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::EclipticCoords{e}, b::EclipticCoords{e}; kwargs...) where {e} = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
Base.isapprox(a::AltAzCoords, b::AltAzCoords; kwargs...) = isapprox(SVector(lon(a), lat(a)), SVector(lon(b), lat(b)); kwargs...)
