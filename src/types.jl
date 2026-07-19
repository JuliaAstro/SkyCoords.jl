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
    FK4Coords{e, T <: Real} <: AbstractSkyCoords
    FK4Coords{equinox}(ra, dec)

[B1950 Equatorial Coordinate System](https://en.wikipedia.org/wiki/Epoch_(astronomy)#B1950.0)

The predecessor to [`FK5Coords`](@ref). Like `FK5Coords`, this system is defined relative to the mean equator and equinox of a given epoch, but that epoch is Besselian rather than Julian (traditionally 1950, i.e. B1950). Unlike every other coordinate system in this package, `FK4Coords` also carries the E-terms of elliptic aberration, a small (~20″) position-dependent correction that was folded directly into the original FK4 catalog positions. Because this correction is not a rotation, `FK4Coords` implements its own `frame_transform` methods rather than providing a `rotmat`; all conversions — including through [`CartesianCoords`](@ref) — otherwise work like any other system. See [`FK4NoETerms`](@ref) for the E-terms-free, purely rotational variant.

### Coordinates
- `ra` - Right ascension in radians (0, 2π)
- `dec` - Declination in radians (-π/2, π/2)
"""
struct FK4Coords{e, T <: Real} <: AbstractSkyCoords
    ra::T
    dec::T
    FK4Coords{e, T}(ra, dec) where {T <: Real, e} = new(mod2pi(ra), dec)
end
FK4Coords{e}(ra::T, dec::T) where {e, T <: Real} = FK4Coords{e, float(T)}(ra, dec)
FK4Coords{e}(ra::Real, dec::Real) where {e} = FK4Coords{e}(promote(ra, dec)...)
FK4Coords{e}(c::T) where {e, T <: AbstractSkyCoords} = convert(FK4Coords{e}, c)
FK4Coords{e, F}(c::T) where {e, F, T <: AbstractSkyCoords} = convert(FK4Coords{e, F}, c)
constructorof(::Type{<:FK4Coords{e}}) where {e} = FK4Coords{e}


"""
    FK4NoETerms{e, T <: Real} <: AbstractSkyCoords
    FK4NoETerms{equinox}(ra, dec)

[`FK4Coords`](@ref) with the E-terms of elliptic aberration removed, leaving a purely rotational equatorial system defined relative to the mean equator and equinox of a given Besselian epoch (traditionally 1950, i.e. B1950). This is astropy's `FK4NoETerms` frame.

### Coordinates
- `ra` - Right ascension in radians (0, 2π)
- `dec` - Declination in radians (-π/2, π/2)
"""
struct FK4NoETerms{e, T <: Real} <: AbstractSkyCoords
    ra::T
    dec::T
    FK4NoETerms{e, T}(ra, dec) where {T <: Real, e} = new(mod2pi(ra), dec)
end
FK4NoETerms{e}(ra::T, dec::T) where {e, T <: Real} = FK4NoETerms{e, float(T)}(ra, dec)
FK4NoETerms{e}(ra::Real, dec::Real) where {e} = FK4NoETerms{e}(promote(ra, dec)...)
FK4NoETerms{e}(c::T) where {e, T <: AbstractSkyCoords} = convert(FK4NoETerms{e}, c)
FK4NoETerms{e, F}(c::T) where {e, F, T <: AbstractSkyCoords} = convert(FK4NoETerms{e, F}, c)
constructorof(::Type{<:FK4NoETerms{e}}) where {e} = FK4NoETerms{e}


"""
    AltAzCoords <: AbstractSkyCoords
    AltAzCoords(alt, az)

[Horizontal Coordinate System](https://en.wikipedia.org/wiki/Horizontal_coordinate_system)

This coordinate system uses the observer's local horizon as the fundamental plane to define an object in the local sky, making it useful for planning observations with azimuthal mount telescopes.

Unlike the other coordinate systems, mapping to and from horizontal coordinates requires knowing where and when the observation takes place. Conversions are therefore not available through `convert`. Instead, load [SOFA.jl](@extref SOFA :doc:`index`), describe the observing context with an [`AltAzFrame`](@ref), and pass the frame to the coordinate constructors:

```jldoctest
julia> using SkyCoords, SOFA

julia> m13 = ICRSCoords(4.3713, 0.6364)
ICRSCoords{Float64}(4.3713, 0.6364)

julia> mt_wilson = Observer(deg2rad(34.2247), deg2rad(-118.0572), 1742)
Observer{Float64}(0.5973337005073033, -2.0604868456854497, 1742.0)

julia> frame = AltAzFrame(mt_wilson, 2459526.6667)  # 2021-11-08T04:00 UTC
AltAzFrame(Observer{Float64}(0.5973337005073033, -2.0604868456854497, 1742.0), 2.4595266667e6)

julia> altaz = AltAzCoords(m13, frame)
AltAzCoords{Float64}(0.23340717371373332, 5.326718442771835)

julia> icrs = ICRSCoords(altaz, frame)
ICRSCoords{Float64}(4.371299999999999, 0.6364000000000003)
```

`AltAzCoords(c, observer, jd; kwargs...)` and `ICRSCoords(c, observer, jd; kwargs...)` are shorthands that construct the [`AltAzFrame`](@ref) in place, and `AltAzCoords(c, from, to)` maps a coordinate between two horizontal frames. The transforms use the IAU SOFA observed-place algorithms in their two-stage form: the star-independent astrometry context ([`SOFA.apco13`](@extref)) is computed once per frame and cached there, and each coordinate is then transformed against it ([`SOFA.atciqz`](@extref)/[`SOFA.atioq`](@extref) in, [`SOFA.atoiq`](@extref)/[`SOFA.aticq`](@extref) out), so converting many coordinates against the same frame is cheap.

The source is treated as a fixed point at infinity: proper motion, parallax, and radial velocity are assumed to be zero, matching astropy's behavior for position-only coordinates. The resulting error stays below an arcsecond for almost all stars (bounded by the parallax, plus the proper motion accumulated since the catalog epoch), but reaches arcminutes for the fastest-moving nearby stars. Solar-system objects are out of scope; their topocentric parallax alone (up to about a degree for the Moon) is not modeled.

### Coordinates
- `alt`: Altitude (elevation) angle above the observer's local horizon in radians (-π/2, π/2)
- `az`: Azimuth angle of the object around the horizon, increasing eastward from true north in radians (0, 2π)
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

The geodetic (WGS84) location of an observer on Earth. Together with an observation time, it forms an [`AltAzFrame`](@ref) for conversions to and from [`AltAzCoords`](@ref).

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

"""
    AltAzFrame(observer, jd; kwargs...)

The observing context of a horizontal coordinate: an [`Observer`](@ref) location together with a UTC Julian date and optional Earth orientation and atmospheric parameters. Load [SOFA.jl](@extref SOFA :doc:`index`) and pass a frame to the coordinate constructors to convert into and out of the horizontal system it describes. See [`AltAzCoords`](@ref) for examples.

### Parameters

- `observer`: [`Observer`](@ref) location on Earth
- `jd`: UTC Julian date of the observation
- `dut1 = 0`: UT1 -- UTC in seconds (from IERS bulletins)
- `xp = 0`, `yp = 0`: Polar motion coordinates in radians (from IERS bulletins)
- `pressure = 0`: Atmospheric pressure at the observer in hPa. The default of `0` disables atmospheric refraction. Set to the ambient pressure (~1000 hPa at sea level) to include it.
- `temperature = 0`: Ambient temperature at the observer in °C (used for refraction)
- `relative_humidity = 0`: Relative humidity at the observer, 0-1 (used for refraction)
- `wavelength = 0.55`: Observing wavelength in μm (used for refraction)

All parameters are stored as `Float64`. Note the accuracy consequences of the defaults: `dut1 = 0` can misplace a coordinate by up to ~13″ of Earth rotation (|UT1 -- UTC| is kept below 0.9 s), `xp = yp = 0` by a few tenths of an arcsecond, and `pressure = 0` ignores refraction entirely, which amounts to ~1′ at 45° altitude and over half a degree at the horizon. With refraction enabled, the SOFA model is accurate to ~0.05″ (optical) above 20° altitude, degrading to arcminutes at the horizon.

The first conversion against a frame computes SOFA's star-independent astrometry context ([`SOFA.apco13`](@extref)) and caches it within the frame, so converting many coordinates against the same frame, e.g., `AltAzCoords.(catalog, frame)`, performs the expensive per-time setup only once.
"""
struct AltAzFrame
    observer::Observer{Float64}
    jd::Float64
    dut1::Float64
    xp::Float64
    yp::Float64
    pressure::Float64
    temperature::Float64
    relative_humidity::Float64
    wavelength::Float64
    # Star-independent astrometry context, computed and read only by the
    # SOFA.jl extension; `nothing` until first use. Untyped because the
    # value's type belongs to the extension.
    cache::Base.RefValue{Any}
end

function AltAzFrame(
        observer::Observer, jd::Real;
        dut1 = 0, xp = 0, yp = 0,
        pressure = 0, temperature = 0, relative_humidity = 0, wavelength = 0.55,
    )
    return AltAzFrame(
        Observer{Float64}(observer.latitude, observer.longitude, observer.altitude),
        jd, dut1, xp, yp, pressure, temperature, relative_humidity, wavelength,
        Ref{Any}(nothing),
    )
end

# Every field except the astrometry cache, which is an internal memo rather
# than part of the frame's value: frames describe the same observing context
# whether or not it has been computed yet.
_frame_fields(f::AltAzFrame) = (
    f.observer, f.jd, f.dut1, f.xp, f.yp,
    f.pressure, f.temperature, f.relative_humidity, f.wavelength,
)
Base.:(==)(a::AltAzFrame, b::AltAzFrame) = _frame_fields(a) == _frame_fields(b)
Base.hash(f::AltAzFrame, h::UInt) = hash(_frame_fields(f), hash(AltAzFrame, h))

# Print the constructor form, with only the non-default keywords and without
# the internal cache
function Base.show(io::IO, f::AltAzFrame)
    print(io, "AltAzFrame(", f.observer, ", ", f.jd)
    sep = "; "
    for (name, default) in pairs((
            dut1 = 0.0, xp = 0.0, yp = 0.0,
            pressure = 0.0, temperature = 0.0, relative_humidity = 0.0, wavelength = 0.55,
        ))
        value = getfield(f, name)
        if value != default
            print(io, sep, name, " = ", value)
            sep = ", "
        end
    end
    return print(io, ")")
end

# Context objects are scalars under broadcasting, so `AltAzFrame.(observer, jds)`
# and `AltAzCoords.(coords, frame)` work without a `Ref`
Base.broadcastable(o::Observer) = Ref(o)
Base.broadcastable(f::AltAzFrame) = Ref(f)

# Scalar coordinate conversions
Base.convert(::Type{T}, c::T) where {T <: AbstractSkyCoords} = c

function Base.convert(::Type{T}, c::S) where {T <: AbstractSkyCoords, S <: AbstractSkyCoords}
    r = frame_transform(T, S, coords2cart(c))
    lon, lat = cart2coords(r)
    return fromlonlat(T, lon, lat)
end

# Two coordinates are equal when they are in the same frame, identified by the
# element-type-free frame tag, `constructorof`, and their angles are equal.
# Element types need not match, just like `1.0 == 1.0f0`. `hash` mirrors this
# so that value-equal coordinates of different element types collide in
# `Dict`/`Set` as required by the `hash` contract.
Base.:(==)(a::AbstractSkyCoords, b::AbstractSkyCoords) = constructorof(typeof(a)) == constructorof(typeof(b)) && lonlat(a) == lonlat(b)
Base.hash(c::AbstractSkyCoords, h::UInt) = hash(lonlat(c), hash(constructorof(typeof(c)), h))

# Same-frame approximate equality compares (lon, lat), with the longitude
# difference taken in (-π, π] so that nearly-equal points on either side of
# the lon = 0 wrap still compare ≈, e.g., `ICRSCoords(eps(), 1) ≈ ICRSCoords(-eps(), 1)`.
_isapprox_lonlat(a, b; kwargs...) = isapprox(
    SVector(lon(a), lat(a)),
    SVector(lon(a) + rem2pi(lon(b) - lon(a), RoundNearest), lat(b));
    kwargs...
)
Base.isapprox(a::ICRSCoords, b::ICRSCoords; kwargs...) = _isapprox_lonlat(a, b; kwargs...)
Base.isapprox(a::GalCoords, b::GalCoords; kwargs...) = _isapprox_lonlat(a, b; kwargs...)
Base.isapprox(a::SuperGalCoords, b::SuperGalCoords; kwargs...) = _isapprox_lonlat(a, b; kwargs...)
Base.isapprox(a::FK4Coords{e}, b::FK4Coords{e}; kwargs...) where {e} = _isapprox_lonlat(a, b; kwargs...)
Base.isapprox(a::FK4NoETerms{e}, b::FK4NoETerms{e}; kwargs...) where {e} = _isapprox_lonlat(a, b; kwargs...)
Base.isapprox(a::FK5Coords{e}, b::FK5Coords{e}; kwargs...) where {e} = _isapprox_lonlat(a, b; kwargs...)
Base.isapprox(a::EclipticCoords{e}, b::EclipticCoords{e}; kwargs...) where {e} = _isapprox_lonlat(a, b; kwargs...)
Base.isapprox(a::AltAzCoords, b::AltAzCoords; kwargs...) = _isapprox_lonlat(a, b; kwargs...)
