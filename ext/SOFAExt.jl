module SOFAExt

using SkyCoords

# The observed-place transforms come from SOFA.jl, used in their two-stage form:
# `apco13` prepares the star-independent astrometry context for an
# observing frame, and the quick functions transform individual coordinates
# against it. SOFA's one-shot `atco13`/`atoc13` are exactly these compositions.
using SOFA: apco13, atciqz, atioq, atoiq, aticq

const ALTAZ_TO_ALTAZ_MSG = "Converting between two horizontal frames with a single `AltAzFrame` is ambiguous. Use `AltAzCoords(c, from, to)` to name the source and destination frames."

# The astrometry context for a frame, computed on first use and cached so that
# converting many coordinates against the same frame pays for `apco13` only once.
# Concurrent first uses recompute the same value, so the unsynchronized write is benign.
function _astrom(frame::AltAzFrame)
    if frame.cache[] === nothing
        observer = frame.observer
        frame.cache[] = first(apco13(
            frame.jd, 0.0, frame.dut1,
            observer.longitude, observer.latitude, observer.altitude,
            frame.xp, frame.yp,
            frame.pressure, frame.temperature, frame.relative_humidity, frame.wavelength,
        ))
    end
    return frame.cache[]
end

# Any celestial frame --> observed horizontal coordinates
function SkyCoords.AltAzCoords(c::AbstractSkyCoords, frame::AltAzFrame)
    icrs = convert(ICRSCoords, c)
    astrom = _astrom(frame)
    ri, di = atciqz(icrs.ra, icrs.dec, astrom)
    observed = atioq(ri, di, astrom)
    return AltAzCoords(π / 2 - observed.zen, observed.azi)
end

# A horizontal direction entering the reverse transform, in either
# representation; the Cartesian form converts out through its spherical one
const HorizontalCoords = Union{AltAzCoords, CartesianCoords{<:AltAzCoords}}

# Observed horizontal coordinates --> any celestial frame
function _altaz_to(::Type{T}, c::HorizontalCoords, frame::AltAzFrame) where {T <: AbstractSkyCoords}
    sph = spherical(c)
    astrom = _astrom(frame)
    ri, di = atoiq('A', sph.az, π / 2 - sph.alt, astrom)
    ra, dec = aticq(ri, di, astrom)
    return convert(T, ICRSCoords(ra, dec))
end

(::Type{T})(c::HorizontalCoords, frame::AltAzFrame) where {T <: AbstractSkyCoords} =
    _altaz_to(T, c, frame)

# The two-argument inner constructors (`ICRSCoords{T}(ra, dec)` and friends)
# accept any two arguments, making each of them ambiguous with the generic
# reverse above. Cover every concrete family explicitly.
(::Type{ICRSCoords{F}})(c::HorizontalCoords, frame::AltAzFrame) where {F <: Real} =
    _altaz_to(ICRSCoords{F}, c, frame)
(::Type{GalCoords{F}})(c::HorizontalCoords, frame::AltAzFrame) where {F <: Real} =
    _altaz_to(GalCoords{F}, c, frame)
(::Type{SuperGalCoords{F}})(c::HorizontalCoords, frame::AltAzFrame) where {F <: Real} =
    _altaz_to(SuperGalCoords{F}, c, frame)
(::Type{FK4Coords{e, F}})(c::HorizontalCoords, frame::AltAzFrame) where {e, F <: Real} =
    _altaz_to(FK4Coords{e, F}, c, frame)
(::Type{FK4NoETerms{e, F}})(c::HorizontalCoords, frame::AltAzFrame) where {e, F <: Real} =
    _altaz_to(FK4NoETerms{e, F}, c, frame)
(::Type{FK5Coords{e, F}})(c::HorizontalCoords, frame::AltAzFrame) where {e, F <: Real} =
    _altaz_to(FK5Coords{e, F}, c, frame)
(::Type{EclipticCoords{e, F}})(c::HorizontalCoords, frame::AltAzFrame) where {e, F <: Real} =
    _altaz_to(EclipticCoords{e, F}, c, frame)
# Projected targets are rejected inside `convert` with an informative error
(::Type{ProjectedCoords{TC, F}})(c::HorizontalCoords, frame::AltAzFrame) where {TC <: AbstractSkyCoords, F <: Real} =
    _altaz_to(ProjectedCoords{TC, F}, c, frame)

# Between two horizontal frames, through the celestial sphere
SkyCoords.AltAzCoords(c::HorizontalCoords, from::AltAzFrame, to::AltAzFrame) =
    AltAzCoords(ICRSCoords(c, from), to)

# Shorthands constructing the frame in place. These build a fresh frame per call,
# so batch conversions should construct one `AltAzFrame` and reuse it.
SkyCoords.AltAzCoords(c::AbstractSkyCoords, observer::Observer, jd::Real; kwargs...) =
    AltAzCoords(c, AltAzFrame(observer, jd; kwargs...))
(::Type{T})(c::AltAzCoords, observer::Observer, jd::Real; kwargs...) where {T <: AbstractSkyCoords} =
    T(c, AltAzFrame(observer, jd; kwargs...))

# A single observing context cannot define both ends of an AltAz --> AltAz
# conversion. These gates also disambiguate the `Type{AltAzCoords}`
# constructor methods above from the generic `(::Type{T})` reverses.
SkyCoords.AltAzCoords(c::HorizontalCoords, frame::AltAzFrame) =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))
(::Type{AltAzCoords{F}})(c::HorizontalCoords, frame::AltAzFrame) where {F <: Real} =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))
SkyCoords.AltAzCoords(c::AltAzCoords, observer::Observer, jd::Real; kwargs...) =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))
(::Type{T})(c::AltAzCoords, observer::Observer, jd::Real; kwargs...) where {T <: AltAzCoords} =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))

end # module
