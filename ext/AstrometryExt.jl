module AstrometryExt

using SkyCoords
using Astrometry.SOFA: atco13, atoc13

const ALTAZ_TO_ALTAZ_MSG = "Converting between two sets of horizontal coordinates requires the observer location and time of both; convert through `ICRSCoords` instead."

# Any celestial frame --> observed horizontal coordinates
function SkyCoords.AltAzCoords(
        c::AbstractSkyCoords, observer::Observer, jd::Real;
        dut1 = 0, xp = 0, yp = 0,
        pressure = 0, temperature = 0, relative_humidity = 0, wavelength = 0.55,
    )
    icrs = convert(ICRSCoords, c)
    obs = atco13(
        float(icrs.ra), float(icrs.dec), 0.0, 0.0, 0.0, 0.0,
        float(jd), 0.0, float(dut1),
        float(observer.longitude), float(observer.latitude), float(observer.altitude),
        float(xp), float(yp),
        float(pressure), float(temperature), float(relative_humidity), float(wavelength),
    )
    return AltAzCoords(π / 2 - obs.zen, obs.azi)
end

# Observed horizontal coordinates --> any celestial frame
function (::Type{T})(
        c::AltAzCoords, observer::Observer, jd::Real;
        dut1 = 0, xp = 0, yp = 0,
        pressure = 0, temperature = 0, relative_humidity = 0, wavelength = 0.55,
    ) where {T <: AbstractSkyCoords}
    ra, dec = atoc13(
        'A', float(c.az), float(π / 2 - c.alt),
        float(jd), 0.0, float(dut1),
        float(observer.longitude), float(observer.latitude), float(observer.altitude),
        float(xp), float(yp),
        float(pressure), float(temperature), float(relative_humidity), float(wavelength),
    )
    return convert(T, ICRSCoords(ra, dec))
end

(::Type{T})(c::AltAzCoords, observer::Observer, jd::Real; kwargs...) where {T <: AltAzCoords} =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))
SkyCoords.AltAzCoords(c::AltAzCoords, observer::Observer, jd::Real; kwargs...) =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))

end
