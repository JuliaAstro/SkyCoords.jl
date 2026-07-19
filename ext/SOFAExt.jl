module SOFAExt

using SkyCoords

# The observed-place transforms currently come from SOFA.jl
using SOFA: atco13, atoc13

const ALTAZ_TO_ALTAZ_MSG = "Converting between two sets of horizontal coordinates requires the observer location and time of both; convert through `ICRSCoords` instead."

# ICRS (ra, dec) --> observed (azimuth, zenith distance), all in radians
function _icrs_to_observed(ra, dec, observer, jd; dut1, xp, yp, pressure, temperature, relative_humidity, wavelength)
    obs = atco13(
        float(ra), float(dec), 0.0, 0.0, 0.0, 0.0,
        float(jd), 0.0, float(dut1),
        float(observer.longitude), float(observer.latitude), float(observer.altitude),
        float(xp), float(yp),
        float(pressure), float(temperature), float(relative_humidity), float(wavelength),
    )
    return obs.azi, obs.zen
end

# Observed (azimuth, zenith distance) --> ICRS (ra, dec), all in radians
function _observed_to_icrs(az, zen, observer, jd; dut1, xp, yp, pressure, temperature, relative_humidity, wavelength)
    return atoc13(
        'A', float(az), float(zen),
        float(jd), 0.0, float(dut1),
        float(observer.longitude), float(observer.latitude), float(observer.altitude),
        float(xp), float(yp),
        float(pressure), float(temperature), float(relative_humidity), float(wavelength),
    )
end

# Any celestial frame --> observed horizontal coordinates
function SkyCoords.AltAzCoords(
        c::AbstractSkyCoords, observer::Observer, jd::Real;
        dut1 = 0, xp = 0, yp = 0,
        pressure = 0, temperature = 0, relative_humidity = 0, wavelength = 0.55,
    )
    icrs = convert(ICRSCoords, c)
    azi, zen = _icrs_to_observed(
        icrs.ra, icrs.dec, observer, jd;
        dut1, xp, yp, pressure, temperature, relative_humidity, wavelength,
    )
    return AltAzCoords(π / 2 - zen, azi)
end

# Observed horizontal coordinates --> any celestial frame
function (::Type{T})(
        c::AltAzCoords, observer::Observer, jd::Real;
        dut1 = 0, xp = 0, yp = 0,
        pressure = 0, temperature = 0, relative_humidity = 0, wavelength = 0.55,
    ) where {T <: AbstractSkyCoords}
    ra, dec = _observed_to_icrs(
        c.az, π / 2 - c.alt, observer, jd;
        dut1, xp, yp, pressure, temperature, relative_humidity, wavelength,
    )
    return convert(T, ICRSCoords(ra, dec))
end

(::Type{T})(c::AltAzCoords, observer::Observer, jd::Real; kwargs...) where {T <: AltAzCoords} =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))
SkyCoords.AltAzCoords(c::AltAzCoords, observer::Observer, jd::Real; kwargs...) =
    throw(ArgumentError(ALTAZ_TO_ALTAZ_MSG))

end # module
