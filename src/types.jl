# -----------------------------------------------------------------------------
# Types

abstract type AbstractSkyCoords end

struct ICRSCoords{T <: AbstractFloat} <: AbstractSkyCoords
    ra::T
    dec::T
    ICRSCoords{T}(ra, dec) where {T <: AbstractFloat} = new(mod2pi(ra), dec)
end
ICRSCoords(ra::T, dec::T) where {T <: AbstractFloat} = ICRSCoords{T}(ra, dec)
ICRSCoords(ra::Real, dec::Real) = ICRSCoords(promote(float(ra), float(dec))...)
ICRSCoords(ra::Quantity, dec::Quantity) = ICRSCoords(ustrip(u"rad", ra), ustrip(u"rad", dec))

struct GalCoords{T <: AbstractFloat} <: AbstractSkyCoords
    l::T
    b::T
    GalCoords{T}(l, b) where {T <: AbstractFloat} = new(mod2pi(l), b)
end
GalCoords(l::T, b::T) where {T <: AbstractFloat} = GalCoords{T}(l, b)
GalCoords(l::Real, b::Real) = GalCoords(promote(float(l), float(b))...)
GalCoords(ra::Quantity, dec::Quantity) = GalCoords(ustrip(u"rad", ra), ustrip(u"rad", dec))

# FK5 is parameterized by equinox (e)
struct FK5Coords{e,T <: AbstractFloat} <: AbstractSkyCoords
    ra::T
    dec::T
    FK5Coords{e,T}(ra, dec) where {T <: AbstractFloat,e} = new(mod2pi(ra), dec)
end
FK5Coords{e}(ra::T, dec::T) where {e,T <: AbstractFloat} = FK5Coords{e,T}(ra, dec)
FK5Coords{e}(ra::Real, dec::Real) where {e} =
   FK5Coords{e}(promote(float(ra), float(dec))...)
FK5Coords{e}(ra::Quantity, dec::Quantity) where {e} = FK5Coords{e}(ustrip(u"rad", ra), ustrip(u"rad", dec))

# Scalar coordinate conversions
convert(::Type{T}, c::T) where {T <: AbstractSkyCoords} = c
function convert(::Type{T}, c::S) where {T <: AbstractSkyCoords,S <: AbstractSkyCoords}
    r = rotmat(T, S) * coords2cart(c)
    lon, lat = cart2coords(r)
    T(lon, lat)
end