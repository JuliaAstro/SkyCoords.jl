module SkyCoords
using StaticArrays
using Unitful

export AbstractSkyCoords,
       ICRSCoords,
       GalCoords,
       FK5Coords,
       separation,
       @coord_str

include("types.jl")
include("utils.jl")

# ------------------------------------------------------------------------------
# generate coordinates given string

"""
    @coord_str input

Parses strings that specify common astronomical coordinates into radians. The two supported formats are `XhYmZs` and `X°Y'Z"`

# Examples
```jldoctest
julia> coord"12h52m64.300s"
3.3731614843033575

julia> coord"0°12'5\\\"" # Note that you have to escape the "
0.003514899188044136

julia> coord"-  10° 02 '  10.885 \\\"" # Whitespace does not matter
-0.17389837681291273
```
"""
macro coord_str(input::String)
    input = replace(strip(input), r"\s" => "")
    ha_r = r"([+-]?\d+)h(\d+)m(\d+\.?\d*)s"
    deg_r = r"([+-]?\d+)°(\d+)'(\d+\.?\d*)\""
    if occursin(ha_r, input)
        m = match(ha_r, input)
        h, m, s = parse.(Float64, m.captures)
        rad = h * 2π / 24
        rad += m * 2π / 24 / 60
        rad += s * 2π / 24 / 60 / 60
        return rad
    elseif occursin(deg_r, input)
        m = match(deg_r, input)
        d, m, s = parse.(Float64, m.captures)
        rad = deg2rad(d)
        rad += deg2rad(m / 60)
        rad += deg2rad(s / 60 / 60)
        return rad
    else
        error("Could not parse $input to sky coordinates")
    end
end

# ------------------------------------------------------------------------------
# Distance between coordinates

function _separation(λ_1, ϕ_1, λ_2, ϕ_2)
    Δλ = λ_2 - λ_1
    sin_Δλ, cos_Δλ = sincos(Δλ)
    sin_ϕ1, cos_ϕ1 = sincos(ϕ_1)
    sin_ϕ2, cos_ϕ2 = sincos(ϕ_2)
    return atan(hypot(cos_ϕ2 * sin_Δλ,
                      cos_ϕ1 * sin_ϕ2 - sin_ϕ1 * cos_ϕ2 * cos_Δλ),
                sin_ϕ1 * sin_ϕ2 + cos_ϕ1 * cos_ϕ2 * cos_Δλ)
end

"""
    separation(c1::AbstractSkyCoords, c2::AbstractSkyCoords) -> distance

Return angular separation between two sky coordinates, in radians.

The angular separation is calculated using the [Vincenty formula](http://en.wikipedia.org/wiki/Great-circle_distance), which is slightly more
complex and computationally expensive than some alternatives, but is stable at
at all distances, including the poles and antipodes.
"""
separation(c1::T, c2::T) where {T <: AbstractSkyCoords} =
    _separation(lon(c1), lat(c1), lon(c2), lat(c2))

separation(c1::T1, c2::T2) where {T1 <: AbstractSkyCoords,T2 <: AbstractSkyCoords} =
    separation(c1, convert(T1, c2))

end # module
