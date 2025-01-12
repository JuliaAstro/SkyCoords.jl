module MakieExt

using Makie
using SkyCoords: AbstractSkyCoords, lonlat

Makie.convert_arguments(ct::PointBased, c::AbstractSkyCoords) = convert_arguments(ct, [c])

Makie.convert_arguments(ct::PointBased, cs::AbstractVector{<:AbstractSkyCoords}) =
    convert_arguments(ct, map(c -> Point(lonlat(c)...), cs))

end
