abstract type AbstractProjectedCoords <: AbstractSkyCoords end

""" Projected coordinates with origin and offset.
Approximation for small FoV. """
struct ProjectedCoords{TC,T} <: AbstractProjectedCoords
    origin::TC
    offset::SVector{2,T}
end

origin(c::ProjectedCoords) = c.origin
lon(c::AbstractProjectedCoords) = lon(origin(c)) + c.offset[1] / cos(lat(origin(c)))
lat(c::AbstractProjectedCoords) = lat(origin(c)) + c.offset[2]

Base.convert(::Type{T}, c::T) where {T<:AbstractProjectedCoords} = c
Base.convert(::Type{TCto}, c::AbstractProjectedCoords) where {TCto<:AbstractSkyCoords} =
    convert(TCto,
        constructorof(typeof(origin(c)))(
            lon(c), lat(c)
        )
    )

Base.isapprox(a::ProjectedCoords, b::ProjectedCoords; kwargs...) = isapprox(origin(a), origin(b); kwargs...) && isapprox(a.offset, b.offset; kwargs...)

function project(origin::AbstractSkyCoords, c::AbstractSkyCoords)
    cc = convert(typeof(origin), c)
    Δlon = rem2pi(lon(cc) - lon(origin), RoundNearest)
    offset = SVector(Δlon * cos(lat(origin)), lat(cc) - lat(origin))
    ProjectedCoords(origin, offset)
end

