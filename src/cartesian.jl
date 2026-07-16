"""
    CartesianCoords{TC <: AbstractSkyCoords, TF} <: AbstractSkyCoords
    CartesianCoords{TC}(x, y, z)
    CartesianCoords(c::AbstractSkyCoords)

Cartesian representation of arbitrary coordinate systems.
The type paramter `TC` identifies the coordinate scheme being encoded.
Instances of this type should be created using the [`cartesian`](@ref) function.

Note: since all sky coordinates are angle measures, the Cartesian coordinates
are assumed to lie on the unit sphere.
"""
struct CartesianCoords{TC <: AbstractSkyCoords, TF <: Real} <: AbstractSkyCoords
    vec::SVector{3, TF}
end

CartesianCoords{TC}(args::Real...) where {TC} = CartesianCoords{TC}(SVector(float.(args)...))
CartesianCoords{TC}(vec::AbstractVector{TF}) where {TC, TF} = CartesianCoords{TC, TF}(vec)
CartesianCoords(c::AbstractSkyCoords) = convert(CartesianCoords, c)
CartesianCoords{TC}(c::AbstractSkyCoords) where {TC} = convert(CartesianCoords{TC}, c)
CartesianCoords{TC, TF}(c::AbstractSkyCoords) where {TC <: AbstractSkyCoords, TF} = convert(CartesianCoords{TC, TF}, c)
constructorof(::Type{<:CartesianCoords{TC}}) where {TC} = CartesianCoords{TC}

Base.vec(c::CartesianCoords) = c.vec

coords2cart(c::CartesianCoords) = vec(c)

"""
    cartesian(c::AbstractSkyCoords)

Returns a Cartesian representation of the coordinate `c` projected onto the unit sphere as a [`CartesianCoords`](@ref).

### See also
[`spherical`](@ref)
"""
cartesian(c::CartesianCoords) = c
cartesian(c::T) where {T <: AbstractSkyCoords} = CartesianCoords{T}(coords2cart(lon(c), lat(c)))

"""
    spherical(c::AbstractSkyCoords)

Returns the spherical representation of the coordinate `c`. Coordinate arguments that are already in spherical coordinates are simply returned. If `c` is a [`CartesianCoords`](@ref) object, the Cartesian representation is converted to spherical and returned as a `T <: AbstractSkyCoords` where `c::CartesianCoords{T}`. 

### See also
[`cartesian`](@ref)
"""
spherical(c::AbstractSkyCoords) = c
spherical(c::CartesianCoords{T}) where {T <: AbstractSkyCoords} = T(cart2coords(vec(c))...)

Base.convert(::Type{T}, c::CartesianCoords{S}) where {T <: AbstractSkyCoords, S <: AbstractSkyCoords} =
    spherical(convert(CartesianCoords{T}, c))
Base.convert(::Type{<:CartesianCoords{T}}, c::S) where {T <: AbstractSkyCoords, S <: AbstractSkyCoords} =
    convert(CartesianCoords{T}, cartesian(c))
Base.convert(::Type{<:CartesianCoords{T}}, c::CartesianCoords{S}) where {T <: AbstractSkyCoords, S <: AbstractSkyCoords} =
    CartesianCoords{T}(rotmat(T, S) * vec(c))


Base.:(==)(a::T, b::T) where {T <: CartesianCoords} = vec(a) == vec(b)
Base.isapprox(a::CartesianCoords{TC}, b::CartesianCoords{TC}; kwargs...) where {TC} =
    isapprox(vec(a), vec(b); kwargs...)
Base.isapprox(a::CartesianCoords{TCa}, b::CartesianCoords{TCb}; kwargs...) where {TCa, TCb} =
    isapprox(vec(a), vec(convert(CartesianCoords{TCa}, b)); kwargs...)
