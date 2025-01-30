struct CartesianCoords{TC<:AbstractSkyCoords, TF<:Real} <: AbstractSkyCoords
    vec::SVector{3,TF}
end

CartesianCoords{TC}(args::Real...) where {TC} = CartesianCoords{TC}(SVector(float.(args)...))
CartesianCoords{TC}(vec::AbstractVector{TF}) where {TC, TF} = CartesianCoords{TC, TF}(vec)
CartesianCoords(c::AbstractSkyCoords) = convert(CartesianCoords, c)
CartesianCoords{TC}(c::AbstractSkyCoords) where {TC} = convert(CartesianCoords{TC}, c)
CartesianCoords{TC,TF}(c::AbstractSkyCoords) where {TC<:AbstractSkyCoords,TF} = convert(CartesianCoords{TC,TF}, c)
constructorof(::Type{<:CartesianCoords{TC}}) where {TC} = CartesianCoords{TC}

Base.vec(c::CartesianCoords) = c.vec

cartesian(c::CartesianCoords) = c
spherical(c::AbstractSkyCoords) = c
cartesian(c::T) where {T <: AbstractSkyCoords} = CartesianCoords{T}(coords2cart(lon(c), lat(c)))
spherical(c::CartesianCoords{T}) where {T <: AbstractSkyCoords} = T(cart2coords(vec(c))...)

Base.convert(::Type{<:T}, c::CartesianCoords{S}) where {T<:AbstractSkyCoords,S<:AbstractSkyCoords} =
    spherical(convert(CartesianCoords{T}, c))
Base.convert(::Type{<:CartesianCoords{T}}, c::S) where {T<:AbstractSkyCoords,S<:AbstractSkyCoords} =
    convert(CartesianCoords{T}, cartesian(c))
Base.convert(::Type{<:CartesianCoords{T}}, c::CartesianCoords{S}) where {T<:AbstractSkyCoords,S<:AbstractSkyCoords} =
    CartesianCoords{T}(rotmat(T, S) * vec(c))


Base.:(==)(a::T, b::T) where {T<:CartesianCoords} = vec(a) == vec(b)
Base.isapprox(a::CartesianCoords{TC}, b::CartesianCoords{TC}; kwargs...) where {TC} =
    isapprox(vec(a), vec(b); kwargs...)
Base.isapprox(a::CartesianCoords{TCa}, b::CartesianCoords{TCb}; kwargs...) where {TCa,TCb} =
    isapprox(vec(a), vec(convert(CartesianCoords{TCa}, b)); kwargs...)
