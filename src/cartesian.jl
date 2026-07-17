"""
    CartesianCoords{TC <: AbstractSkyCoords, TF <: Real} <: AbstractSkyCoords
    CartesianCoords{TC}(x, y, z)
    CartesianCoords{TC}(vec::AbstractVector)
    CartesianCoords(c::AbstractSkyCoords)

Cartesian representation of a sky coordinate as a 3-vector on the unit sphere.

The type parameter `TC` is a *frame tag*: the element-type-free constructor of
the coordinate system being represented, e.g., [`ICRSCoords`](@ref),
[`GalCoords`](@ref), or [`FK5Coords{2000}`](@ref).

Instances are typically created with [`cartesian`](@ref). The element type of
the stored vector is carried by `TF` alone, so `cartesian(ICRSCoords{Float32}(0.1, 0.2))`,
would return a `CartesianCoords{ICRSCoords, Float32}`.

A fully parameterized tag such as `CartesianCoords{ICRSCoords{Float16}}` is also accepted.
Its element type then determines the element type of the stored vector,
so `CartesianCoords{ICRSCoords{Float16}}(c)` holds `Float16` data.
Combining a parameterized tag with a conflicting explicit `TF` throws an `ArgumentError`,
so the frame tag and the data can never disagree.
"""
struct CartesianCoords{TC <: AbstractSkyCoords, TF <: Real} <: AbstractSkyCoords
    vec::SVector{3, TF}

    function CartesianCoords{TC, TF}(vec) where {TC <: AbstractSkyCoords, TF <: Real}
        TCF = _eltype(TC)
        if !(TCF === TF || TCF === nothing)
            throw(ArgumentError("Element type $TF conflicts with the frame tag $TC. Use the element-type-free tag $(constructorof(TC)) or the matching element type $TCF"))
        end
        return new{TC, TF}(vec)
    end
end

# The element type implied by a coordinate type:
#   - `Float32` for `ICRSCoords{Float32}`
#   - `nothing` for an element-type-free frame tag such as
#     `ICRSCoords` or `FK5Coords{2000}`
_eltype(::Type{TC}) where {TC <: AbstractSkyCoords} = isconcretetype(TC) ? fieldtype(TC, 1) : nothing
_eltype(::Type{CartesianCoords{TC, TF}}) where {TC <: AbstractSkyCoords, TF <: Real} = TF

# Construction from raw components. As for the spherical types, an unspecified
# element type is chosen from the frame tag if it carries one, otherwise from
# the input (and floated); an explicit `TF` is honored exactly.
CartesianCoords{TC}(args::Real...) where {TC <: AbstractSkyCoords} = CartesianCoords{TC}(SVector(float.(args)...))
CartesianCoords{TC}(vec::AbstractVector{TF}) where {TC <: AbstractSkyCoords, TF <: Real} =
    CartesianCoords{TC, something(_eltype(TC), float(TF))}(vec)
CartesianCoords{TC, TF}(args::Real...) where {TC <: AbstractSkyCoords, TF <: Real} = CartesianCoords{TC, TF}(SVector{3, TF}(args))

# Construction from another coordinate: delegate to `convert`,
# which treats unspecified type parameters as "infer from the input"
CartesianCoords(c::AbstractSkyCoords) = convert(CartesianCoords, c)
CartesianCoords{TC}(c::AbstractSkyCoords) where {TC <: AbstractSkyCoords} = convert(CartesianCoords{TC}, c)
CartesianCoords{TC, TF}(c::AbstractSkyCoords) where {TC <: AbstractSkyCoords, TF <: Real} = convert(CartesianCoords{TC, TF}, c)
constructorof(::Type{<:CartesianCoords{TC}}) where {TC} = CartesianCoords{TC}

Base.vec(c::CartesianCoords) = c.vec

coords2cart(c::CartesianCoords) = vec(c)

"""
    cartesian(c::AbstractSkyCoords)

Returns a Cartesian representation of the coordinate `c` projected onto the unit sphere as a [`CartesianCoords`](@ref).

The result is tagged with the element-type-free frame of `c`, e.g.
`cartesian(ICRSCoords(0.1, 0.2)) isa CartesianCoords{ICRSCoords, Float64}`.

### See also
[`spherical`](@ref)
"""
cartesian(c::CartesianCoords) = c
cartesian(c::AbstractSkyCoords) = CartesianCoords{constructorof(typeof(c))}(coords2cart(c))

"""
    spherical(c::AbstractSkyCoords)

Returns the spherical representation of the coordinate `c`. Coordinate arguments that are already in spherical coordinates are simply returned. If `c` is a [`CartesianCoords`](@ref) object, the Cartesian representation is converted to spherical and returned as `TC(lon, lat)` where `c::CartesianCoords{TC}`.

### See also
[`cartesian`](@ref)
"""
spherical(c::AbstractSkyCoords) = c
spherical(c::CartesianCoords{TC}) where {TC <: AbstractSkyCoords} = TC(cart2coords(vec(c))...)

# `convert` handles target type parameters uniformly: parameters that are
# specified are honored exactly (the `convert` contract requires returning the
# requested type); unspecified ones are inferred from the input. Every method
# reduces to the two orthogonal primitives: representation change
# (`cartesian`/`spherical`) and frame rotation (`rotmat`).

# No parameters: keep the input frame, change representation only.
Base.convert(::Type{CartesianCoords}, c::AbstractSkyCoords) = cartesian(c)
Base.convert(::Type{CartesianCoords}, c::CartesianCoords) = c

# Frame tag specified: rotate. Element type inferred from the rotated vector.
Base.convert(::Type{CartesianCoords{TC}}, c::AbstractSkyCoords) where {TC <: AbstractSkyCoords} = convert(CartesianCoords{TC}, cartesian(c))
Base.convert(::Type{CartesianCoords{TC}}, c::CartesianCoords{SC}) where {TC <: AbstractSkyCoords, SC <: AbstractSkyCoords} = CartesianCoords{TC}(rotmat(TC, SC) * vec(c))
Base.convert(::Type{CartesianCoords{TC}}, c::CartesianCoords{TC}) where {TC <: AbstractSkyCoords} = c

# Frame tag and element type specified: both honored exactly.
Base.convert(::Type{CartesianCoords{TC, TF}}, c::AbstractSkyCoords) where {TC <: AbstractSkyCoords, TF <: Real} = convert(CartesianCoords{TC, TF}, cartesian(c))
Base.convert(::Type{CartesianCoords{TC, TF}}, c::CartesianCoords{SC}) where {TC <: AbstractSkyCoords, TF <: Real, SC <: AbstractSkyCoords} = CartesianCoords{TC, TF}(rotmat(TC, SC) * vec(c))
Base.convert(::Type{CartesianCoords{TC, TF}}, c::CartesianCoords{TC, TF}) where {TC <: AbstractSkyCoords, TF <: Real} = c

# Spherical target from a Cartesian source: rotate, then change representation.
function Base.convert(::Type{T}, c::CartesianCoords{SC}) where {T <: AbstractSkyCoords, SC <: AbstractSkyCoords}
    r = rotmat(T, SC) * vec(c)
    return T(cart2coords(r)...)
end

Base.:(==)(a::CartesianCoords, b::CartesianCoords) = constructorof(typeof(a)) == constructorof(typeof(b)) && vec(a) == vec(b)
Base.hash(c::CartesianCoords, h::UInt) = hash(vec(c), hash(constructorof(typeof(c)), h))
Base.isapprox(a::CartesianCoords{TC}, b::CartesianCoords{TC}; kwargs...) where {TC} =
    isapprox(vec(a), vec(b); kwargs...)
Base.isapprox(a::CartesianCoords{TCa}, b::CartesianCoords{TCb}; kwargs...) where {TCa, TCb} =
    isapprox(vec(a), vec(convert(CartesianCoords{TCa}, b)); kwargs...)
